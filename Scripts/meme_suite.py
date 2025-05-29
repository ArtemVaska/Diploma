import os
import subprocess
import re

import pandas as pd

from fasta_processing import read_fasta, read_single_fasta
from tg_logger import telegram_logger


def create_input_meme(df: pd.DataFrame, phylum: str) -> None:
    introns_dict = {}

    for sub_phylum in df.sub_phylum.unique():
        introns = read_fasta(f"../rnafold/{sub_phylum}/{sub_phylum}_introns.fa")
        introns_dict.update(introns)  # создаем словарь для всех интронов у лучеперых рыб

    os.makedirs(f"../meme/{phylum}", exist_ok=True)
    with open(f"../meme/{phylum}/{phylum}_introns.fa", "w") as outfile:
        for header, seq in introns_dict.items():
            outfile.write(f">{header}\n{seq}\n")  # сохраняем его в папку для meme инпута в отдельную папку


@telegram_logger(chat_id=611478740)
def run_meme(
        input_file: str, output_dir: str,
        nmotifs: int = 5, minw: int = 10, maxw: int = 60,
) -> None:

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    os.makedirs(output_dir, exist_ok=True)

    meme_path = "/home/artemvaska/meme/bin/meme"  # ← полный путь, т.к. мы запускаем в мамбе
    cmd = [
        meme_path, input_file,
        "-dna",
        "-oc", f"../meme/{output_dir}",
        "-mod", "zoops",
        "-nmotifs", str(nmotifs),
        "-minw", str(minw),
        "-maxw", str(maxw),
        "-objfun", "classic",
        "-markov_order", "0"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MEME failed:\n{result.stderr}")
    else:
        print(f"Meme finished successfully!\n"
              f"Input: {input_file}\n"
              f"Output: {output_dir}\n")


def parse_meme_block(file_path: str, motif: str) -> pd.DataFrame | None:
    with open(file_path) as f:
        content = f.read()

    # Найдём весь нужный блок по мотиву
    pattern = (
        r"-{80,}\n"                             # Верхняя линия
        r"\s*Motif " + re.escape(motif) + r" MEME-\d+ in BLOCKS format\n"  # Название мотива
        r"-{80,}\n"                             # Нижняя линия
        r"BL\s+MOTIF\s+" + re.escape(motif) +
        r"\s+width=(\d+)\s+seqs=(\d+)\n"        # Строка с BL MOTIF
        r"(.*?)(?://)"                         # Все строки до // (ленивый захват)
    )

    match = re.search(pattern, content, re.DOTALL)
    if not match:
        print("Motif not found.")
        return None

    width, seqs, block = match.groups()
    print(f"MOTIF: {motif}, width={width}, seqs={seqs}")

    # Парсинг таблицы
    rows = []
    for line in block.strip().split("\n"):
        m = re.match(r"^(\S+)\s+\(\s*(\d+)\)\s+([A-Z]+)", line)
        if m:
            org, length, site_seq = m.groups()
            rows.append({
                "org_name": org,
                "intron_length": int(length),
                "site_seq": site_seq
            })
    df = pd.DataFrame(rows)

    return df


@telegram_logger(chat_id=611478740)
def run_tomtom(
        input_file: str, output_dir: str,
        e_thresh: float = 10.0, mi: int = None
) -> None:

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    os.makedirs(output_dir, exist_ok=True)

    tomtom_path = "/home/artemvaska/meme/bin/tomtom"  # ← полный путь, т.к. мы запускаем в мамбе
    cmd = [
        tomtom_path,
        "-no-ssc",
        "-oc", output_dir,
        "-min-overlap", "5",
        "-dist", "pearson",
        "-evalue", "-thresh", str(e_thresh),
        input_file,  # результат из meme, который лежит в папке с таким же названием по филе
        "/home/artemvaska/meme-5.5.8/dbs/JASPAR2024_ALL.meme"  # ← полный путь, название БД
    ]

    # с помощью mi можно указать номер мотива, который мы хотим сравнить с базой данных
    if mi is not None:
        cmd.extend(["-mi", str(mi)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Tomtom failed:\n{result.stderr}")
    else:
        print(f"Tomtom finished successfully!\n"
              f"Input: {input_file}\n"
              f"Output: {output_dir}\n")


def add_highlight_coords(df_info: pd.DataFrame, df_to_upd: pd.DataFrame) -> pd.DataFrame:
    coords_list = []
    # нужно обрезать значение в этом столбце, т.к. meme не может вывести больше)) FIXME
    df_info_org_name_24_len = df_info.copy()
    df_info_org_name_24_len["org_name"] = df_info_org_name_24_len["org_name"].str[:24]
    # в изначальной таблице выбираем только те виды, для которых нашлись мотивы
    df_subset = df_info_org_name_24_len[df_info_org_name_24_len.org_name.isin(df_to_upd.org_name)]

    for row in df_subset.iterrows():
        sub_phylum = row[1].sub_phylum  # достаем названия сабфилы == название папок для вторичных структур
        org_name = row[1].org_name
        org_name_protein_id = row[1].org_name_protein_id  # название папки с организмом
        source = row[1].source
        match source:
            case "datasets":
                cds_cassette_seq = read_single_fasta(f"../Datasets/{sub_phylum}/{org_name_protein_id}/ncbi_dataset/data/cds_cassette_plain.fa")
            case "psi_blast":
                cds_cassette_seq = read_single_fasta(f"../Sequences_protein_id/{sub_phylum}/{org_name_protein_id.split("__")[1]}/cds_cassette_plain.fa")

        site_seq = df_to_upd[df_to_upd["org_name"] == org_name].site_seq.values[0]
        start = cds_cassette_seq.find(site_seq)
        end = cds_cassette_seq.find(site_seq) + len(site_seq)
        coords_list.append((start, end))

    df_to_upd["org_name"] = df_info.org_name  # фиксим имена организмов
    df_to_upd["highlight_coords"] = coords_list

    return df_to_upd


def highlight_coords(df_info: pd.DataFrame, df_meme: pd.DataFrame) -> None:
    df_subset = df_info[df_info.org_name.isin(df_meme.org_name)]

    for row in df_subset.iterrows():
        sub_phylum = row[1].sub_phylum
        org_name = row[1].org_name
        start, end = df_meme[df_meme["org_name"] == org_name].highlight_coords.values[0]

        organism_dir = Path(f"../rnafold/{sub_phylum}/{sub_phylum}_structure/{org_name}").resolve()

        for fold_type in ["mfe", "centroid"]:  # выделяем для обоих вариантов структур

            # достаем координаты для выделения интрона из уже существующего файла
            with open(f"{organism_dir}/{fold_type}_highlight.eps") as structure_txt:
                lines = structure_txt.readlines()
                intron_highlight_cmd = lines[lines.index("% Start Annotations\n")+1].rstrip()

            highlight_cmd = ["--pre", f"{intron_highlight_cmd} {start} {end} 8 0 10 0 omark"]  # выделяем желтым интрон и зеленым мотив по координатам
            with open(f"{organism_dir}/{fold_type}.fold") as infile:
                subprocess.run(
                    ["RNAplot", "--jobs=16"] + highlight_cmd,
                    stdin=infile,
                    cwd=organism_dir,
                )
            motif_highlight_eps = organism_dir / "rna.eps"
            if motif_highlight_eps.exists():
                motif_highlight_eps.rename(organism_dir / f"{fold_type}_motif_highlight.eps")
            else:
                print(f"⚠️ Highlight file for {sub_phylum}/{org_name} not created")
