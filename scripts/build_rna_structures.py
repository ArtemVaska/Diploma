import os
import subprocess
from pathlib import Path
import argparse
import sys

import pandas as pd

from tg_logger import telegram_logger
from fasta_processing import read_single_fasta


def read_fasta(file_path: Path) -> dict:
    sequences = {}
    current_id = None
    with open(file_path) as infile:
        for line in infile:
            if line.startswith(">"):
                current_id = line[1:].strip()
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.strip()
    return sequences


def create_input_files(df: pd.DataFrame, sub_phylum: str) -> None:

    dir = "../rnafold"
    file_dict = {}

    df_sub_phylum = df[df["sub_phylum"] == sub_phylum].copy()
    for row in df_sub_phylum.iterrows():
        if row[1].source == "datasets":
            prefix = "../datasets"
            postfix = "ncbi_dataset/data/"  # здесь обязательно еще один слеш прописываем
            intron_name = "cassette"  # FIXME
            org_name = row[1].org_name_protein_id
            new_org_name = org_name.rsplit("_", 1)[0].capitalize()
        elif row[1].source == "psi_blast":
            prefix = "../sequences_protein_id"
            postfix = ""  # а тут нет, потому что да
            intron_name = "cassette_intron"  # FIXME
            org_name = row[1].org_name_protein_id.split("__")[1]
            new_org_name = row[1].org_name_protein_id.split("__")[0]
        else:
            raise ValueError("Unknown source")

        input_file = f"{prefix}/{sub_phylum}/{org_name}/{postfix}cds_cassette.fa"  # тут не ставим доп. слеш
        cds_cassette_seq = read_single_fasta(input_file)
        intron_seq = read_fasta(f"{prefix}/{sub_phylum}/{org_name}/{postfix}cassette.fa")[intron_name]  # и тут

        file_dict[new_org_name] = [cds_cassette_seq, intron_seq]

    os.makedirs(f"{dir}/{sub_phylum}", exist_ok=True)

    with open(f"{dir}/{sub_phylum}/{sub_phylum}_cds_cassette.fa", "w") as cds_cassette_outfile:
        for org_name, (cds_cassette_seq, _) in file_dict.items():
            cds_cassette_outfile.write(f">{org_name}\n{cds_cassette_seq}\n")

    with open(f"{dir}/{sub_phylum}/{sub_phylum}_introns.fa", "w") as introns_outfile:
        for org_name, (_, intron) in file_dict.items():
            introns_outfile.write(f">{org_name}\n{intron}\n")


def find_intron_coordinates(full_seq: str, intron_seq: str) -> tuple[int, int] | None:
    start = full_seq.find(intron_seq)
    if start == -1:
        return None
    return (start + 1, start + len(intron_seq))  # 1-based inclusive


# @telegram_logger(chat_id=611478740)
def run_rnafold_with_highlight(
    fasta_path: Path,
    paint_path: Path | None = None,
    output_dir: Path | None = None
) -> None:
    if not fasta_path.exists():
        print(f"❌ Error: input file {fasta_path} not found.")
        sys.exit(1)

    input_sequences = read_fasta(fasta_path)
    intron_sequences = read_fasta(paint_path) if paint_path else {}

    # Determine output root directory
    if output_dir:
        output_root = output_dir.resolve()
    else:
        project_root = Path(__file__).resolve().parents[1]
        output_root = project_root
    output_root.mkdir(parents=True, exist_ok=True)

    for seq_id, full_seq in input_sequences.items():
        print(f"→ Processing {seq_id}...")
        organism_dir = output_root / seq_id
        organism_dir.mkdir(parents=True, exist_ok=True)

        # проверка на то, что организм уже обработан
        fold_file = organism_dir / "mfe.fold"
        if fold_file.exists():
            print(f"Files for {seq_id} already exist\n")
            continue

        # Run RNAfold (creates only dot.ps)
        result = subprocess.run(
            ["RNAfold", "--jobs=16", "-p", "-d2", "--noLP", "--noPS"],
            input=full_seq.strip(),
            cwd=organism_dir,
            capture_output=True,
            text=True
        )
        lines = result.stdout.strip().splitlines()
        if len(lines) < 5:
            print(f"⚠️ Error: RNAfold failed for {seq_id}")
            continue

        # Save result to .fold file
        fold_file = organism_dir / f"all.fold"
        with open(fold_file, "w") as f:
            for line in lines:
                f.write(f"{line.strip()}\n") # delete space in line-5

        # Create .fold files and plots for mfe and centroid
        for fold_type in ["mfe", "centroid"]:
            match fold_type:
                case "mfe":
                    fold_line = lines[1]
                case "centroid":
                    fold_line = lines[3]
            with open(f"{organism_dir}/{fold_type}.fold", "w") as f:
                f.write(f"{full_seq}\n{fold_line}")

            with open(f"{organism_dir}/{fold_type}.fold") as infile:
                subprocess.run(
                    ["RNAplot", "--jobs=16"],
                    stdin=infile,
                    cwd=organism_dir,
                )
            rna_eps = organism_dir / "rna.eps"
            if rna_eps.exists():
                rna_eps.rename(organism_dir / f"{fold_type}_rna.eps")

            # Run RNAplot with intron highlighting if available
            if seq_id in intron_sequences:
                coords = find_intron_coordinates(full_seq, intron_sequences[seq_id])
                if coords:
                    start, end = coords
                    highlight_cmd = ["--pre", f"{start} {end} 8 10 10 0 omark"]
                    with open(f"{organism_dir}/{fold_type}.fold") as infile:
                        subprocess.run(
                            ["RNAplot", "--jobs=16"] + highlight_cmd,
                            stdin=infile,
                            cwd=organism_dir,
                        )
                    highlight_eps = organism_dir / "rna.eps"
                    if highlight_eps.exists():
                        highlight_eps.rename(organism_dir / f"{fold_type}_highlight.eps")
                    else:
                        print(f"⚠️ Highlight file not created for {seq_id}")
                else:
                    print(f"⚠️ Warning: intron not found in {seq_id}. Skipping highlight")

        print(f"✔ Done: {seq_id}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RNAfold with intron highlighting")
    parser.add_argument("--input", type=str, required=True, help="FASTA with full exon-intron sequences")
    parser.add_argument("--paint", type=str, required=False, help="FASTA with intron sequences to highlight")
    parser.add_argument("--output", type=str, required=False, help="Directory to store output files")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_path = Path.cwd() / args.input
    paint_path = Path.cwd() / args.paint if args.paint else None
    output_path = Path.cwd() / args.output if args.output else None
    run_rnafold_with_highlight(input_path, paint_path, output_path)
