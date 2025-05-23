import os
import time
from typing import Any, List, Dict

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqFeature import SeqFeature

import pandas as pd

from data_processing import find_codon
from fasta_processing import plain_to_fasta, read_single_fasta, read_fasta
from tg_logger import telegram_logger


def extract_accession(hit_id: str) -> str:
    return hit_id.rsplit("|", 2)[-2] if "|" in hit_id else hit_id


def merge_intervals(intervals: List[tuple[int, int]]) -> List[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:  # перекрытие
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged


def parse_psiblast_xml(xml_path: str, max_len: int = 5000) -> pd.DataFrame:
    with open(xml_path) as f:
        try:
            record = next(NCBIXML.parse(f))  # Берём первую итерацию
        except StopIteration:
            return pd.DataFrame()  # Если файл пуст

    query_len = record.query_length
    hits = []

    for alignment in record.alignments:
        hsps = [hsp for hsp in alignment.hsps if hsp.align_length <= max_len]
        if not hsps:
            continue

        # Сбор покрытых интервалов на query
        query_intervals = [(min(hsp.query_start, hsp.query_end),
                            max(hsp.query_start, hsp.query_end)) for hsp in hsps]
        merged_query_intervals = merge_intervals(query_intervals)
        covered_len = sum(end - start + 1 for start, end in merged_query_intervals)
        qc = round(covered_len / query_len, 2)
        qc = min(qc, 1.0)  # Защита от >1

        # Координаты на subject
        sbjct_starts = [min(hsp.sbjct_start, hsp.sbjct_end) for hsp in hsps]
        sbjct_ends = [max(hsp.sbjct_start, hsp.sbjct_end) for hsp in hsps]
        sbjct_start = min(sbjct_starts)
        sbjct_end = max(sbjct_ends)

        # Средний процент идентичности
        identities = sum(hsp.identities for hsp in hsps)
        alignment_len = sum(hsp.align_length for hsp in hsps)
        per_ident = round(identities / alignment_len, 2) if alignment_len else 0.0

        hits.append({
            "hit_id": alignment.hit_id,
            "protein_id": extract_accession(alignment.hit_id),
            "qc": qc,
            "per_ident": per_ident,
            "sbjct_start": sbjct_start,
            "sbjct_end": sbjct_end,
            "sbjct_len": sbjct_end - sbjct_start + 1,
            "query_len": query_len,
        })

    df = pd.DataFrame(hits)
    if not df.empty:
        df.set_index("hit_id", inplace=True)
    return df


def filter_psiblast_hits(
        df: pd.DataFrame,
        min_qc: float = 0.8,
        min_ident: float = 0.8,
        min_sbjct_len: int = 500
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()

    filtered_df = df[
        (df["qc"] >= min_qc) &
        (df["per_ident"] >= min_ident) &
        (df["sbjct_len"] >= min_sbjct_len)
        ]

    return filtered_df


def add_gene_id_locus_tag(df: pd.DataFrame) -> pd.DataFrame:
    gene_id_list = []
    locus_tag_list = []

    for protein_id in df.protein_id:
        time.sleep(0.333333334)
        with Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
        gene_id = record.annotations["db_source"].split()[-1]

        if gene_id.startswith("XM"):
            cds_feature = [f for f in record.features if f.type == "CDS"][0]
            try:
                locus_tag = cds_feature.qualifiers["locus_tag"][0]
            except KeyError:
                locus_tag = cds_feature.qualifiers["gene"][0]
        else:
            for feature in record.features:
                if "locus_tag" in feature.qualifiers.keys():
                    locus_tag = feature.qualifiers["locus_tag"][0]
        gene_id_list.append(gene_id)
        locus_tag_list.append(locus_tag)

    df["gene_id"] = gene_id_list
    df["locus_tag"] = locus_tag_list

    return df


def parse_location(gene: SeqFeature, cds: SeqFeature):
    gene_location = {
        "start": int(gene.location.start) + 1,
        "end": int(gene.location.end),
        "strand": 1 if gene.location.strand == 1 else 2,
    }
    gene_len = int(gene.location.end) - int(gene.location.start)

    cds_location = [
        {
            "start": int(part.start) + 1,
            "end": int(part.end),
            "strand": 1 if part.strand == 1 else 2,
        } for part in cds.location.parts
    ]
    return gene_location, gene_len, cds_location


def obtain_gene_cds_location_locus_tag(record: SeqIO.SeqRecord, locus_tag: str):
    features_with_locus_tag = [f for f in record.features if "locus_tag" in f.qualifiers]
    features_locus_tag = [f for f in features_with_locus_tag if locus_tag in f.qualifiers["locus_tag"]]

    gene = [f for f in features_locus_tag if f.type == "gene"][0]
    cds = [f for f in features_locus_tag if f.type == "CDS"][0]

    gene_location, gene_len, cds_location = parse_location(gene, cds)

    return gene_location, gene_len, cds_location


def obtain_gene_cds_location_xm(record: SeqIO.SeqRecord, locus_tag: str):
    # по XM достаем Gene Accession, из которого и будем доставать координаты
    time.sleep(0.333333334)
    with Entrez.efetch(db="gene", id=locus_tag, retmode="xml") as handle:
        record = Entrez.read(handle)

    gene_acc = record[0]["Entrezgene_locus"][0]["Gene-commentary_accession"]
    seq_interval = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
    start = int(seq_interval["Seq-interval_from"]) + 1
    end = int(seq_interval["Seq-interval_to"]) + 1
    strand = 1 if seq_interval["Seq-interval_strand"]["Na-strand"].attributes["value"] == "plus" else 2

    # по найденному Accession достаем координаты для gene и cds
    time.sleep(0.333333334)
    with Entrez.efetch(db="nucleotide", id=gene_acc, idtype="acc",
                       seq_start=start, seq_stop=end, strand=strand,
                       rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
    cds = [f for f in record.features if f.type == "CDS"][0]

    gene_location = {"start": start, "end": end, "strand": strand}
    gene_len = end - start + 1

    cds_location = [
        {
            "start": int(part.start) + 1,
            "end": int(part.end),
            "strand": 1 if part.strand == 1 else 2,
        } for part in cds.location.parts
    ]

    return gene_acc, gene_location, gene_len, cds_location


def add_org_name_gene_cds_location(df: pd.DataFrame) -> pd.DataFrame:
    org_name_list = []
    gene_location_list = []
    gene_len_list = []
    gene_id_list = []
    cds_location_list = []
    protein_id_delete_list = []

    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        gene_id = df_subset.gene_id.iloc[0]
        locus_tag = df_subset.locus_tag.iloc[0]

        time.sleep(0.333333334)
        with Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")

        if gene_id.startswith("XM"):
            # обновляем gene_id, т.к. координаты для XM будут доставать через Gene Accession
            try:
                gene_id, gene_location, gene_len, cds_location = obtain_gene_cds_location_xm(record, locus_tag)
                # добавляем org_name только если получилось извлечь координаты
                org_name = [f for f in record.features if f.type == "source"][0].qualifiers["organism"][0].replace(
                    " ", "_")
                org_name_list.append(org_name)
            except KeyError:
                print(f"KeyError: ProteinID {protein_id} GeneID {gene_id} -> skipping...")
                protein_id_delete_list.append(protein_id)
                continue
        else:
            try:
                gene_location, gene_len, cds_location = obtain_gene_cds_location_locus_tag(record, locus_tag)
                # добавляем org_name только если получилось извлечь координаты
                org_name = [f for f in record.features if f.type == "source"][0].qualifiers["organism"][0].replace(
                    " ", "_")
                org_name_list.append(org_name)
            except IndexError:
                print(f"IndexError: ProteinID {protein_id} GeneID not found -> skipping...")
                protein_id_delete_list.append(protein_id)
                continue

        gene_id_list.append(gene_id)
        gene_location_list.append(gene_location)
        gene_len_list.append(gene_len)
        cds_location_list.append(cds_location)

    df = df[~df["protein_id"].isin(protein_id_delete_list)].copy()
    df["gene_id"] = gene_id_list
    df["org_name"] = org_name_list
    df["gene_location"] = gene_location_list
    df["gene_len"] = gene_len_list
    df["cds_location"] = cds_location_list

    return df


@telegram_logger(chat_id=611478740)
def update_df(df: pd.DataFrame) -> pd.DataFrame:
    df = add_gene_id_locus_tag(df)
    df = add_org_name_gene_cds_location(df)
    df["protein_id"] = df["protein_id"].str.replace("_", "")  # если XP в начале -> загрузка файлов
    return df


def save_genes(df, dir: str = "../Sequences_protein_id") -> None:
    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        gene_id = df_subset.gene_id.iloc[0]
        gene_location = df_subset.gene_location.iloc[0]
        start, end, strand = gene_location["start"], gene_location["end"], gene_location["strand"]

        os.makedirs(f"{dir}/{protein_id}", exist_ok=True)

        time.sleep(0.333333334)
        with Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text",
                           seq_start=start, seq_stop=end, strand=strand) as handle:
            with open(f"{dir}/{protein_id}/gene.fna", "w") as outfile:
                outfile.write(handle.read())


def save_cdss_exons(df, dir: str = "../Sequences_protein_id") -> None:
    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        org_name = df_subset.org_name.iloc[0]
        gene_id = df_subset.gene_id.iloc[0]
        cds_location = df_subset.cds_location.iloc[0]

        os.makedirs(f"{dir}/{protein_id}", exist_ok=True)

        cds = ""
        exons_dict = {}
        for i, exon_location in enumerate(cds_location):
            start, end, strand = exon_location["start"], exon_location["end"], exon_location["strand"]

            time.sleep(0.333333334)
            with Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text",
                               seq_start=start, seq_stop=end, strand=strand) as handle:
                record = SeqIO.read(handle, "fasta")

                exon_seq = str(record.seq)
                cds += exon_seq
                exons_dict[f">{i}:{start}-{end}"] = exon_seq

        with open(f"{dir}/{protein_id}/cds_plain.fna", "w") as outfile:
            outfile.write(f">{org_name}\n{cds}\n")
        plain_to_fasta(f"{dir}/{protein_id}/cds_plain.fna")

        with open(f"{dir}/{protein_id}/exons.fa", "w") as outfile:
            for header, exon_seq in exons_dict.items():
                outfile.write(f"{header}\n{exon_seq}\n")


def save_proteins(df, dir: str = "../Sequences_protein_id") -> None:
    for protein_id in df.protein_id:
        os.makedirs(f"{dir}/{protein_id}", exist_ok=True)

        time.sleep(0.333333334)
        with Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text") as handle:
            with open(f"{dir}/{protein_id}/protein.faa", "w") as outfile:
                outfile.write(handle.read())


# @telegram_logger(chat_id=611478740)
def save_files(df: pd.DataFrame, dir: str = "../Sequences_protein_id") -> None:
    save_genes(df, dir)
    save_cdss_exons(df, dir)
    save_proteins(df, dir)


def create_cassette(path: str, df_exons: pd.DataFrame, exons_i: list) -> dict:
    gene = read_single_fasta(f"{path}/gene.fna")

    exon_110 = df_exons.loc[exons_i[0]]
    exon_37 = df_exons.loc[exons_i[1]]

    cassette_start = gene.find(exon_110["sequence"]) + len(exon_110["sequence"])
    cassette_end = gene.find(exon_37["sequence"])
    cassette_intron = gene[cassette_start:cassette_end]

    cassette_dict = {
        f"{exons_i[0]}:{exon_110.coords}": exon_110.sequence,
        "cassette_intron": cassette_intron,
        f"{exons_i[1]}:{exon_37.coords}": exon_37.sequence,
    }
    with open(f"{path}/cassette.fa", "w") as outfile:
        for header, seq in cassette_dict.items():
            outfile.write(f">{header}\n"
                          f"{seq}\n")

    with open(f"{path}/cds.fna") as infile:
        line = infile.readline()
        cds_header = line.rstrip()
    cds = read_single_fasta(f"{path}/cds.fna")
    cds_cassette = cds.replace(f"{exon_110.sequence}{exon_37.sequence}",
                               f"{exon_110.sequence}{cassette_intron}{exon_37.sequence}")
    with open(f"{path}/cds_cassette_plain.fa", "w") as outfile:
        outfile.write(f"{cds_header} [+cassette_intron]\n")
        outfile.write(f"{cds_cassette}\n")
    plain_to_fasta(f"{path}/cds_cassette_plain.fa", 70)

    return cassette_dict


def concat_cassette(cassette_dict: dict, concat_type: str) -> str | None:
    match concat_type:
        case "i":
            return "".join([seq for header, seq in cassette_dict.items() if header == 'cassette_intron'])
        case "ee":
            return "".join([seq for header, seq in cassette_dict.items() if header != 'cassette_intron'])
        case "eie":
            return "".join([seq for header, seq in cassette_dict.items()])
    return None


def create_many_cassettes(dir: str, data: dict) -> dict:
    introns = {}
    for protein_id_org_name, (df, exons_i) in data.items():
        protein_id = protein_id_org_name.split("_")[-1]
        cassette = create_cassette(f"{dir}/{protein_id}", df, exons_i=exons_i)
        introns[protein_id_org_name] = concat_cassette(cassette, "i")
    return introns


def dict_align_create(df: pd.DataFrame, align_type: str, dir: str = "../Sequences_protein_id") -> dict:
    align_types = ["gene", "rna", "protein", "cds", "cassette", "cds_cassette"]
    if align_type not in align_types:
        raise ValueError(f"Unknown alignment type: {align_type}")
    else:
        match align_type:
            case "gene":
                ext = "fna"
            case "rna":
                ext = "fna"
            case "protein":
                ext = "faa"
            case "cds":
                ext = "fna"
            case "cassette":
                ext = "fa"
            case "cds_cassette":
                ext = "fa"

    filename = f"{align_type}.{ext}"
    dict_align = {}

    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        org_name = df_subset.org_name.iloc[0]
        dict_align[f"{org_name}_{protein_id}"] = read_single_fasta(f"{dir}/{protein_id}/{filename}")

    return dict_align


def dict_align_info_analyze(df: pd.DataFrame, feature: str, dir: str = "../Sequences_protein_id") -> (pd.DataFrame,
                                                                                                      dict):
    rows = []
    dict_align = dict_align_create(df, feature, dir)

    for org_name_protein_id, cds_seq in dict_align.items():
        protein_id = org_name_protein_id.split("_")[-1]
        stop_codon_pos = find_codon(cds_seq, which="stop")
        cassette_intron = read_fasta(f"{dir}/{protein_id}/cassette.fa")["cassette_intron"]
        cassette_intron_start = cds_seq.find(cassette_intron)
        rows.append(
            {
                "org_name_protein_id": org_name_protein_id,
                "stop_codon_pos": stop_codon_pos,
                "equal_to_cds": stop_codon_pos + 3 == len(cds_seq),
                "cassette_intron_start": cassette_intron_start,
                "intron_length_to_stop_codon": stop_codon_pos - cassette_intron_start,
                "intron_length": len(cassette_intron)
            }
        )
    df = pd.DataFrame(rows)

    return df, dict_align
