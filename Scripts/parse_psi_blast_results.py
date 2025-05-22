import os
import time
from typing import Any, List, Dict

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
import pandas as pd

from fasta_processing import plain_to_fasta


def extract_accession(hit_id: str) -> str:
    return hit_id.rsplit("|", 2)[-2] if "|" in hit_id else hit_id

def merge_intervals(intervals: List[tuple[int, int]]) -> List[tuple[int, int]]:
    """Объединяет перекрывающиеся интервалы."""
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
    hits = []  # type: List[Dict[str, Any]]

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
        for feature in record.features:
            if "locus_tag" in feature.qualifiers.keys():
                locus_tag = feature.qualifiers["locus_tag"][0]
                gene_id_list.append(gene_id)
                locus_tag_list.append(locus_tag)

    df["gene_id"] = gene_id_list
    df["locus_tag"] = locus_tag_list

    return df


def add_org_name_gene_location(df: pd.DataFrame) -> pd.DataFrame:
    org_name_list = []
    gene_location_list = []
    gene_len_list = []

    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        gene_id = df_subset.gene_id.iloc[0]
        locus_tag = df_subset.locus_tag.iloc[0]

        time.sleep(0.333333334)
        with Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")

        genes = [f for f in record.features if f.type == "gene"]
        for gene in genes:
            if locus_tag in gene.qualifiers["locus_tag"]:
                org_name = [f for f in record.features if f.type == "source"][0].qualifiers["organism"][0]
                org_name_list.append(org_name)
                # + 1 -> Entrez 1-based
                gene_location = {
                    "start": int(gene.location.start)+1,
                    "end": int(gene.location.end),
                    "strand": 1 if gene.location.strand == 1 else 2,
                }
                gene_location_list.append(gene_location)
                gene_len_list.append(int(gene.location.end) - int(gene.location.start))
                break  # чтобы сохранялся только первый ген и не было переизбытка значений


    df["org_name"] = org_name_list
    df["gene_location"] = gene_location_list
    df["gene_len"] = gene_len_list

    return df


def add_cds_location(df: pd.DataFrame) -> pd.DataFrame:
    cds_location_list = []

    for protein_id in df.protein_id:
        df_subset = df[df["protein_id"] == protein_id]
        gene_id = df_subset.gene_id.iloc[0]
        start = df_subset.gene_location.iloc[0]["start"]
        end = df_subset.gene_location.iloc[0]["end"]
        strand = df_subset.gene_location.iloc[0]["strand"]

        time.sleep(0.333333334)
        with Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text",
                           seq_start=start, seq_stop=end, strand=strand) as handle:
            record = SeqIO.read(handle, "genbank")

        cds = [f for f in record.features if f.type == "CDS"]
        cds_location = [
            {
                "start": int(part.start)+1,
                "end": int(part.end),
                "strand": part.strand
            } for part in cds[0].location.parts
        ]
        cds_location_list.append(cds_location)

    df["cds_location"] = cds_location_list

    return df


def update_df(df: pd.DataFrame) -> pd.DataFrame:
    df = add_gene_id_locus_tag(df)
    df = add_org_name_gene_location(df)
    df = add_cds_location(df)
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
                exons_dict[f">{i}:{start}:{end}"] = exon_seq

        with open(f"{dir}/{protein_id}/cds_plain.fna", "w") as outfile:
            outfile.write(f">{org_name}\n{cds}\n")
        plain_to_fasta(f"{dir}/{protein_id}/cds_plain.fna")

        with open(f"{dir}/{protein_id}/exons.fa", "w") as outfile:
            for header, exon_seq in exons_dict.items():
                outfile.write(f">{header}\n{exon_seq}\n")


def save_proteins(df, dir: str = "../Sequences_protein_id") -> None:
    for protein_id in df.protein_id:
        os.makedirs(f"{dir}/{protein_id}", exist_ok=True)

        time.sleep(0.333333334)
        with Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text") as handle:
            with open(f"{dir}/{protein_id}/protein.faa", "w") as outfile:
                outfile.write(handle.read())
