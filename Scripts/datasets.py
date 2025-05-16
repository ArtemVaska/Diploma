import os
import time
import urllib.error
from subprocess import check_call, DEVNULL, STDOUT, CalledProcessError

import pandas as pd
from Bio import Entrez, SeqIO

from fasta_processing import read_single_fasta
from tg_logger import telegram_logger

Entrez.email = "artemvaskaa@gmail.com"


@telegram_logger(chat_id=611478740)
def select_all_phylas(df: pd.DataFrame) -> (dict, list):
    """
    Fetches taxonomy info for every taxon ID in the table.

    :param df: pd.DataFrame
    :return: dict, list
    """
    phylas = {}

    taxids_error = []
    for index, row in df.iterrows():
        try:
            stream = Entrez.efetch(db="Taxonomy", id=str(index), retmode="xml")
        except urllib.error.HTTPError:
            print(f"ERROR taxid: {index}, {row['org_name']}")
            taxids_error.append(index)
            continue
        records = Entrez.read(stream)
        phylas[index] = f"{records[0]['Lineage']}; {records[0]['ScientificName']}"
        time.sleep(0.333333334)
    return phylas, taxids_error


@telegram_logger(chat_id=611478740)
def select_phyla(df: pd.DataFrame, phylas: list, taxids: list = None) -> (dict, list):
    """
    Selects all taxon_id from the table that match the given phylas.

    :param df:
    :param phylas:
    :return:
    """
    phyla_taxids = {}
    for phyla in phylas:
        phyla_taxids[phyla] = []

    if taxids is not None:
        df = df.loc[taxids]

    taxids_error = []
    for index, row in df.iterrows():
        try:
            stream = Entrez.efetch(db="Taxonomy", id=str(index), retmode="xml")
        except urllib.error.HTTPError:
            print(f"ERROR taxid: {index}, {row['org_name']}")
            taxids_error.append(index)
            continue
        records = Entrez.read(stream)
        for phyla in phylas:
            if phyla in records[0]["Lineage"]:
                phyla_taxids[phyla].append(index)
                print(f"{phyla}: {index}, {row['org_name']}")
        time.sleep(0.333333334)
    return phyla_taxids, taxids_error


def download_subset_df_datasets(df: pd.DataFrame, phyla: str = "") -> list:
    """
    Downloads gene, rna and protein for every GeneID in the provided dataframe
    via ncbi-datasets-cli.
    Additional function for download_all_files_ncbi.

    :param df: Pandas DataFrame
    :param phyla: Phylum
    :return: list of org_names for next analysis
    """
    if phyla != "":
        phyla = f"{phyla}/"
        if not os.path.exists(f"../Datasets/{phyla}"):
            os.mkdir(f"../Datasets/{phyla}")

    org_names = []
    for i, (index, row) in enumerate(df.iterrows()):
        gene_id = str(row["gene_id"])
        org_name = row["org_name"].lower().replace(" ", "_").replace("/", "_").replace("-", "_") + f"_{i}"
        shell_commands = [
            ["datasets", "download", "gene", "gene-id", gene_id,
             "--include", "gene,cds,rna,protein",
             "--filename", f"../Datasets/{phyla}{org_name}.zip"],
            ["unzip", f"../Datasets/{phyla}{org_name}.zip",
             "-d", f"../Datasets/{phyla}{org_name}"],
            ["rm", "-r", f"../Datasets/{phyla}{org_name}.zip"]
        ]
        for command in shell_commands:
            try:
                check_call(command, stdout=DEVNULL, stderr=STDOUT)
            except CalledProcessError:
                print(f"ERROR command: {command}")
                continue
        print(f"Gene, mRNA, protein for {phyla}{org_name} downloaded successfully")
        org_names.append(org_name)
    return org_names


def download_gene_gb(phyla: str, org_names: list) -> None:
    """
    Downloads gene in GenBank format for provided organisms.
    Additional function for download_all_files_ncbi.

    :param phyla:
    :param org_names:
    :return:
    """
    # Obtain mRNA Accession and ranges of the gene
    for org_name in org_names:
        try:
            with open(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/gene.fna") as infile:
                line = infile.readline().rstrip()
                gene_acc = line.split(":")[0][1:]
                gene_ranges = line.split()[0].split(":")[1]
                gene_range_1 = gene_ranges.split("-")[0]
                gene_range_2 = gene_ranges.split("-")[1]
        except FileNotFoundError:
            continue
        if "c" in gene_range_1:
            gene_range_1 = gene_range_1.replace("c", "")
            strand = 2
        else:
            strand = 1
        gene_range_1, gene_range_2 = int(gene_range_1), int(gene_range_2)

        # Obtain and download mRNA in GenBank format
        stream = Entrez.efetch(db="nucleotide", id=gene_acc, idtype="acc",
                               seq_start=gene_range_1, seq_stop=gene_range_2, strand=strand,
                               rettype="gb", retmode="text")
        with open(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/gene.gb", "w") as outfile:
            outfile.write(stream.read())

        print(f"Gene GenBank for {phyla}/{org_name} downloaded successfully")

        time.sleep(0.333333334)


def parse_exon_ranges(phyla: str, org_names: list, feature_type: str = "mRNA") -> dict:
    """
    Parses exon ranges for provided organisms to select exons for provided organisms.
    Additional function for download_all_files_ncbi.

    :param phyla:
    :param org_names:
    :param feature_type:
    :return:
    """
    exons_dict = {}

    feature_types = ["mRNA", "CDS"]
    if feature_type not in feature_types:
        raise Exception(f"Feature type {feature_type} is not supported")

    for org_name in org_names:
        try:
            records = SeqIO.parse(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/gene.gb", "genbank")
        except FileNotFoundError:
            continue
        for record in records:
            mrnas = [feature for feature in record.features if feature.type == "mRNA"]
            parts = {
                i: [int(part.start), int(part.end)] for i, part in enumerate(mrnas[0].location.parts)
            }
        exons_dict[org_name] = parts

    return exons_dict


def create_exons(phyla: str, exon_ranges: dict) -> None:
    """
    Creates exons for provided organisms by provided exon ranges.
    Additional function for download_all_files_ncbi.

    :param phyla:
    :param exon_ranges:
    :return:
    """
    for org_name in exon_ranges.keys():
        gene_seq = read_single_fasta(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/gene.fna")
        with open(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/exons.fa", "w") as outfile:
            for exon_i, exon_range in exon_ranges[org_name].items():
                outfile.write(f">{exon_i}:{exon_range[0]}-{exon_range[1]}\n"
                              f"{gene_seq[exon_range[0]:exon_range[1]]}\n")
        print(f"Exons for {phyla}/{org_name} created successfully")


@telegram_logger(chat_id=611478740)
def download_all_files_ncbi(df: pd.DataFrame,
                            phyla_taxids: dict,
                            phylas: list,
                            feature_type: str = "mRNA") -> None:
    """
    Downloads gene.fna, rna.fna, protein.faa and gene.gb
    for every GeneID in the provided dataframe of selected phylas.

    :param df:
    :param phyla_taxids:
    :param phylas:
    :param feature_type:
    :return:
    """
    feature_types = ["mRNA"]
    if feature_type not in feature_types:
        raise Exception(f"Feature type {feature_type} is not supported")

    for phyla in phylas:
        if os.path.exists(f"../Datasets/{phyla}"):
            print(f"Files for {phyla} already downloaded")
            continue

        df_subset = df.loc[phyla_taxids[phyla]]
        org_names = download_subset_df_datasets(df_subset, phyla)
        download_gene_gb(phyla, org_names)
        exon_ranges = parse_exon_ranges(phyla, org_names, feature_type)
        create_exons(phyla, exon_ranges)
        print()


def check_transcript_count(phylas: list):
    message = False
    for phyla in phylas:
        for org_name in os.listdir(f"../Datasets/{phyla}"):
            with open(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/cds.fna") as infile:
                lines = "".join(infile.readlines())
                transcript_count = lines.count(">")
                if transcript_count > 1:
                    message = True
                    print(f"{phyla}/{org_name}: {transcript_count} transcripts")
    if message:
        print(f"\nDelete other transcripts from cds.fna, protein.faa and rna.fna")
