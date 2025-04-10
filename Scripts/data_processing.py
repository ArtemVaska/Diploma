import os
import time
import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from Bio import Entrez

Entrez.email = "artemvaskaa@gmail.com"


def cluster_analysis_preview(df: pd.DataFrame) -> None:
    """
    Prints preview cluster analysis (DBSCAN) based on QC

    :param df: dataframe with parameters for analysis
    """
    qcs_values = df.QC.values.tolist()

    for eps in range(1, 11):
        df_copy = df.copy()
        eps = eps / 100  # to get 0.01 etc.
        dbscan = DBSCAN(eps=eps, min_samples=1)
        clustering = dbscan.fit(np.asarray(qcs_values).reshape(-1, 1))
        labels = clustering.labels_

        n_clusters = len(np.unique(labels))
        df_copy["Cluster"] = labels

        qcs_range_dict = {}
        for cluster in range(n_clusters):
            qc_min, qc_max = (float(df_copy[df_copy["Cluster"] == cluster].QC.min()),
                              float(df_copy[df_copy["Cluster"] == cluster].QC.max()))
            qcs_range_dict[cluster] = qc_min, qc_max

        print(f"eps: {eps}, n_clusters: {n_clusters}")
        for cluster, qcs_range in sorted(qcs_range_dict.items(), key=lambda value: value[1]):
            items_in_cluster = len(df_copy[df_copy["Cluster"] == cluster])
            print(f"cluster: {cluster}, qcs_range: {qcs_range}, items_in_cluster: {items_in_cluster}")
        print()  # empty space for separator


def cluster_analysis(df: pd.DataFrame, eps: float) -> pd.DataFrame:
    """
    Performs cluster analysis (DBSCAN) based on defined eps

    :param df:
    :param eps:
    :return: table updated with clusters
    """
    qcs_values = df.QC.values.tolist()

    dbscan = DBSCAN(eps=eps, min_samples=1)
    clustering = dbscan.fit(np.asarray(qcs_values).reshape(-1, 1))
    labels = clustering.labels_

    df_copy = df.copy()
    df_copy["Cluster"] = labels

    return df_copy


def save_seq(save_dir: str, folder: str, subfolder: str, acc: str, strand: int, seq_start: int, seq_stop: int) -> None:
    """
    Saves sequence obtained with Entrez.efetch in the specified folder
    Additional function for save_seqs

    :param save_dir:
    :param folder:
    :param subfolder:
    :param acc:
    :param strand:
    :param seq_start:
    :param seq_stop:
    :return:
    """

    save_path = os.path.join(save_dir, folder, subfolder)

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(f"{save_path}/{acc}.fa", "w") as outfile:
        seq = Entrez.efetch(db="nucleotide", id=acc, strand=strand, seq_start=seq_start, seq_stop=seq_stop,
                            rettype="fasta").read()
        outfile.write(seq)


def save_single_seq(save_dir, filename: str,
                    acc: str, strand: int, seq_start: int, seq_stop: int) -> None:
    """
    Saves single sequence obtained with Entrez.efetch in the specified folder

    :param save_dir:
    :param filename:
    :param acc:
    :param strand:
    :param seq_start:
    :param seq_stop:
    :return:
    """

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    with open(f"{save_dir}/{filename}.fa", "w") as outfile:
        seq = Entrez.efetch(db="nucleotide", id=acc, strand=strand, seq_start=seq_start, seq_stop=seq_stop,
                            rettype="fasta").read()
        outfile.write(seq)


def find_codon(seq: str, which: str, frame_shift: int = 0,
               sep: str = "-", print_seq: bool = False) -> int:
    """
    Finds start or stop codon position of given sequence.
    Additional function for choose_best_frameshift

    :param seq:
    :param which:
    :param frame_shift:
    :param sep:
    :param print_seq:
    :return:
    """
    codons = {
        "start": ["ATG"],
        "stop": ["TAA", "TAG", "TGA"]
    }
    seq_list = []
    codon_pos = None

    for i_nt in range(frame_shift, len(seq), 3):
        codon = seq[i_nt:i_nt + 3]
        seq_list.append(codon)
        if codon in codons[which]:
            codon_pos = i_nt - frame_shift
            break

    if print_seq:
        print(f"{sep}".join(seq_list))

    return codon_pos


def choose_best_frameshift(seq: str) -> list:
    """
    Chooses best frameshifts for start and stop codons of given sequence based on total sequence length

    :param seq:
    :return:
    """
    start_codons = {}
    stop_codons = {}

    for frame_shift in range(0, 3):
        start_codon = find_codon(seq, which="start", frame_shift=frame_shift)
        start_codons[frame_shift] = start_codon + frame_shift  # !!! важный момент

    i = 0
    for start_codon in start_codons.values():
        for frame_shift in range(0, 3):
            i += 1
            stop_codon = find_codon(seq[start_codon:], which="stop", frame_shift=frame_shift)
            total_shift = start_codon + stop_codon+ frame_shift + 3
            seq_slice_len = len(seq[start_codon:total_shift])
            stop_codons[i] = [start_codon, total_shift, seq_slice_len]

    # print
    for key, value in stop_codons.items():
        print(f"Variant {key}: START: {value[0]}, STOP: {value[1]}, LEN: {value[2]}")

    max_length = 0
    best_variant = 1
    for variant, (start, stop, length) in stop_codons.items():
        if length > max_length:
            max_length = length
            best_variant = variant

    return stop_codons[best_variant]


def save_seqs(df: pd.DataFrame, folder: str,
              save_dir: str = "../Sequences") -> None:
    """
    Saves sequences from dataframe

    :param df:
    :param folder:
    :param save_dir:
    :return:
    """
    for cluster in range(len(df.Cluster.unique())):
        subset_df = df.loc[df["Cluster"] == cluster]
        min_value = str(df.loc[df["Cluster"] == cluster].min().QC.round(2))
        max_value = str(df.loc[df["Cluster"] == cluster].max().QC.round(2))
        subfolder = f"{min_value}-{max_value}"

        for hit in subset_df.index:
            line = df.loc[hit]
            acc, strand, seq_start, seq_stop = line.Acc, line.Strand, line.Start, line.Stop
            if seq_stop - seq_start >= 10_000:  # FIXME
                time.sleep(0.333333334)
                save_seq(save_dir, folder, subfolder, acc, strand, seq_start, seq_stop)


def extract_genome_coverage(acc: str) -> str:
    """
    Extracts genome coverage with Entrez.efetch using accession number from NCBI.
    Additional function for extract_genome_coverages

    :param acc:
    :return:
    """
    with Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text") as handle:
        genbank_data = handle.read()

    try:
        genome_coverage = [line for line in genbank_data.split("\n") if "Genome Coverage" in line][0]
    except IndexError:
        genome_coverage = "            Genome Coverage        :: 0x"

    return genome_coverage


def extract_genome_coverages(df: pd.DataFrame) -> list:
    """
    Extracts genome coverages for Acc column in the provided dataframe

    :param df:
    :return:
    """
    accessions = df.Acc.values.tolist()
    genome_coverages = []

    for acc in accessions:
        time.sleep(0.333333334)
        genome_coverage = extract_genome_coverage(acc)
        genome_coverages.append(genome_coverage)

    return genome_coverages


def add_genome_coverages(genome_coverages: list, df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds genome coverages to the provided dataframe parsing coverages to float format

    :param genome_coverages:
    :param df:
    :return:
    """
    for i in range(len(genome_coverages)):
        genome_coverages[i] = float(genome_coverages[i].split(":: ")[-1].split("x")[0])

    df["Genome_Coverage"] = genome_coverages

    return df


def select_max_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Updates dataframe with maximum genome coverage accessions - deletes duplicated species

    :param df:
    :return:
    """
    n_clusters = df.Cluster.nunique()
    max_ids = []

    for cluster in range(n_clusters):
        cluster_species = df[df["Cluster"] == cluster].Species_name.unique().tolist()
        for species in cluster_species:
            species_max_genome_coverage = df.query(
                "Cluster == @cluster & Species_name == @species").Genome_Coverage.idxmax()
            max_ids.append(species_max_genome_coverage)

    df = df.loc[max_ids]

    return df


def filter_genome_coverages(df: pd.DataFrame, genome_coverage_threshold: int = 50) -> pd.DataFrame:
    """
    Filters genome coverages based on genome coverage threshold

    :param df:
    :param genome_coverage_threshold:
    :return:
    """
    df = df[df["Genome_Coverage"] >= genome_coverage_threshold]

    return df
