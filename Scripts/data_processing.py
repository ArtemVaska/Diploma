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
        for cluster, qcs_range in sorted(qcs_range_dict.items(), key=lambda value:value[1]):
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


def save_seqs(folder: str, df: pd.DataFrame,
              save_dir: str = "../Sequences") -> None:
    """
    Saves sequences from dataframe

    :param folder:
    :param df:
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
            if seq_stop - seq_start >= 10_000:  #FIXME
                time.sleep(1)
                save_seq(save_dir, folder, subfolder, acc, strand, seq_start, seq_stop)
