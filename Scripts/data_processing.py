import os
import time
import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from Bio import Entrez

Entrez.email = "artem_vasilev_01@list.ru"


def cluster_analysis(qcs: dict, qc_threshold=0.1) -> pd.DataFrame:
    """
    Performs cluster analysis based on QC

    :param qcs:
    :param qc_threshold:
    :return: table with Acc., QC, Cluster
    """
    qcs_filtered = {acc: qc for acc, qc in qcs.items() if qc > qc_threshold}
    qcs_values = list(qcs_filtered.values())

    dbscan = DBSCAN(eps=0.04, min_samples=1)  # TODO eps automate calculate
    clustering = dbscan.fit(np.asarray(qcs_values).reshape(-1, 1))
    labels = clustering.labels_

    qcs_df = pd.DataFrame.from_dict(qcs_filtered, orient="index", columns=["QC"])
    qcs_df["Cluster"] = labels

    return qcs_df


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
              save_dir: str = "/home/artemvaska/Master_degree/Diploma/Sequences") -> None:
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
            time.sleep(1)
            save_seq(save_dir, folder, subfolder, acc, strand, seq_start, seq_stop)
            