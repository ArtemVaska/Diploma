import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN


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
