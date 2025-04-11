import os
import time
import urllib.error

import pandas as pd

from Bio import Entrez
from tg_logger import telegram_logger


Entrez.email = "artemvaskaa@gmail.com"


def nucl_search(query: str) -> list:
    """
    Searches smth via Entrez.esearch in the nucleotide database

    query example: (Drosophilidae[ORGN] NOT Drosophila melanogaster[ORGN]) AND \
    (chromosome X[WORD] NOT PREDICTED[WORD] NOT gene[WORD]) AND 15000000:75000000[SLEN]
    """
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=100)
    record = Entrez.read(handle)
    return record["IdList"]


def save_esearch_results(id_list: list, folder: str) -> None:
    """
    Saves .fa files in specified folder from esearch results
    """
    SAVE_DIR = "/home/artemvaska/Master_degree/Diploma/References"
    save_path = os.path.join(SAVE_DIR, folder)

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    for seq_id in id_list:
        with open(f"{save_path}/{seq_id}.fa", "w") as outfile:
            seq = Entrez.efetch(db="nucleotide", id=seq_id, retmode="text", rettype="fasta").read()
            outfile.write(seq)


@telegram_logger(chat_id=611478740)
def select_phyla(df: pd.DataFrame, phyla: str) -> dict:
    """
    Selects all taxon_id from the table that match the given phylum.

    :param df:
    :param phyla:
    :return:
    """
    phyla_taxid = {phyla: []}

    for index, row in df.iterrows():
        time.sleep(0.333333334)
        try:
            stream = Entrez.efetch(db="Taxonomy", id=str(index), retmode="xml")
        except urllib.error.HTTPError:
            print(f"ERROR taxid: {index}, {row['org_name']}")
            continue
        records = Entrez.read(stream)
        if phyla in records[0]["Lineage"]:
            phyla_taxid[phyla].append(index)
            print(index, row["org_name"])

    return phyla_taxid
