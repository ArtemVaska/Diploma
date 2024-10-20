import os

from Bio import Entrez

Entrez.email = "artem_vasilev_01@list.ru"


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
