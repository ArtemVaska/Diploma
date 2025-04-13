import os
import time
import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from subprocess import DEVNULL, STDOUT, check_call

from fasta_processing import read_single_fasta

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


def choose_best_frameshift(seq: str, translate: bool = False) -> str | list:
    """
    Chooses best frameshifts for start and stop codons of given sequence based on total sequence length.
    Also slices mRNA by found codon positions and translates to protein sequence if translate flag is True.

    :param seq:
    :param translate:
    :return:
    """
    codons = {}

    for frame_shift in range(0, 3):
        start_codon = find_codon(seq, which="start", frame_shift=frame_shift)
        codons[frame_shift] = [start_codon + frame_shift]  # !!! важный момент

    for frame_shift in range(0, 3):
        start_codon = codons[frame_shift][0]
        stop_codon = find_codon(seq[start_codon:], which="stop")
        codons[frame_shift].append(start_codon + stop_codon + 3)

    for key, value in codons.items():
        print(f"Frameshift {key}: Start: {value[0]}, Stop: {value[1]}, Length: {value[1]-value[0]}")

    max_length = 0
    best_variant = 1
    for variant, (start, stop) in codons.items():
        length = stop - start + 1
        if length > max_length:
            max_length = length
            best_variant = variant

    if translate:
        start_pos, stop_pos = codons[best_variant][0], codons[best_variant][1]
        sliced_seq = seq[start_pos:stop_pos]
        protein = str(Seq(sliced_seq).translate())
        return protein

    return codons[best_variant]


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


def download_transcripts_by_geneid(gene_id: str, org_name: str, max_results: int = 5) -> str:
    """
    Downloads XM_... mRNA transcripts from Entrez via GeneID to the specified folder.
    Additional function for save_subset_df_transcripts.

    :param gene_id:
    :param org_name:
    :param max_results:
    :return:
    """
    # Step 1: Search for the GeneID and get linked mRNA nucleotide IDs
    # search_term = f"{gene_id}[GeneID]"  # Search query using GeneID
    search_term = f"{gene_id}[GeneID]"
    handle = Entrez.esearch(db="nucleotide", term=search_term,  idtype="acc", retmax=max_results)
    search_results = Entrez.read(handle)
    handle.close()

    # Get list of Nucleotide IDs (mRNA sequences)
    ids = []
    for acc in search_results["IdList"]:
        if "XM_" in acc:
            ids.append(acc)
    print(f"Found {len(ids)} XM_mRNA sequences for GeneID {gene_id}.")

    save_path = f"../References/{org_name}"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Step 2: Fetch the nucleotide sequences (mRNA) from NCBI
    seq_records = []
    for seq_id in ids:
        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(fetch_handle, "fasta")
        seq_records.append(seq_record.id)
        fetch_handle.close()

        filename = f"{seq_record.id}.fasta"
        with open(f"{save_path}/{filename}", "w") as out_file:
            SeqIO.write(seq_record, out_file, "fasta")
        print(f"Saved transcript to {save_path}/{filename}")

    return seq_records[0]


def save_subset_df_transcripts(df: pd.DataFrame) -> dict:
    """
    Saves mRNA transripts by given GeneID.

    :param df: pd.DataFrame
    :return: dict org_name:seq
    """
    seq_dict = {}

    for index, row in df.iterrows():
        folder_name = row["org_name"].lower().replace(" ", "_")
        filename = download_transcripts_by_geneid(str(row["gene_id"]), folder_name, max_results=3)
        print()
        seq = read_single_fasta(f"../References/{folder_name}/{filename}.fasta")
        seq_dict[folder_name] = seq

    return seq_dict


def download_subset_df_datasets(df: pd.DataFrame) -> None:
    """
    Downloads gene, rna and protein for every GeneID in the provided dataframe
    via ncbi-datasets-cli.

    :param df: Pandas DataFrame
    :return: None
    """

    for index, row in df.iterrows():
        gene_id = str(row["gene_id"])
        org_name = row["org_name"].lower().replace(" ", "_")
        shell_commands = [
            ["datasets", "download", "gene", "gene-id", gene_id,
             "--include", "gene,cds,rna,protein",
             "--filename", f"../Datasets/{org_name}.zip"],
            ["unzip", f"../Datasets/{org_name}.zip",
             "-d", f"../Datasets/{org_name}"],
            ["rm", "-r", f"../Datasets/{org_name}.zip"]
        ]
        for command in shell_commands:
            check_call(command, stdout=DEVNULL, stderr=STDOUT)
        print(f"Files for {org_name} downloaded successfully")
