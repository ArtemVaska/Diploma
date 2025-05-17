import os
import time

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from sklearn.cluster import DBSCAN

from Scripts.fasta_processing import plain_to_fasta
from fasta_processing import read_fasta, read_single_fasta

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
        print(f"Frameshift {key}: Start: {value[0]}, Stop: {value[1]}, Length: {value[1] - value[0]}")

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
    handle = Entrez.esearch(db="nucleotide", term=search_term, idtype="acc", retmax=max_results)
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


def dict_align_create(phyla: str, org_names: list, align_type: str) -> dict:
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
    for org_name in org_names:
        dict_align[f"{org_name}"] = read_single_fasta(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/{filename}")

    return dict_align


def analyze_exons(path: str) -> pd.DataFrame:
    exons = read_fasta(path)

    # >0:0-183
    ranges = []
    for header in exons.keys():
        coords = header.split(":")[1]
        ranges.append(coords)

    lengths = [len(seq) for seq in exons.values()]
    seqs = exons.values()

    df = pd.DataFrame(
        {
            "length": lengths,
            "coords": ranges,
            "sequence": seqs,
        }
    )
    return df


def concat_2_exons(df: pd.DataFrame, org_name: str, indices: list):
    path = f"../Datasets/{org_name}/ncbi_dataset/data/"

    seq_0 = df.iloc[indices[0]].sequence
    seq_1 = df.iloc[indices[1]].sequence

    filename = f"{org_name}_{indices[0]}-{len(seq_0)}-{indices[1]}_{len(seq_1)}.fa"

    with open(f"{path}/{filename}", "w") as handle:
        handle.write(f">{org_name}\n")
        handle.write(f"{''.join([seq_0, seq_1])}\n")


def concat_cassette(cassette_dict: dict, concat_type: str) -> str | None:
    """
    Concatenates 2 exons from conservative cassette or a conservative cassette from dict to a single string.
    :param cassette_dict:
    :param concat_type:
    :return:
    """
    match concat_type:
        case "i":
            return "".join([seq for header, seq in cassette_dict.items() if header == 'cassette'])
        case "ee":
            return "".join([seq for header, seq in cassette_dict.items() if header != 'cassette'])
        case "eie":
            return "".join([seq for header, seq in cassette_dict.items()])
    return None


def create_cassette(phyla: str, org_name: str, df_exons: pd.DataFrame, exons_i: list) -> dict:
    path = f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data"
    gene_fna = read_single_fasta(f"{path}/gene.fna")

    exon_0 = df_exons.loc[exons_i[0]]
    exon_1 = df_exons.loc[exons_i[1]]

    cassette_start = int(exon_0.coords.split("-")[1]) + 1
    cassette_end = int(exon_1.coords.split("-")[0])
    cassette = gene_fna[cassette_start:cassette_end]

    cassette_dict = {
        f"{exons_i[0]}:{exon_0.coords}": exon_0.sequence,
        "cassette": cassette,
        f"{exons_i[1]}:{exon_1.coords}": exon_1.sequence,
    }
    with open(f"{path}/cassette.fa", "w") as outfile:
        for header, seq in cassette_dict.items():
            outfile.write(f">{header}\n"
                          f"{seq}\n")

    with open(f"{path}/cds.fna") as infile:
        line = infile.readline()
        cds_header = line.rstrip()
    cds = read_single_fasta(f"{path}/cds.fna")
    cds_cassette = cds.replace(f"{exon_0.sequence}{exon_1.sequence}",
                               f"{exon_0.sequence}{cassette}{exon_1.sequence}")
    with open(f"{path}/cds_cassette_plain.fa", "w") as outfile:
        outfile.write(f"{cds_header} [+cassette_intron]\n")
        outfile.write(f"{cds_cassette}\n")
    plain_to_fasta(f"{path}/cds_cassette_plain.fa", 70)

    return cassette_dict


def dict_align_info_analyze(dict_align_info: dict, feature: str) -> (pd.DataFrame, dict):
    dict_align = {}
    rows = []

    for phyla, org_names in dict_align_info.items():
        dict_align[phyla] = dict_align_create(phyla, org_names, feature)

    for phyla, org_name_seq in dict_align.items():
        for org_name, seq in org_name_seq.items():
            stop_codon_pos = find_codon(seq, which="stop")
            cassette_intron = read_fasta(f"../Datasets/{phyla}/{org_name}/ncbi_dataset/data/cassette.fa")["cassette"]
            cassette_intron_start = seq.find(cassette_intron)
            rows.append(
                {
                    "phylum": phyla,
                    "org_name": org_name,
                    "stop_codon_pos": stop_codon_pos,
                    "equal_to_cds": stop_codon_pos + 3 == len(seq),
                    "cassette_start": cassette_intron_start,
                    "intron_len_to_stop_codon": stop_codon_pos - cassette_intron_start,
                }
            )
    df = pd.DataFrame(rows)

    # formatted org_name for dict_align
    new_dict = {}
    for phyla, org_name_seq in dict_align.items():
        new_subdict = {}
        for org_name, seq in org_name_seq.items():
            formatted_name = "_".join(org_name.split("_")[:-1]).capitalize()
            new_subdict[formatted_name] = seq
        new_dict[phyla] = new_subdict

    return df, new_dict


def dict_align_update_keys(dict_align: dict):
    new_org_names = []
    for org_name, seq in dict_align.items():
        new_org_names.append("_".join(org_name.split("_")[:-1]).capitalize())
    new_dict_align = dict(zip(new_org_names, dict_align.values()))
    return new_dict_align
