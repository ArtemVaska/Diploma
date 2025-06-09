import os
import subprocess

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from fasta_processing import read_single_fasta, read_fasta


def run_maxentscan(phylum: str, org_names: list) -> pd.DataFrame:
    prefix = "../datasets"
    postfix = "ncbi_dataset/data"
    site_path = "../maxentscan_output"

    rows = []
    for org_name in org_names:
        path_file = f"{prefix}/{phylum}/{org_name}/{postfix}"
        sequence = read_single_fasta(f"{path_file}/cds_cassette.fa")
        intron = read_fasta(f"{path_file}/cassette.fa")["cassette"]
        print(f">{org_name.rsplit("_", 1)[0].capitalize()}\n{intron}")

        # exon 3 | 6 intron ... intron 20 | 3 exon
        donor_site = sequence[sequence.find(intron) - 3:sequence.find(intron) + 6]  # 9, 5'
        acceptor_site = sequence[
                        sequence.find(intron) + len(intron) - 20:sequence.find(intron) + len(intron) + 3]  # 23, 3'

        for site_type in ["donor", "acceptor"]:
            match site_type:
                case "donor":
                    site = donor_site
                    cmd = "maxentscan_score5.pl"
                case "acceptor":
                    site = acceptor_site
                    cmd = "maxentscan_score3.pl"

            site_dir = f"{site_path}/{phylum}/{org_name}"
            if not os.path.exists(site_dir):
                os.makedirs(site_dir, exist_ok=True)

            with open(f"{site_dir}/{site_type}.fa", "w") as outfile:
                outfile.write(f"{site}\n")

            result = subprocess.run(
                [cmd] + [f"{site_dir}/{site_type}.fa"],
                capture_output=True,
            )
            lines = result.stdout.decode("utf-8").strip().split("\t")
            rows.append(
                {
                    "org_name": org_name,
                    "site_type": site_type,
                    "site_seq": lines[0],
                    "maxentscore": lines[1],
                }
            )
            df = pd.DataFrame(rows)
            df["org_name"] = df["org_name"].apply(
                lambda name: name.rsplit("_", 1)[0].capitalize()
            )
            df.to_csv(f"{site_path}/{phylum}/result.tsv", sep="\t", index=False)
    return df


def run_maxentscan_all_introns(phylum: str, org_name: str) -> pd.DataFrame:
    full_seq = read_single_fasta(f"../datasets/{phylum}/{org_name}/ncbi_dataset/data/gene.fna")
    exons_dict = read_fasta(f"../datasets/{phylum}/{org_name}/ncbi_dataset/data/exons.fa")

    # Creating list of exons and exons_coords for introns parsing
    exons = []
    exons_coords = []
    for key, exon_seq in exons_dict.items():
        _, exon_coords = key.split(":")
        exons.append(exon_seq)
        exons_coords.append([int(coord) for coord in exon_coords.split("-")])

    # Creates list of introns from parsed coords
    introns = []
    for intron_i in range(len(exons_coords)-1):
        intron_seq = full_seq[exons_coords[intron_i][1] : exons_coords[intron_i+1][0]]
        introns.append(intron_seq)

    print(f"{org_name}: {[len(intron) for intron in introns]}")

    # Creates lists of donors and acceptors for maxentscan processing
    donors = []
    acceptors = []

    for i in range(len(exons)-1):
        donor = exons[i][-3:] + introns[i][:6]
        acceptor = introns[i][-20:] + exons[i+1][:3]
        donors.append(donor)
        acceptors.append(acceptor)

    # Saves donors and acceptors to files
    dir_path = f"../maxentscan_output/{phylum}_full_gene/{org_name}"
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)

    with open(f"{dir_path}/donors.fa", "w") as outfile:
        for donor in donors:
            outfile.write(f"{donor}\n")

    with open(f"{dir_path}/acceptors.fa", "w") as outfile:
        for acceptor in acceptors:
            outfile.write(f"{acceptor}\n")

    # Runs maxentscan
    result = subprocess.run(
        ["maxentscan_score5.pl"] + [f"{dir_path}/donors.fa"],
        capture_output=True,
    )
    lines_donors = result.stdout.decode("utf-8").strip().split("\n")
    lines_donors_parsed = [score.split("\t") for score in lines_donors]

    result = subprocess.run(
        ["maxentscan_score3.pl"] + [f"{dir_path}/acceptors.fa"],
        capture_output=True,
    )
    lines_acceptors = result.stdout.decode("utf-8").strip().split("\n")
    lines_acceptors_parsed = [score.split("\t") for score in lines_acceptors]

    # Creates pd.DataFrame with results
    rows = []
    for (donor_seq, donor_score), (acceptor_seq, acceptor_score) in zip(lines_donors_parsed, lines_acceptors_parsed):
        rows.append(
            {
                "donor_seq": donor_seq,
                "donor_score": donor_score,
                "acceptor_score": acceptor_score,
                "acceptor_seq": acceptor_seq,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(f"{dir_path}/result.tsv", sep="\t", index=False)

    return df


def maxentscan_boxplot(phylum: str, df: pd.DataFrame):
    plt.figure(figsize=(8, 6))

    # Boxplot без hue, чтобы всё было по центру
    sns.boxplot(data=df, x="site_type", y="maxentscore",
                hue="site_type", palette="Set2", legend=False)

    # Добавляем точки поверх боксплота, строго по центру каждой категории
    sns.stripplot(data=df, x="site_type", y="maxentscore",
                  color="black", size=6, jitter=False, marker="o", edgecolor="white", linewidth=0.5)

    plt.title(f"MaxEntScore по сайтам сплайсинга у {phylum}", fontsize=14)
    plt.ylabel("MaxEntScore")
    plt.xlabel("Тип сайта")
    plt.tight_layout()
    plt.show()
