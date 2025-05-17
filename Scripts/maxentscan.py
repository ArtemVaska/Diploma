import os
import subprocess

import pandas as pd

from fasta_processing import read_single_fasta, read_fasta


def run_maxentscan(phylum: str, org_names: list) -> pd.DataFrame:
    prefix = "../Datasets"
    postfix = "ncbi_dataset/data"
    site_path = "../maxentscan"

    rows = []
    for org_name in org_names:
        path_file = f"{prefix}/{phylum}/{org_name}/{postfix}"
        sequence = read_single_fasta(f"{path_file}/cds_cassette.fa")
        intron = read_fasta(f"{path_file}/cassette.fa")["cassette"]

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
    return df
