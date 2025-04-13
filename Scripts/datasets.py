import time

from Bio import Entrez, SeqIO
from fasta_processing import read_single_fasta


def download_gene_gb(org_names: list) -> None:
    # Obtain mRNA Accession and ranges of the gene
    for org_name in org_names:
        with open(f"../Datasets/{org_name}/ncbi_dataset/data/gene.fna") as infile:
            line = infile.readline().rstrip()
            gene_acc = line.split(":")[0][1:]
            gene_ranges = line.split()[0].split(":")[1]
            gene_range_1 = gene_ranges.split("-")[0]
            gene_range_2 = gene_ranges.split("-")[1]

        if "c" in gene_range_1:
            gene_range_1 = gene_range_1.replace("c", "")
            strand = 2
        else:
            strand = 1
        gene_range_1, gene_range_2 = int(gene_range_1), int(gene_range_2)

        # Obtain and download mRNA in GenBank format
        stream = Entrez.efetch(db="nucleotide", id=gene_acc, idtype="acc",
                               seq_start=gene_range_1, seq_stop=gene_range_2, strand=strand,
                               rettype="gb", retmode="text")
        with open(f"../Datasets/{org_name}/ncbi_dataset/data/gene.gb", "w") as outfile:
            outfile.write(stream.read())

        print(f"{org_name}.gb has been downloaded")

        time.sleep(0.333333334)


def parse_exon_ranges(org_names: list, feature_type: str = "mRNA") -> dict:
    exons_dict = {}

    feature_types = ["mRNA", "CDS"]
    if feature_type not in feature_types:
        raise Exception(f"Feature type {feature_type} is not supported")

    for org_name in org_names:
        exons_range = {}
        records = SeqIO.parse(f"../Datasets/{org_name}/ncbi_dataset/data/gene.gb", "genbank")
        for record in records:
            for feature in record.features:
                if feature.type == feature_type:
                    for part_i, part in enumerate(feature.location.parts):
                        start, end = int(part.start), int(part.end)
                        exons_range[part_i] = [start, end]
        exons_dict[org_name] = exons_range

    return exons_dict


def create_exons(exon_ranges: dict) -> None:
    for org_name in exon_ranges.keys():
        gene_seq = read_single_fasta(f"../Datasets/{org_name}/ncbi_dataset/data/gene.fna")
        with open(f"../Datasets/{org_name}/ncbi_dataset/data/exons.fa", "w") as outfile:
            for exon_i, exon_range in exon_ranges[org_name].items():
                outfile.write(f">{exon_i}\n{gene_seq[exon_range[0]:exon_range[1]]}\n")
        print(f"{org_name} exons has been created")
