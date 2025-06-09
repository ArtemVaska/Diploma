import os
import shutil

import pandas as pd


def extract_species_names(dir_name: str, folder_name: str) -> dict:
    """
    Creates dictionary subfolder: [species_name_acc] from specified folder in Sequences dir.
    Additional function for group_species function

    :param folder_name:
    :param dir_name:
    :return:
    """
    path = os.path.join(dir_name, folder_name)

    species_dict = {}
    for subfolder in os.listdir(path):
        subfolder_path = os.path.join(path, subfolder)
        species_list = []
        for file in os.listdir(subfolder_path):
            acc = file.replace(".fa", "")
            with open(os.path.join(subfolder_path, file), "r") as infile:
                for line in infile:
                    species_name = "_".join(line.split(",")[0].split()[1:])
                    species_name_acc = f"{species_name}_{acc}"
                    species_list.append(species_name_acc)
                    break
        species_dict[subfolder] = species_list

    return species_dict


def extract_unique_species(species_dict: dict, subfolder_name: str) -> list:
    """
    Extracts unique species from a dictionary subfolder.
    Additional function for group_species function

    :param species_dict:
    :param subfolder_name:
    :return:
    """
    unique_species = []
    for species_name_acc in species_dict[subfolder_name]:
        species_name = "_".join(species_name_acc.split("_")[:-1])
        if species_name not in unique_species:
            unique_species.append(species_name)

    return unique_species


def group_species(df: pd.DataFrame, folder_name: str, new_folder_name: str, dir_name: str = "../Sequences") -> None:
    """
    Groups species from specified folder by species names to the new specified folder.
    Also adds new column to df with species_name

    :param df:
    :param folder_name:
    :param new_folder_name:
    :param dir_name:
    :return:
    """
    species_dict = extract_species_names(dir_name, folder_name)
    df["Species_name"] = ""

    for subfolder in species_dict:
        unique_species = extract_unique_species(species_dict, subfolder)

        for species in unique_species:
            species_path = os.path.join(dir_name, new_folder_name, subfolder, species)
            if not os.path.exists(species_path):
                os.makedirs(species_path)

            for species_name_acc in species_dict[subfolder]:
                species_name = "_".join(species_name_acc.split("_")[:-1])
                if species_name == species:
                    acc = species_name_acc.split("_")[-1]
                    df.loc[df["Acc"] == acc, "Species_name"] = species_name

                    file_name = f"{acc}.fa"
                    file_to_copy = os.path.join(dir_name, folder_name, subfolder, file_name)
                    shutil.copy(file_to_copy, species_path)


def create_cluster_subfolder_name(df: pd.DataFrame, cluster: int) -> str:
    """
    Creates cluster subfolder name based on cluster number from provided dataframe
    Additional function for group_species_genome_coverage function

    :param df:
    :param cluster:
    :return:
    """
    min_value = str(df.loc[df["Cluster"] == cluster].min().QC.round(2))
    max_value = str(df.loc[df["Cluster"] == cluster].max().QC.round(2))
    subfolder_name = f"{min_value}-{max_value}"

    return subfolder_name


def group_species_genome_coverage(df: pd.DataFrame, folder_name: str, new_folder_name: str,
                                  dir_name: str = "../Sequences", ) -> None:
    """
    Groups species from folder to a new folder by genome coverage based on provided dataframe.
    Also renames filenames for convenience

    :param df:
    :param folder_name:
    :param new_folder_name:
    :param dir_name:
    :return:
    """
    folder_path = os.path.join(dir_name, folder_name)  # ex: ../Sequences/Drosophilidae
    subfolders = os.listdir(folder_path)

    n_clusters = df.Cluster.nunique()

    for cluster in range(n_clusters):
        new_subfolder_name = create_cluster_subfolder_name(df, cluster)  # name of cluster, ex: 0.16-0.39

        # path to filtered cluster, ex: ../Sequences/Drosophilidae_filtered/0.16-0.39
        new_folder_path = os.path.join(dir_name, new_folder_name, new_subfolder_name)
        os.makedirs(new_folder_path, exist_ok=True)

        accessions = df[df["Cluster"] == cluster].Acc.tolist()
        for accession in accessions:
            filename = f"{accession}.fa"

            for subfolder in subfolders:
                file_path = os.path.join(folder_path, subfolder, filename)
                if os.path.isfile(file_path):
                    shutil.copy(file_path, new_folder_path)

            species_name = df[df["Acc"] == accession].Species_name.iloc[0]
            new_filename = f"{species_name}_{accession}.fa"

            source = os.path.join(new_folder_path, filename)
            destination = os.path.join(new_folder_path, new_filename)
            os.rename(source, destination)
