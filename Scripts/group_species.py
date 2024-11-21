import os
import shutil


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
                    species_name = "_".join(line.split()[1:3])
                    species_name_acc = "_".join([species_name, acc])
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
        species_name = "_".join(species_name_acc.split("_")[:2])
        if species_name not in unique_species:
            unique_species.append(species_name)

    return unique_species


def group_species(folder_name: str, new_folder_name: str, dir_name: str = "../Sequences") -> None:
    """
    Groups species from specified folder by species names to the new specified folder

    :param folder_name:
    :param new_folder_name:
    :param dir_name:
    :return:
    """
    species_dict = extract_species_names(dir_name, folder_name)

    for subfolder in species_dict:
        unique_species = extract_unique_species(species_dict, subfolder)

        for species in unique_species:
            species_path = os.path.join(dir_name, new_folder_name, subfolder, species)
            if not os.path.exists(species_path):
                os.makedirs(species_path)

            for species_name_acc in species_dict[subfolder]:
                species_name = "_".join(species_name_acc.split("_")[:2])
                if species_name == species:
                    file_name = f"{species_name_acc.split("_")[-1]}.fa"
                    file_to_copy = os.path.join(dir_name, folder_name, subfolder, file_name)
                    shutil.copy(file_to_copy, species_path)
