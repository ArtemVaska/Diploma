import shutil
import os


def edit_names_for_alignment(folder_name: str, dir_name: str = "../Sequences") -> None:
    """
    Renames the 1st line of genome.fa in specified folder to paste in multiple alignment tools faster

    :param folder_name:
    :param dir_name:
    :return:
    """
    subfolders = os.listdir(os.path.join(dir_name, folder_name))

    for subfolder in subfolders:
        source = os.path.join(dir_name, folder_name, subfolder)
        new_subfolder = f"{subfolder}_for_alignment"
        destination = os.path.join(dir_name, folder_name, new_subfolder)
        shutil.copytree(source, destination)

        files = os.listdir(destination)
        for file in files:
            file_path = os.path.join(dir_name, folder_name, new_subfolder, file)
            with open(file_path, "r") as infile:  # FIXME errors
                lines = infile.readlines()

            line_to_edit = lines[0].strip()
            species_name = "_".join(line_to_edit.split()[1:3])
            line_edit = f">{species_name}\n"
            lines[0] = line_edit

            with open(file_path, "w") as outfile:
                outfile.writelines(lines)
