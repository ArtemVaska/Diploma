import os
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from fasta_processing import read_single_fasta, read_fasta


def run_maxentscan(df: pd.DataFrame, sub_phylum: str) -> pd.DataFrame:

    # путь для всего для maxentscan'а
    site_path = "../maxentscan"

    # строки для заполнения результатов
    rows_maxentscan = []

    # сразу создаем папку, куда будем сохранять доноров и акцепторов
    os.makedirs(f"{site_path}/{sub_phylum}", exist_ok=True)

    # создаем доноров и акцепторов
    for row in df.iterrows():

        # определяем пути для получения данных из файлов в зависимости от источника
        if row[1].source == "datasets":
            prefix = "../Datasets"
            postfix = "ncbi_dataset/data/"  # здесь обязательно еще один слеш прописываем
            intron_name = "cassette"  # FIXME
            org_name = row[1].org_name_protein_id
            new_org_name = org_name.rsplit("_", 1)[0].capitalize()
        elif row[1].source == "psi_blast":
            prefix = "../Sequences_protein_id"
            postfix = ""  # а тут нет, потому что да
            intron_name = "cassette_intron"  # FIXME
            org_name = row[1].org_name_protein_id.split("__")[1]
            new_org_name = row[1].org_name_protein_id.split("__")[0]
        else:
            raise ValueError("Unknown source")


        input_file = f"{prefix}/{sub_phylum}/{org_name}/{postfix}cds_cassette.fa"  # тут не ставим доп. слеш
        cds_cassette_seq = read_single_fasta(input_file)
        intron = read_fasta(f"{prefix}/{sub_phylum}/{org_name}/{postfix}cassette.fa")[intron_name]  # и тут

        # exon 3 | 6 intron ... intron 20 | 3 exon
        donor_site = cds_cassette_seq[cds_cassette_seq.find(intron) - 3:cds_cassette_seq.find(intron) + 6]  # 9, 5'
        acceptor_site = cds_cassette_seq[
                        cds_cassette_seq.find(intron) + len(intron) - 20:cds_cassette_seq.find(intron) + len(intron) + 3]  # 23, 3'

        # определяем команды для запуска maxentscan
        for site_type in ["donor", "acceptor"]:
            match site_type:
                case "donor":
                    site = donor_site
                    cmd = "maxentscan_score5.pl"
                case "acceptor":
                    site = acceptor_site
                    cmd = "maxentscan_score3.pl"

            # записываем донора и акцептора в fasta-файл, т.к. maxentscan не читает строки в качестве инпута
            with open(f"{site_path}/{sub_phylum}/{new_org_name}_{site_type}.fa", "w") as outfile:
                outfile.write(f"{site}\n")

            result = subprocess.run(
                [cmd] + [f"{site_path}/{sub_phylum}/{new_org_name}_{site_type}.fa"],
                capture_output=True,
            )
            lines = result.stdout.decode("utf-8").strip().split("\t")
            rows_maxentscan.append(
                {
                    "sub_phylum": sub_phylum,
                    "org_name": new_org_name,
                    "site_type": site_type,
                    "site_seq": lines[0],
                    "maxentscore": lines[1],
                }
            )
            df = pd.DataFrame(rows_maxentscan)
            # сохраняем результат в папку с сабфилой, потом все прочитаем и будем анализировать
            df.to_csv(f"{site_path}/{sub_phylum}/{sub_phylum}_maxentscan.tsv", sep="\t", index=False)

    # возвращаем таблицу с результатом, если захочется посмотреть 1 филу
    return df


def maxentscan_boxplot_summary(df: pd.DataFrame, figsize: tuple = (14, 7)):
    # группировка по sub_phylum и подсчёт количества уникальных видов
    group_counts = df.groupby("sub_phylum")["org_name"].nunique()
    # переименовываем sub_phylum в формате "Имя (n = число)"
    df = df.copy()  # чтобы не менять исходный DataFrame
    df["sub_phylum_labeled"] = df["sub_phylum"].apply(
        lambda sp: f"{sp} (n = {group_counts.get(sp, 0)})"
    )

    # автоматическое определение пределов и тиков
    y_min = np.floor(df["maxentscore"].min() / 3) * 3
    y_max = np.ceil(df["maxentscore"].max() / 3) * 3
    yticks = np.arange(y_min, y_max + 1, 3)

    plt.figure(figsize=figsize)

    sns.boxplot(
        data=df,
        x="sub_phylum_labeled",
        y="maxentscore",
        hue="site_type"
    )

    plt.xlabel("Таксономическая группа", fontsize=12)
    plt.ylabel("MaxEntScan Score", fontsize=12)
    plt.title("Распределение MaxEntScan score по типам сайтов и таксономическим группам", fontsize=14)
    plt.xticks(rotation=45)
    plt.legend(title="Тип сайта")

    # применяем шаг 3 по оси Y
    plt.yticks(yticks)

    plt.tight_layout()
    plt.show()


def maxentscan_boxplot_single_group(
        df: pd.DataFrame,
        sub_phylum: str,
        figsize: tuple = (8, 6)
):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # фильтруем по выбранной группе
    df_group = df[df["sub_phylum"] == sub_phylum].copy()

    if df_group.empty:
        print(f"Нет данных для подтипа '{sub_phylum}'")
        return

    # добавим метку с числом уникальных видов
    n_species = df_group["org_name"].nunique()
    label = f"{sub_phylum} (n = {n_species})"
    df_group["sub_phylum_labeled"] = label

    # границы и тики по Y
    y_min = np.floor(df_group["maxentscore"].min() / 3) * 3
    y_max = np.ceil(df_group["maxentscore"].max() / 3) * 3
    yticks = np.arange(y_min, y_max + 1, 3)

    # рисуем график
    plt.figure(figsize=figsize)

    sns.boxplot(
        data=df_group,
        x="sub_phylum_labeled",
        y="maxentscore",
        hue="site_type"
    )

    plt.xlabel("Таксономическая группа", fontsize=12)
    plt.ylabel("MaxEntScan Score", fontsize=12)
    plt.title(f"Распределение MaxEntScan score в группе {sub_phylum}", fontsize=14)
    plt.yticks(yticks)
    plt.legend(title="Тип сайта")
    plt.tight_layout()
    plt.show()
