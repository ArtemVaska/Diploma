import time
import Bio.Blast
import pandas as pd


def filter_df(df, qc_threshold=0.1) -> pd.DataFrame:
    """
    Filters dataframe based on QS threshold

    :param df:
    :param qc_threshold:
    :return: filtered dataframe based on qc_threshold
    """
    df_filtered = df.query("QC >= @qc_threshold")

    return df_filtered


def extract_target_range(hit: Bio.Blast.Hit) -> tuple:
    """
    Extracts hit range in target

    :param hit:
    :return: list with coords [start, end]
    """
    target_range_list = []

    for hsp in hit:
        target_range_list.append(int(hsp.coordinates[0][0]))
        target_range_list.append(int(hsp.coordinates[0][-1]))

    start, stop = min(target_range_list), max(target_range_list)
    if (stop - start) > 40000:
        print("Target range cannot be calculated automatically. Please enter coordinates manually from the list below:")
        time.sleep(0.5)
        print(sorted(target_range_list))
        time.sleep(0.5)
        start, stop = [int(coord) for coord in input("Enter start and stop coordinates' indices separated by space: ").split()]

    return start, stop


def calculate_qc(blast_record: Bio.Blast.Record) -> pd.DataFrame:
    """
    Calculates query coverage for every hit in blast_record and creates table

    :param blast_record:
    :return: table target.id, QC
    """
    query_length = len(blast_record.query)
    qcs = {}

    for hit in blast_record:
        target_id = hit.target.id
        query_total_range = 0
        for hsp in hit:
            query_total_range += int(hsp.coordinates[1][-1]) - int(hsp.coordinates[1][0])
        qc = query_total_range / query_length
        qcs[target_id] = round(qc, 4)

    qcs_df = pd.DataFrame.from_dict(qcs, orient="index", columns=["QC"])

    return qcs_df


def check_strands(hit: Bio.Blast.Hit) -> (int, int):
    """
    Checks the strand of the target sequence
    Additional function for update_df

    :param hit:
    :return:
    """
    hsp = hit[0]

    query_strand = (
        1 if hsp.coordinates[1][0] <= hsp.coordinates[1][-1] else 2
    )
    target_strand = (
        1 if hsp.coordinates[0][0] <= hsp.coordinates[0][-1] else 2
    )
    return query_strand, target_strand


def update_df(df: pd.DataFrame, blast_record: Bio.Blast.Record) -> pd.DataFrame:
    """
    Adds additional columns to the dataframe (Acc, Strand, Start, Stop) based on blast_record

    :param df:
    :param blast_record:
    :return: updated dataframe
    """
    df["Acc"] = [acc.split("|")[-2] for acc in list(df.index)]

    strands = []
    starts = []
    stops = []

    for hit_id in df.index:
        hit = blast_record[hit_id]
        query_strand, target_strand = check_strands(hit)
        strands.append(target_strand)

        start, stop = extract_target_range(hit)
        starts.append(start)
        stops.append(stop)

    df["Strand"] = strands
    df["Start"] = starts
    df["Stop"] = stops

    return df
