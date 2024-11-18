import Bio.Blast


def extract_target_range(hit: Bio.Blast.Hit) -> list:
    """
    Extracts hit range in target

    :param hit:
    :return: list with coords [start, end]
    """
    target_range_list = []

    for hsp in hit:
        target_range_list.append(int(hsp.coordinates[0][0]))
        target_range_list.append(int(hsp.coordinates[0][-1]))

    target_range = [min(target_range_list), max(target_range_list)]

    return target_range


def calculate_qc(blast_record: Bio.Blast.Record) -> dict:
    """
    Calculates query coverage for every hit in blast_record

    :param blast_record:
    :return: dict target.id:QC
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

    return qcs


def check_strands(hsp: Bio.Blast.HSP) -> (int, int):
    """
    Checks the strand of the target sequence
    Additional function for update_df

    :param hsp:
    :return:
    """
    query_strand = (
        1 if hsp.coordinates[1][0] <= hsp.coordinates[1][-1] else 2
    )
    target_strand = (
        1 if hsp.coordinates[0][0] <= hsp.coordinates[0][-1] else 2
    )
    return query_strand, target_strand
