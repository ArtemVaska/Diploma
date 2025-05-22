from typing import Any, List, Dict
from Bio.Blast import NCBIXML
import pandas as pd

def extract_accession(hit_id: str) -> str:
    return hit_id.rsplit("|", 2)[-2] if "|" in hit_id else hit_id

def merge_intervals(intervals: List[tuple[int, int]]) -> List[tuple[int, int]]:
    """Объединяет перекрывающиеся интервалы."""
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:  # перекрытие
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def parse_psiblast_xml(xml_path: str, max_len: int = 5000) -> pd.DataFrame:
    with open(xml_path) as f:
        try:
            record = next(NCBIXML.parse(f))  # Берём первую итерацию
        except StopIteration:
            return pd.DataFrame()  # Если файл пуст

    query_len = record.query_length
    hits = []  # type: List[Dict[str, Any]]

    for alignment in record.alignments:
        hsps = [hsp for hsp in alignment.hsps if hsp.align_length <= max_len]
        if not hsps:
            continue

        # Сбор покрытых интервалов на query
        query_intervals = [(min(hsp.query_start, hsp.query_end),
                            max(hsp.query_start, hsp.query_end)) for hsp in hsps]
        merged_query_intervals = merge_intervals(query_intervals)
        covered_len = sum(end - start + 1 for start, end in merged_query_intervals)
        qc = round(covered_len / query_len, 2)
        qc = min(qc, 1.0)  # Защита от >1

        # Координаты на subject
        sbjct_starts = [min(hsp.sbjct_start, hsp.sbjct_end) for hsp in hsps]
        sbjct_ends = [max(hsp.sbjct_start, hsp.sbjct_end) for hsp in hsps]
        sbjct_start = min(sbjct_starts)
        sbjct_end = max(sbjct_ends)

        # Средний процент идентичности
        identities = sum(hsp.identities for hsp in hsps)
        alignment_len = sum(hsp.align_length for hsp in hsps)
        per_ident = round(identities / alignment_len, 2) if alignment_len else 0.0

        hits.append({
            "Hit_id": alignment.hit_id,
            "Acc": extract_accession(alignment.hit_id),
            "QC": qc,
            "Per_Ident": per_ident,
            "Start": sbjct_start,
            "End": sbjct_end,
            "Sbjct_Len": sbjct_end - sbjct_start + 1,
        })

    df = pd.DataFrame(hits)
    if not df.empty:
        df.set_index("Hit_id", inplace=True)
    return df


def filter_psiblast_hits(
    df: pd.DataFrame,
    min_qc: float = 0.8,
    min_ident: float = 0.8,
    min_sbjct_len: int = 500
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()

    filtered_df = df[
        (df["QC"] >= min_qc) &
        (df["Per_Ident"] >= min_ident) &
        (df["Sbjct_Len"] >= min_sbjct_len)
    ]

    return filtered_df
