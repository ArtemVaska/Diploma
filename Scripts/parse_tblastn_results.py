from Bio.Blast import NCBIXML
import pandas as pd
from collections import Counter
from typing import List


def calculate_query_coverage(hsps: list) -> int:
    covered_ranges = []
    for hsp in hsps:
        covered_ranges.append((hsp.query_start, hsp.query_end))

    merged = []
    for start, end in sorted(covered_ranges):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    total_covered = sum(end - start + 1 for start, end in merged)
    return total_covered


def extract_accession(index_str: str) -> str:
    return index_str.rsplit("|", 2)[-2]


def parse_tblastn_xml(xml_path: str) -> pd.DataFrame:
    records: List[dict] = []
    indices: List[str] = []

    with open(xml_path) as handle:
        blast_records = NCBIXML.parse(handle)
        for record in blast_records:
            query_len = record.query_length

            for alignment in record.alignments:
                hsps = alignment.hsps
                if not hsps:
                    continue

                qc = calculate_query_coverage(hsps) / query_len

                starts = [min(h.sbjct_start, h.sbjct_end) for h in hsps]
                ends = [max(h.sbjct_start, h.sbjct_end) for h in hsps]
                start = min(starts)
                end = max(ends)

                frames = []
                for h in hsps:
                    f = h.frame
                    if isinstance(f, tuple):
                        frames.append(f[1])
                    elif f is not None:
                        frames.append(f)

                strand = 1 if frames and Counter(frames).most_common(1)[0][0] > 0 else 2

                index = alignment.hit_id
                acc = extract_accession(index)

                indices.append(index)
                records.append({
                    "QC": round(qc, 2),
                    "Acc": acc,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "Length": end-start+1,
                })

    return pd.DataFrame(records, index=indices)


def fix_anomalous_hits(df: pd.DataFrame, xml_path: str, max_len: int = 10000) -> pd.DataFrame:
    updated = df.copy()
    long_hits = df[(df["End"] - df["Start"]) > max_len]
    if long_hits.empty:
        return updated

    with open(xml_path) as handle:
        blast_records = NCBIXML.parse(handle)
        for record in blast_records:
            query_len = record.query_length

            for alignment in record.alignments:
                index = alignment.hit_id
                if index not in long_hits.index:
                    continue

                hsps = alignment.hsps
                if not hsps:
                    continue

                all_coords = [min(h.sbjct_start, h.sbjct_end) for h in hsps] + \
                             [max(h.sbjct_start, h.sbjct_end) for h in hsps]
                median_center = sorted(all_coords)[len(all_coords) // 2]

                filtered_hsps = [
                    h for h in hsps
                    if abs(h.sbjct_start - median_center) <= max_len // 2 and
                       abs(h.sbjct_end - median_center) <= max_len // 2
                ]

                if not filtered_hsps:
                    continue

                qc = calculate_query_coverage(filtered_hsps) / query_len
                starts = [min(h.sbjct_start, h.sbjct_end) for h in filtered_hsps]
                ends = [max(h.sbjct_start, h.sbjct_end) for h in filtered_hsps]
                start = min(starts)
                end = max(ends)

                frames = []
                for h in filtered_hsps:
                    f = h.frame
                    if isinstance(f, tuple):
                        frames.append(f[1])
                    elif f is not None:
                        frames.append(f)
                strand = 1 if frames and Counter(frames).most_common(1)[0][0] > 0 else 2

                updated.loc[index, "QC"] = round(qc, 2)
                updated.loc[index, "Start"] = start
                updated.loc[index, "End"] = end
                updated.loc[index, "Strand"] = strand
                updated.loc[index, "Acc"] = extract_accession(index)
                updated.loc[index, "Length"] = end-start+1

    return updated


def filter_hits(df: pd.DataFrame, min_qc: float = 0.0, max_len: int = 10000) -> pd.DataFrame:
    return df[(df["QC"] >= min_qc) & ((df["End"] - df["Start"]) <= max_len)]

