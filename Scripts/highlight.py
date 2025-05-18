from IPython.display import HTML

from fasta_processing import read_fasta, read_single_fasta
from data_processing import find_codon


def highlight_and_wrap(seq: str, highlight: str, width: int = 70) -> str:
    start = seq.find(highlight)
    if start == -1:
        raise ValueError("Интрон не найден в последовательности.")
    end = start + len(highlight)

    # Подсветка по позициям
    parts = []
    for i, nt in enumerate(seq):
        if i < 3:
            parts.append(f'<span style="color:green;">{nt}</span>')
        elif i >= len(seq) - 3:
            parts.append(f'<span style="color:red;">{nt}</span>')
        elif start <= i < end:
            parts.append(f'<span style="color:yellow;">{nt}</span>')
        else:
            parts.append(nt)

    highlighted_seq = ''.join(parts)

    # Разбиение по 70 символов с учётом HTML
    result_lines = []
    count = 0
    html_index = 0

    while count < len(seq):
        line_html = ""
        seen = 0
        while seen < width and html_index < len(highlighted_seq):
            if highlighted_seq[html_index] == '<':
                tag_start = html_index
                while highlighted_seq[html_index] != '>':
                    html_index += 1
                html_index += 1
                line_html += highlighted_seq[tag_start:html_index]
            else:
                line_html += highlighted_seq[html_index]
                html_index += 1
                seen += 1
        result_lines.append(line_html)
        count += width

    html_block = "<br>".join(result_lines)
    return f'<div style="white-space: pre-wrap; font-family: monospace;">{html_block}</div>'


def highlight_intron_in_seq(phylum: str, org_name: str):
    prefix = "../Datasets"
    postfix = "ncbi_dataset/data"

    cds = read_single_fasta(f"{prefix}/{phylum}/{org_name}/{postfix}/cds_cassette.fa")
    intron = read_fasta(f"{prefix}/{phylum}/{org_name}/{postfix}/cassette.fa")["cassette"]
    stop_codon_pos = find_codon(cds, which="stop")

    output = HTML(highlight_and_wrap(cds[:stop_codon_pos + 3], intron))

    return output
