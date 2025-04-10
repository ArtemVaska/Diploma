def calculate_n_fasta_lines(line: str, fasta_line_length) -> (int, int):
    """
    Additional function for plain_to_fasta function
    """
    line_length = len(line)
    n_fasta_lines = line_length // fasta_line_length
    if line_length % fasta_line_length == 0:
        extra_line = 0
    else:
        extra_line = line_length - (fasta_line_length * n_fasta_lines)
    return n_fasta_lines, extra_line


def plain_to_fasta(file: str,
                   fasta_line_length: int = 80,
                   uppercase: bool = False) -> None:
    """
    Converts plain file to .fa format
    """
    with open(file, "r") as infile:
        with open(file.replace("_plain", ""), "w") as outfile:
            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    header = line
                else:
                    if uppercase:
                        line = line.upper()
                    n_fasta_lines, extra_line = calculate_n_fasta_lines(line, fasta_line_length)

            outfile.write(header + "\n")
            if n_fasta_lines == 0:
                outfile.write(line + "\n")
            else:
                for n_fasta_line in range(n_fasta_lines):
                    outfile.write(line[n_fasta_line * fasta_line_length: (n_fasta_line + 1) * fasta_line_length] + "\n")
                outfile.write(line[-extra_line:] + "\n")


def read_fasta(file: str) -> dict:
    """
    Reads one- or multiline fasta-file

    :param file: path to file
    :return: dictionary id:seq
    """
    seqs = {}
    seq = []
    with open(file) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                seq.append(line)
    seqs[header] = "".join(seq)
    return seqs


def read_single_fasta(file: str) -> str:
    """
    Reads fasta-file with 1 sequence

    :param file: path to file
    :return: sequence
    """
    seq_list = []
    with open(file) as infile:
        for line in infile:
            line = line.strip()
            if not line.startswith(">"):
                seq_list.append(line)
    seq_str = "".join(seq_list)
    return seq_str


def dict_align_to_fasta(dict_align: dict, filename: str) -> None:
    """
    Creates multiline fasta file with seqs for next alignment.

    :param dict_align:
    :param filename:
    :return:
    """
    with open(filename, "w") as outfile:
        for header, seq in dict_align.items():
            outfile.write(">" + header + "\n")
            outfile.write(seq + "\n")


def exons_to_cds_plain(infilename: str,
                       outfilename: str,
                       header: str):
    """
    Converts exons.fa to CDS oneline fasta
    """
    seq = []
    with open(infilename, "r") as infile:
        with open(outfilename, "w") as outfile:
            outfile.write(">" + header + "\n")
            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    continue
                else:
                    seq.append(line)
            outfile.write("".join(seq))
