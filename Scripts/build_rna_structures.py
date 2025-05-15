import subprocess
from pathlib import Path
import argparse
import sys


def read_fasta(file_path: Path) -> dict:
    """
    Parses a FASTA file into a dictionary.

    :param file_path: Path to the FASTA file
    :return: Dictionary {header: sequence}
    """
    sequences = {}
    current_id = None
    with open(file_path) as f:
        for line in f:
            if line.startswith(">"):
                current_id = line[1:].strip().replace(" ", "_")
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.strip()
    return sequences


def run_rnafold(fasta_path: Path) -> None:
    """
    Runs RNAfold on sequences in the provided FASTA file.

    :param fasta_path: Path to the input FASTA file (relative to project root)
    :return: None
    """
    if not fasta_path.exists():
        print(f"❌ Error: file {fasta_path} not found.")
        sys.exit(1)

    project_root: Path = Path(__file__).resolve().parents[1]
    output_root: Path = project_root / "rnafold_output"
    output_root.mkdir(exist_ok=True)

    input_sequences = read_fasta(fasta_path)

    for seq_id, sequence in input_sequences.items():
        print(f"→ Processing {seq_id}...")

        organism_dir: Path = output_root / seq_id
        organism_dir.mkdir(parents=True, exist_ok=True)

        # Run RNAfold to get structure only, suppressing rna.ps
        result = subprocess.run(
            ["RNAfold", "--noLP", "--noPS"],
            input=sequence.strip(),
            capture_output=True,
            text=True
        )

        fold_file: Path = organism_dir / f"{seq_id}.fold.txt"
        with open(fold_file, "w") as f:
            f.write(result.stdout)

        # Generate dot plot (rna.ps will be created here)
        subprocess.run(
            ["RNAfold", "-p", "-d2", "--noLP"],
            input=sequence.strip(),
            cwd=organism_dir,
            text=True
        )
        dotplot_src: Path = organism_dir / "rna.ps"
        if dotplot_src.exists():
            dotplot_dst: Path = organism_dir / f"{seq_id}_dotplot.ps"
            dotplot_src.rename(dotplot_dst)

        print(f"✔ Done: {seq_id}")


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments.

    :return: argparse.Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(description="Run RNAfold to generate RNA secondary structures")
    parser.add_argument(
        "-input", "--input", type=str, required=True,
        help="Path to input FASTA file relative to the project root"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args: argparse.Namespace = parse_args()
    fasta_full_path: Path = Path.cwd() / args.input
    run_rnafold(fasta_full_path)
