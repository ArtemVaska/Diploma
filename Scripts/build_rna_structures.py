import subprocess
from pathlib import Path
import argparse
import sys


def read_fasta(file_path: Path) -> dict:
    sequences = {}
    current_id = None
    with open(file_path) as infile:
        for line in infile:
            if line.startswith(">"):
                current_id = line[1:].strip()
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.strip()
    return sequences


def find_intron_coordinates(full_seq: str, intron_seq: str) -> tuple[int, int] | None:
    start = full_seq.find(intron_seq)
    if start == -1:
        return None
    return (start + 1, start + len(intron_seq))  # 1-based inclusive


def run_rnafold_with_highlight(fasta_path: Path, paint_path: Path | None = None) -> None:
    if not fasta_path.exists():
        print(f"❌ Error: input file {fasta_path} not found.")
        sys.exit(1)

    input_sequences = read_fasta(fasta_path)
    intron_sequences = read_fasta(paint_path) if paint_path else {}

    project_root = Path(__file__).resolve().parents[1]
    output_root = project_root / "rnafold_output"
    output_root.mkdir(exist_ok=True)

    for seq_id, full_seq in input_sequences.items():
        print(f"→ Processing {seq_id}...")
        organism_dir = output_root / seq_id
        organism_dir.mkdir(parents=True, exist_ok=True)

        # Run RNAfold without .ps output
        result = subprocess.run(
            ["RNAfold", "--noLP", "--noPS"],
            input=full_seq.strip(),
            capture_output=True,
            text=True
        )
        lines = result.stdout.strip().splitlines()
        if len(lines) < 2:
            print(f"⚠️ Error: RNAfold failed for {seq_id}")
            continue

        structure_line = lines[1]  # dot-bracket + energy
        structure, energy = structure_line.rsplit(' ', 1)

        fold_file = organism_dir / f"{seq_id}.fold"
        with open(fold_file, "w") as f:
            f.write(f"{full_seq}\n{structure} {energy}\n")

        # Generate dot plot (rna.ps will be created here)
        subprocess.run(
            ["RNAfold", "-p", "-d2", "--noLP"],
            input=full_seq.strip(),
            cwd=organism_dir,
            text=True
        )
        rna_ps = organism_dir / "rna.ps"
        if rna_ps.exists():
            rna_ps.rename(organism_dir / f"{seq_id}_dotplot.ps")

        dot_ps = organism_dir / "dot.ps"
        if dot_ps.exists():
            dot_ps.rename(organism_dir / f"{seq_id}_dot.ps")

        # Highlight intron using RNAplot
        if seq_id in intron_sequences:
            coords = find_intron_coordinates(full_seq, intron_sequences[seq_id])
            if coords:
                start, end = coords
                highlight_cmd = ["--pre", f"{start} {end} 8 GREEN omark"]
                with open(fold_file) as infile:
                    subprocess.run(
                        ["RNAplot"] + highlight_cmd,
                        stdin=infile,
                        cwd=organism_dir,
                    )
                highlight_eps = organism_dir / "rna.eps"
                if highlight_eps.exists():
                    highlight_eps.rename(organism_dir / f"{seq_id}_highlight.eps")
                else:
                    print(f"⚠️ Highlight file not created for {seq_id}")
            else:
                print(f"⚠️ Warning: intron not found in {seq_id}. Skipping highlight.")

        print(f"✔ Done: {seq_id}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RNAfold with intron highlighting")
    parser.add_argument("-input", "--input", type=str, required=True, help="FASTA with full exon-intron sequences")
    parser.add_argument("--paint", type=str, required=False, help="FASTA with intron sequences to highlight")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_path = Path.cwd() / args.input
    paint_path = Path.cwd() / args.paint if args.paint else None
    run_rnafold_with_highlight(input_path, paint_path)
