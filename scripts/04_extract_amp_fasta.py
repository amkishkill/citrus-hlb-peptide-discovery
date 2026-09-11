"""
Step 3c: Extract the actual peptide sequences for AMPDiscover-classified AMPs.

AMPDiscover's output CSV only contains sequence IDs and prediction columns
-- not the sequences themselves. This script joins those IDs back against
the cleaned FASTA (from 01_filter_dbsmorf.py) to produce a FASTA of just
the predicted AMPs, ready for MMseqs2 clustering.

Usage:
    python scripts/04_extract_amp_fasta.py \
        --clean-fasta data/dbsmorf_hmp_gut_clean.fasta \
        --amp-csv data/ampdiscover_active_amps.csv \
        --id-col IDs \
        --out data/dbsmorf_amps.fasta
"""

import argparse
from pathlib import Path

import pandas as pd


def read_fasta_dict(path: Path) -> dict:
    """Return {header_id: sequence} where header_id is the '>' line with
    the leading '>' stripped."""
    seqs = {}
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        seqs[header] = "".join(seq_lines)
    return seqs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean-fasta", required=True)
    parser.add_argument("--amp-csv", required=True)
    parser.add_argument(
        "--id-col",
        default="IDs",
        help="Column in the AMPDiscover CSV containing the sequence ID "
        "(default: 'IDs', matching AMPDiscover's own output header)",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    print(f"Loading cleaned FASTA from {args.clean_fasta} ...")
    seq_dict = read_fasta_dict(Path(args.clean_fasta))
    print(f"  {len(seq_dict)} sequences available")

    print(f"Loading AMP predictions from {args.amp_csv} ...")
    df = pd.read_csv(args.amp_csv)
    print(f"  {len(df)} predicted-AMP rows")

    ids = df[args.id_col].astype(str).tolist()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_missing = 0
    with open(out_path, "w") as fh:
        for seq_id in ids:
            seq = seq_dict.get(seq_id)
            if seq is None:
                n_missing += 1
                continue
            fh.write(f">{seq_id}\n{seq}\n")
            n_written += 1

    print(f"\nWrote {n_written} sequences to {out_path}")
    if n_missing:
        print(
            f"WARNING: {n_missing} IDs from the AMP CSV were not found in "
            f"the cleaned FASTA -- check that --id-col matches the FASTA "
            f"header naming convention."
        )


if __name__ == "__main__":
    main()