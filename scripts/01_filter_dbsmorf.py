"""
Step 1-2: Clean the raw DBsmORF (HMP, Gut) download.

Reproduces the filtering described in the Methods of Zhao et al. 2025
(Science 388, 191-198).

This script:
  1. Loads the raw TSV export from the DBsmORF Download tab
     (Source Database = HMP, Body Site = Gut, default SmORFinder
     significance filters, Protein Sequence field included).
  2. Removes exact-duplicate protein sequences.
  3. Removes any sequence containing a character outside the 20
     standard amino acids.
  4. Writes a cleaned FASTA (for AMPDiscover / SmORFinder downstream
     steps) and a cleaned TSV (retaining the original metadata columns).

Usage:
    python 01_filter_dbsmorf.py \
        --input data/dbsmorf_hmp_gut_raw.tsv \
        --seq-col protein_seq \
        --id-col smorf_id \
        --out-prefix data/dbsmorf_hmp_gut_clean
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def has_only_standard_aa(seq: str) -> bool:
    """Return True if every character in seq is a standard amino acid."""
    if not isinstance(seq, str) or len(seq) == 0:
        return False
    return set(seq.upper()) <= STANDARD_AA


def find_seq_column(df: pd.DataFrame, requested: str | None) -> str:
    """Locate the protein sequence column, allowing for naming drift
    in DBsmORF's export (e.g. 'protein_seq', 'protein_sequence')."""
    if requested and requested in df.columns:
        return requested
    candidates = [c for c in df.columns if re.search(r"protein.*seq", c, re.I)]
    if not candidates:
        raise ValueError(
            f"Could not find a protein sequence column. "
            f"Available columns: {list(df.columns)}"
        )
    return candidates[0]


def find_id_column(df: pd.DataFrame, requested: str | None) -> str:
    if requested and requested in df.columns:
        return requested
    candidates = [c for c in df.columns if re.search(r"(smorf.*id|^id$)", c, re.I)]
    if candidates:
        return candidates[0]
    return None  # fall back to a generated index-based id


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Raw DBsmORF TSV export")
    parser.add_argument(
        "--seq-col",
        default=None,
        help="Name of the protein sequence column (auto-detected if omitted)",
    )
    parser.add_argument(
        "--id-col",
        default=None,
        help="Name of the smORF id column (auto-detected if omitted)",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output path prefix; writes <prefix>.fasta and <prefix>.tsv",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        sys.exit(f"Input file not found: {in_path}")

    print(f"Loading {in_path} ...")
    df = pd.read_csv(in_path, sep="\t")
    n_start = len(df)
    print(f"  {n_start} rows loaded")

    seq_col = find_seq_column(df, args.seq_col)
    id_col = find_id_column(df, args.id_col)
    print(f"  Using sequence column: '{seq_col}'")
    print(f"  Using id column: '{id_col if id_col else '(generated)'}'")

    if id_col is None:
        df["_generated_id"] = [f"smorf_{i}" for i in range(len(df))]
        id_col = "_generated_id"

    # Drop rows with missing/empty sequences up front
    df = df[df[seq_col].notna() & (df[seq_col].astype(str).str.len() > 0)].copy()
    n_after_missing = len(df)
    print(f"  Dropped {n_start - n_after_missing} rows with missing sequences")

    # Remove exact-duplicate sequences (keep first occurrence)
    df = df.drop_duplicates(subset=[seq_col], keep="first")
    n_after_dedup = len(df)
    print(f"  Dropped {n_after_missing - n_after_dedup} duplicate sequences")

    # Remove sequences containing non-standard amino acids
    mask_standard = df[seq_col].apply(has_only_standard_aa)
    df_clean = df[mask_standard].copy()
    n_final = len(df_clean)
    print(
        f"  Dropped {n_after_dedup - n_final} sequences with non-standard "
        f"amino acids"
    )
    print(f"Final count: {n_final} sequences")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta_path = out_prefix.with_suffix(".fasta")
    with open(fasta_path, "w") as fh:
        for _, row in df_clean.iterrows():
            fh.write(f">{row[id_col]}\n{row[seq_col]}\n")
    print(f"Wrote {fasta_path}")

    tsv_path = out_prefix.with_suffix(".tsv")
    df_clean.to_csv(tsv_path, sep="\t", index=False)
    print(f"Wrote {tsv_path}")


if __name__ == "__main__":
    main()
