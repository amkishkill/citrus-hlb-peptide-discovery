"""
Step 3a: Split the cleaned DBsmORF FASTA into chunks under AMPDiscover's
1 MB per-file upload limit.

AMPDiscover (https://biocom-ampdiscover.cicese.mx) enforces a 1 MB max
file size per submission. Our cleaned FASTA (122,826 sequences, ~7.1 MB)
must be split into multiple batches, submitted individually through the
web interface, and the results merged afterward (see
03_merge_ampdiscover_results.py).

Usage:
    python scripts/02_split_fasta_for_ampdiscover.py \
        --input data/dbsmorf_hmp_gut_clean.fasta \
        --out-dir data/ampdiscover_batches \
        --max-bytes 900000
"""

import argparse
from pathlib import Path


def read_fasta_records(path: Path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def record_size(header: str, seq: str) -> int:
    """Size in bytes of a FASTA record as written to disk."""
    return len(header) + 1 + len(seq) + 1  # +1 for each newline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Cleaned FASTA file")
    parser.add_argument(
        "--out-dir", required=True, help="Directory to write batch_NNN.fasta files"
    )
    parser.add_argument(
        "--max-bytes",
        type=int,
        default=900_000,
        help="Max size per batch file in bytes (default 900000, safely under "
        "AMPDiscover's 1 MB limit)",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    batch_idx = 1
    current_size = 0
    current_records = []
    manifest_rows = []
    total_records = 0

    def flush_batch():
        nonlocal batch_idx, current_size, current_records
        if not current_records:
            return
        batch_path = out_dir / f"batch_{batch_idx:03d}.fasta"
        with open(batch_path, "w") as fh:
            for header, seq in current_records:
                fh.write(f"{header}\n{seq}\n")
        manifest_rows.append(
            (batch_path.name, len(current_records), current_size)
        )
        print(
            f"  Wrote {batch_path} "
            f"({len(current_records)} sequences, {current_size:,} bytes)"
        )
        batch_idx += 1
        current_size = 0
        current_records = []

    print(f"Reading {in_path} ...")
    for header, seq in read_fasta_records(in_path):
        total_records += 1
        size = record_size(header, seq)

        if size > args.max_bytes:
            # A single sequence alone exceeds the limit -- extremely
            # unlikely for short peptides, but guard against it anyway.
            print(
                f"  WARNING: record '{header}' alone exceeds --max-bytes "
                f"({size} bytes); writing it to its own batch."
            )
            flush_batch()
            current_records = [(header, seq)]
            current_size = size
            flush_batch()
            continue

        if current_size + size > args.max_bytes:
            flush_batch()

        current_records.append((header, seq))
        current_size += size

    flush_batch()

    manifest_path = out_dir / "manifest.tsv"
    with open(manifest_path, "w") as fh:
        fh.write("batch_file\tn_sequences\tsize_bytes\n")
        for row in manifest_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"\nTotal input sequences: {total_records}")
    print(f"Wrote {len(manifest_rows)} batch files to {out_dir}")
    print(f"Manifest: {manifest_path}")
    print(
        "\nNext: submit each batch_NNN.fasta to AMPDiscover "
        "(https://biocom-ampdiscover.cicese.mx) using the ProtDCal-AMP_RF "
        "model, download each result CSV, then run "
        "02_merge_ampdiscover_results.py"
    )


if __name__ == "__main__":
    main()