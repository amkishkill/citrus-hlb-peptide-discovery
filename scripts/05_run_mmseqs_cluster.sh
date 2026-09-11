#!/usr/bin/env bash
# Step 4: Cluster AMPs with MMseqs2 and select representative sequences.
# Usage:
#   bash scripts/05_run_mmseqs_cluster.sh \
#       data/dbsmorf_amps.fasta \
#       data/mmseqs_out
#
# Outputs (written into <out_prefix>_*):
#   <out_prefix>_rep_seq.fasta   -- one representative sequence per cluster
#                                    (this is your "representative AMPs" set)
#   <out_prefix>_cluster.tsv     -- cluster membership (rep_id <TAB> member_id)
#   <out_prefix>_all_seqs.fasta  -- all sequences grouped by cluster

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_fasta> <output_prefix> [min_seq_id] [coverage]"
    echo "  min_seq_id default: 0.5  (50% amino acid identity, per the paper)"
    echo "  coverage default:   0.8  (80% coverage, per the paper)"
    exit 1
fi

INPUT_FASTA="$1"
OUT_PREFIX="$2"
MIN_SEQ_ID="${3:-0.5}"
COVERAGE="${4:-0.8}"
TMP_DIR="$(dirname "$OUT_PREFIX")/mmseqs_tmp"

if ! command -v mmseqs &> /dev/null; then
    echo "ERROR: 'mmseqs' command not found on PATH."
    echo "Install it first -- see comments at the top of this script."
    exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")" "$TMP_DIR"

echo "Running: mmseqs easy-cluster $INPUT_FASTA $OUT_PREFIX $TMP_DIR"
echo "  --min-seq-id $MIN_SEQ_ID -c $COVERAGE --cov-mode 0"
echo ""

mmseqs easy-cluster \
    "$INPUT_FASTA" \
    "$OUT_PREFIX" \
    "$TMP_DIR" \
    --min-seq-id "$MIN_SEQ_ID" \
    -c "$COVERAGE" \
    --cov-mode 0

echo ""
echo "Done. Key outputs:"
echo "  ${OUT_PREFIX}_rep_seq.fasta  <- representative AMPs"
echo "  ${OUT_PREFIX}_cluster.tsv    <- cluster membership table"
echo ""
N_REPS=$(grep -c "^>" "${OUT_PREFIX}_rep_seq.fasta" || echo "0")
echo "Representative sequence count: $N_REPS"
