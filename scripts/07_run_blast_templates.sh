#!/usr/bin/env bash
# Step 5b: batch blastp of the representative AMPs against the PDB database,
# to find homologous structures usable as homology-modelling templates.
#
# The query FASTA is split into chunks so the run is resumable: a chunk whose
# output file already exists is skipped. Outputs are written to a temp file and
# renamed only on success, so an interrupted chunk is never mistaken for done.
#
# Usage:
#   bash scripts/07_run_blast_templates.sh QUERY.fasta DB OUT_DIR
# Example:
#   bash scripts/07_run_blast_templates.sh data/mmseqs_out/<rep_seq>.fasta \
#        data/blastdb/pdbaa data/blast_out
#
# Tunable via environment variables:
#   THREADS=6  CHUNK_SIZE=1000  EVALUE=1000  MAX_TARGETS=20
#
# NOTE: the paper does not state the parameters of its peptide-vs-PDB search.
# These are our design choices (see README).

set -euo pipefail

QUERY="${1:?usage: $0 QUERY.fasta DB OUT_DIR}"
DB="${2:?usage: $0 QUERY.fasta DB OUT_DIR}"
OUT_DIR="${3:?usage: $0 QUERY.fasta DB OUT_DIR}"

THREADS="${THREADS:-6}"
CHUNK_SIZE="${CHUNK_SIZE:-1000}"
EVALUE="${EVALUE:-1000}"
MAX_TARGETS="${MAX_TARGETS:-20}"

# Keep in sync with COLS in 08_filter_blast_hits.py
OUTFMT="6 qaccver saccver pident length mismatch positive ppos gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp stitle"

command -v blastp >/dev/null || { echo "ERROR: blastp not found on PATH" >&2; exit 1; }
[ -s "$QUERY" ] || { echo "ERROR: query file '$QUERY' missing or empty" >&2; exit 1; }

# --- sanity checks on the query file -------------------------------------
n_seq=$(grep -c '^>' "$QUERY")
n_dup=$(grep '^>' "$QUERY" | cut -d' ' -f1 | sort | uniq -d | wc -l)
n_bad=$(grep -v '^>' "$QUERY" | grep -c '[^ACDEFGHIKLMNPQRSTVWY]' || true)
echo "Queries: $n_seq | duplicate IDs: $n_dup | lines with non-standard residues: $n_bad"
[ "$n_dup" -eq 0 ] || { echo "ERROR: duplicate sequence IDs in query file" >&2; exit 1; }
[ "$n_bad" -eq 0 ] || { echo "ERROR: non-standard residues in query file" >&2; exit 1; }

mkdir -p "$OUT_DIR/chunks" "$OUT_DIR/results"

# --- split into chunks (handles multi-line FASTA) -------------------------
rm -f "$OUT_DIR"/chunks/chunk_*.fasta
awk -v n="$CHUNK_SIZE" -v dir="$OUT_DIR/chunks" '
  /^>/ { if (c % n == 0) { if (f) close(f); f = sprintf("%s/chunk_%03d.fasta", dir, c / n) } c++ }
  { print > f }' "$QUERY"

# --- run blastp per chunk -------------------------------------------------
for chunk in "$OUT_DIR"/chunks/chunk_*.fasta; do
  name=$(basename "$chunk" .fasta)
  out="$OUT_DIR/results/$name.tsv"
  if [ -e "$out" ]; then echo "[skip] $name"; continue; fi
  echo "[run ] $name"
  blastp -task blastp-short \
    -query "$chunk" -db "$DB" -out "$out.tmp" \
    -evalue "$EVALUE" \
    -seg no -soft_masking false \
    -comp_based_stats 0 \
    -max_target_seqs "$MAX_TARGETS" \
    -num_threads "$THREADS" \
    -outfmt "$OUTFMT"
  mv "$out.tmp" "$out"
done

cat "$OUT_DIR"/results/chunk_*.tsv > "$OUT_DIR/all_hits.tsv"
echo "Done. $(wc -l < "$OUT_DIR/all_hits.tsv") hits written to $OUT_DIR/all_hits.tsv"
