#!/usr/bin/env bash
# Step 5a: build / download the protein BLAST database of PDB sequences.
#
# Two sources are supported:
#   pdbaa   (default) NCBI's prebuilt, non-redundant PDB protein database.
#           Identical chains are merged; the hit title lists the other PDB IDs.
#   seqres  RCSB pdb_seqres.txt, protein chains only, every chain kept.
#           Useful if one wants to choose templates by resolution/chain manually.
#
# Usage:
#   bash scripts/06_setup_blast_db.sh [pdbaa|seqres] [DB_DIR]
#
# Requires BLAST+ (conda install -c bioconda blast) on PATH.
# The database directory is NOT tracked by git (see .gitignore).

set -euo pipefail

SOURCE="${1:-pdbaa}"
DB_DIR="${2:-data/blastdb}"

command -v makeblastdb >/dev/null || { echo "ERROR: BLAST+ not found on PATH" >&2; exit 1; }
mkdir -p "$DB_DIR"
cd "$DB_DIR"

case "$SOURCE" in
  pdbaa)
    command -v update_blastdb.pl >/dev/null || { echo "ERROR: update_blastdb.pl not found (ships with BLAST+)" >&2; exit 1; }
    update_blastdb.pl --decompress pdbaa
    echo "Done. Use:  -db $(pwd)/pdbaa   (or export BLASTDB=$(pwd))"
    ;;
  seqres)
    URL="https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
    [ -s pdb_seqres.txt ] || { curl -fsSL "$URL" -o pdb_seqres.txt.gz && gunzip -f pdb_seqres.txt.gz; }
    # keep protein chains only (drops DNA/RNA entries, which have mol:na)
    awk '/^>/{p=($0 ~ /mol:protein/)} p' pdb_seqres.txt > pdb_protein.fasta
    makeblastdb -in pdb_protein.fasta -dbtype prot -parse_seqids \
                -out pdb_seqres -title "PDB seqres (protein)"
    echo "Done. Use:  -db $(pwd)/pdb_seqres"
    ;;
  *)
    echo "Unknown source '$SOURCE' (use pdbaa or seqres)" >&2
    exit 1
    ;;
esac
