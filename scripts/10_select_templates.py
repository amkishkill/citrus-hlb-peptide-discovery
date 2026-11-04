#!/usr/bin/env python3
"""Step 6a: choose ONE template per peptide and write the input files for Modeller.

Reads the hits that passed the filters in 08_filter_blast_hits.py
(data/blast_filtered/templates_filtered.tsv), ranks the candidate templates of each
peptide, keeps the best one and writes everything a Modeller batch run needs.

Ranking (--rank):
  identity  (default) 1. fraction of the peptide's residues identical to the template
                         (matches / query length: rewards identity AND coverage together)
                      2. query coverage of the alignment (qcovhsp)
                      3. BLAST bitscore
                      4. resolution, lower first (entries without one, e.g. NMR, come last)
                      5. PDB ID, chain (only to make the result deterministic)
  coverage            qcovhsp first, then % identity of the alignment, then 3-5 above.
                      Favours templates that span the whole peptide even at low identity.
If one query has several alignments to the same chain, only the best one is kept.

Structures (--structures):
  skip      (default) do not touch structure files; pick rank 1 for every peptide.
  local     use the files already in --pdb-dir (<pdbid>.pdb, lower case). A template
            whose file is missing, or lacks the chain / any C-alpha atoms, is rejected and
            the next-best candidate is tried.
  download  like local, but missing files are fetched from RCSB first (PDB format only;
            entries that RCSB does not offer as .pdb are rejected like any other failure).
  Rerunning is cheap: existing structure files are never downloaded again.

Outputs (in --out-dir, default data/modeller):
  targets.fasta         sequences of the peptides that got a template
  targets.ali           the same sequences in PIR format (what Modeller reads)
  selected_templates.tsv  one row per peptide: chosen template and how it was chosen
                          (this is the job list for the Modeller batch run)
  candidates_ranked.tsv   every candidate with its rank (audit trail)
  pdb_ids.txt           unique PDB IDs of the chosen templates (to fetch them elsewhere)
  no_template.tsv       peptides left without a template, and why
  summary.txt           counts
  templates/            structure files (only with --structures local/download)

Usage:
  python scripts/11_select_templates.py
  python scripts/11_select_templates.py --structures download
"""
import argparse
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

RCSB_PDB_URL = "https://files.rcsb.org/download/{pid}.pdb"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
PIR_WIDTH = 75

NUM_COLS = ["pident", "length", "qstart", "qend", "sstart", "send", "bitscore",
            "qlen", "slen", "qcovs", "qcovhsp", "resolution"]
NEEDED = ["qaccver", "saccver", "pident", "length", "qstart", "qend", "sstart", "send",
          "bitscore", "qlen", "slen", "qcovhsp", "pdb_id", "chain", "method", "resolution"]


# ------------------------------------------------------------------ input ---
def read_fasta(path):
    """-> dict id -> upper-case sequence (insertion-ordered)."""
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                words = line[1:].split()
                cur = words[0] if words else None
                if cur is not None:
                    seqs[cur] = ""
            elif cur is not None:
                seqs[cur] += line.strip().upper()
    return seqs


def read_candidates(path):
    # keep_default_na=False: a chain called "NA" must not become a missing value
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in NEEDED if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {path} lacks columns {missing}. Is it the templates_filtered.tsv "
                 "written by 08_filter_blast_hits.py?")
    for c in NUM_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")  # "" / "nan" -> NaN
    return df


# ---------------------------------------------------------------- ranking ---
def rank_candidates(df, rank):
    """Add qident/rank columns; keep the best alignment per (query, template chain)."""
    df = df.copy()
    df["matches"] = (df.pident / 100.0 * df.length).round()
    df["qident"] = df.matches / df.qlen
    df["res_key"] = df.resolution.fillna(float("inf"))
    if rank == "identity":
        keys = [("qident", False), ("qcovhsp", False), ("bitscore", False)]
    else:
        keys = [("qcovhsp", False), ("pident", False), ("bitscore", False)]
    keys += [("res_key", True), ("pdb_id", True), ("chain", True)]
    cols, asc = [k for k, _ in keys], [a for _, a in keys]
    df = df.sort_values(["qaccver"] + cols, ascending=[True] + asc, kind="mergesort")
    df = df.drop_duplicates(["qaccver", "pdb_id", "chain"], keep="first")
    df["rank"] = df.groupby("qaccver").cumcount() + 1
    return df.drop(columns=["res_key"])


# ------------------------------------------------------------- structures ---
def pdb_path(pdb_dir, pid):
    return Path(pdb_dir) / f"{pid.lower()}.pdb"


def count_ca(path, chain):
    """C-alpha residues of `chain` in the first model of a PDB-format file."""
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and len(line) > 26 and line[21] == chain \
                    and line[12:16].strip() == "CA":
                seen.add(line[22:27])  # residue number + insertion code
    return len(seen)


def fetch_pdb(pid, dest, retries=3):
    """Download one entry in PDB format. Returns True on success. Atomic (tmp + rename)."""
    url = RCSB_PDB_URL.format(pid=pid.upper())
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = resp.read()
            if b"ATOM" not in data:
                return False
            tmp = Path(str(dest) + ".part")
            tmp.write_bytes(data)
            tmp.rename(dest)
            return True
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 410):  # no .pdb version of this entry
                return False
            last = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
        time.sleep(2 ** attempt)
    print(f"WARNING: could not download {pid}: {last}", file=sys.stderr)
    return False


def prepare_structures(pids, pdb_dir, download, workers):
    """Make sure the files for `pids` exist (download mode). Returns the set of usable ids."""
    Path(pdb_dir).mkdir(parents=True, exist_ok=True)
    todo = [p for p in pids if not pdb_path(pdb_dir, p).exists()]
    if todo and download:
        print(f"Downloading {len(todo)} structure(s) from RCSB ...")
        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(lambda p: fetch_pdb(p, pdb_path(pdb_dir, p)), todo))
    return {p for p in pids if pdb_path(pdb_dir, p).exists()}


def choose(ranked, structures, pdb_dir, workers, max_tries):
    """Pick one candidate per query. -> (chosen DataFrame, {query: reason} for the others)."""
    if structures == "skip":
        return ranked[ranked["rank"] == 1].copy(), {}

    chosen_idx, ca_counts, why = {}, {}, {}
    cache = {}  # (pdb_id, chain) -> C-alpha count, or None if unusable
    pending = {q: g for q, g in ranked.groupby("qaccver", sort=False)}
    for attempt in range(1, max_tries + 1):
        if not pending:
            break
        round_rows = {}
        for q in list(pending):
            hit = pending[q][pending[q]["rank"] == attempt]
            if len(hit):
                round_rows[q] = hit.iloc[0]
            else:  # candidates exhausted
                why[q] = f"all {len(pending[q])} candidate(s) had no usable structure file/chain"
                del pending[q]
        have = prepare_structures(sorted({r.pdb_id for r in round_rows.values()}),
                                  pdb_dir, structures == "download", workers)
        for q, r in round_rows.items():
            key = (r.pdb_id, r.chain)
            if key not in cache:
                cache[key] = (count_ca(pdb_path(pdb_dir, r.pdb_id), r.chain)
                              if r.pdb_id in have else None)
            if cache[key]:  # None (no file) or 0 (chain absent / no CA) -> try next candidate
                chosen_idx[q], ca_counts[q] = r.name, cache[key]
                del pending[q]
    for q in pending:
        why[q] = f"top {max_tries} candidates had no usable structure file/chain"
    out = ranked.loc[list(chosen_idx.values())].copy()
    out["template_ca_residues"] = [ca_counts[q] for q in out.qaccver]
    return out, why


# ----------------------------------------------------------------- output ---
def wrap(seq, width):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{wrap(s, 80)}\n")


def write_pir(path, seqs):
    """Modeller PIR format, sequence-only entries (the template side is aligned by align2d)."""
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">P1;{k}\nsequence:{k}:::::::0.00: 0.00\n{wrap(s + '*', PIR_WIDTH)}\n\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", default="data/blast_filtered/templates_filtered.tsv",
                    help="filtered hits from 08 (default data/blast_filtered/templates_filtered.tsv)")
    ap.add_argument("--queries", default="data/mmseqs_out_rep_seq.fasta",
                    help="peptide FASTA used for the BLAST search (default data/mmseqs_out_rep_seq.fasta)")
    ap.add_argument("--out-dir", default="data/modeller")
    ap.add_argument("--rank", choices=["identity", "coverage"], default="identity")
    ap.add_argument("--structures", choices=["skip", "local", "download"], default="skip")
    ap.add_argument("--pdb-dir", default=None, help="structure files (default <out-dir>/templates)")
    ap.add_argument("--workers", type=int, default=4, help="parallel downloads (default 4)")
    ap.add_argument("--max-tries", type=int, default=5,
                    help="with --structures local/download: candidates tried per peptide (default 5)")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.pdb_dir or str(out / "templates")

    seqs = read_fasta(args.queries)
    if not seqs:
        sys.exit("ERROR: no sequences found in --queries")

    raw = read_candidates(args.hits)
    n_rows = len(raw)
    bad = raw.pdb_id.eq("") | raw.chain.eq("")
    n_bad = int(bad.sum())
    raw = raw[~bad]
    unknown_q = sorted(set(raw.qaccver) - set(seqs))
    if unknown_q:
        sys.exit(f"ERROR: {len(unknown_q)} query IDs in --hits are not in --queries "
                 f"(e.g. {unknown_q[:3]}). Use the FASTA that was given to BLAST.")
    # BLAST qlen must agree with the FASTA, otherwise hits and sequences are out of sync
    flen = raw.qaccver.map(lambda q: len(seqs[q]))
    n_len_mismatch = int((flen != raw.qlen).sum())
    if n_len_mismatch:
        sys.exit(f"ERROR: query length differs between hits and FASTA for {n_len_mismatch} rows.")

    ranked = rank_candidates(raw, args.rank)
    ranked.to_csv(out / "candidates_ranked.tsv", sep="\t", index=False)

    chosen, why = choose(ranked, args.structures, pdb_dir, args.workers, args.max_tries)

    # peptides we cannot send to Modeller: non-standard residues
    nonstd = {q for q in chosen.qaccver
              if set(seqs[q]) - STANDARD_AA} if len(chosen) else set()
    if nonstd:
        chosen = chosen[~chosen.qaccver.isin(nonstd)]

    # ---- manifest ----
    if len(chosen):
        chosen = chosen.sort_values("qaccver", key=lambda s: s.map({q: i for i, q in enumerate(seqs)}))
        n_cand = ranked.groupby("qaccver").size()
        man = pd.DataFrame({
            "target_id": chosen.qaccver.values,
            "target_length": chosen.qlen.astype(int).values,
            "template_code": [f"{p.lower()}{c}" for p, c in zip(chosen.pdb_id, chosen.chain)],
            "pdb_id": chosen.pdb_id.str.lower().values,
            "chain": chosen.chain.values,
            "template_file": [pdb_path(pdb_dir, p).name for p in chosen.pdb_id],
            "method": chosen.method.values,
            "resolution": chosen.resolution.values,
            "pident": chosen.pident.values,
            "qident": chosen.qident.round(3).values,
            "qcovhsp": chosen.qcovhsp.astype(int).values,
            "bitscore": chosen.bitscore.values,
            "qstart": chosen.qstart.astype(int).values, "qend": chosen.qend.astype(int).values,
            "sstart": chosen.sstart.astype(int).values, "send": chosen.send.astype(int).values,
            "template_length": chosen.slen.astype(int).values,
            "n_candidates": chosen.qaccver.map(n_cand).values,
            "rank_used": chosen["rank"].values,
        })
        if "template_ca_residues" in chosen.columns:
            man["template_ca_residues"] = chosen.template_ca_residues.values
    else:
        man = pd.DataFrame(columns=["target_id"])
    man.to_csv(out / "selected_templates.tsv", sep="\t", index=False)

    final_seqs = {q: seqs[q] for q in man.target_id}
    write_fasta(out / "targets.fasta", final_seqs)
    write_pir(out / "targets.ali", final_seqs)
    (out / "pdb_ids.txt").write_text("".join(f"{p}\n" for p in sorted(set(man.get("pdb_id", [])))))

    # ---- peptides without a template ----
    reasons = []
    has_cand = set(ranked.qaccver)
    done = set(man.target_id)
    for q in seqs:
        if q in done:
            continue
        if q in nonstd:
            reasons.append((q, "non-standard residues in the peptide sequence"))
        elif q in why:
            reasons.append((q, why[q]))
        elif q not in has_cand:
            reasons.append((q, "no hit passed the 08 filters"))
    pd.DataFrame(reasons, columns=["target_id", "reason"]).to_csv(out / "no_template.tsv", sep="\t", index=False)

    n_fallback = int((man["rank_used"] > 1).sum()) if len(man) else 0
    summary = (
        f"Peptides in --queries:                  {len(seqs)}\n"
        f"Candidate rows read from --hits:        {n_rows}"
        + (f"  ({n_bad} dropped: no PDB ID/chain)" if n_bad else "") + "\n"
        f"Peptides with >=1 candidate:            {len(has_cand)}\n"
        f"Ranking:                                {args.rank}\n"
        f"Structure handling:                     {args.structures}\n"
        f"Peptides with a selected template:      {len(man)}\n"
        f"  of which not the rank-1 candidate:    {n_fallback}\n"
        f"Unique template structures needed:      {man.pdb_id.nunique() if len(man) else 0}\n"
        f"Peptides without a template:            {len(seqs) - len(man)}"
        f"  (see no_template.tsv)\n"
    )
    if args.structures == "skip":
        summary += ("NOTE: structure files were not checked. Rerun with --structures download to "
                    "fetch them and fall back to the next candidate when one is unusable.\n")
    (out / "summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
