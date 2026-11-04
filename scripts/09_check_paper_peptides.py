#!/usr/bin/env python3
"""Step 5d (2/2): trace the paper's experimentally tested peptides (APP1-APP22) through
every stage of our pipeline, matching by SEQUENCE.

Input table: data/paper_tested_peptides.tsv from 09_build_paper_peptide_table.py
(needs a 'sequence' column). The FASTA headers use smorf_N identifiers, so peptides are
matched by sequence, not by ID.

Matching:
  exact   the identical sequence is in the cleaned sORF FASTA.
  family  no identical sequence, but one or more sORFs are at least --min-identity
          (default 0.80) identical. Such peptides come in near-identical families (the
          APP3/APP9 family has ~550 ORFs but only 17 distinct sequences), and MMseqs2
          clusters them together, so the family is traced instead. The member that got
          furthest in the pipeline is reported, with its identity to the paper's peptide.

Stages (a missing optional file is skipped with a warning):
  1. clean sORFs (--clean-fasta)   2. AMPs (--amp-fasta)   3. representatives (--rep-fasta)
  4. any BLAST hit (--all-hits)    5. passed identity/coverage (--annotated)
  6. final template (--final)

A sORF that is not itself a representative is followed through the representative of its
MMseqs2 cluster if --cluster-tsv (columns: representative, member) is available;
the output says "via representative X".

Usage (defaults match the repo layout):
  python scripts/10_check_paper_peptides.py
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

STAGES = ["clean sORFs", "AMPs", "representatives", "BLAST hit",
          "identity/coverage", "final template"]
KMER = 5


def read_fasta(path):
    """-> {sequence: [record ids]}, {record id: sequence}"""
    by_seq, by_id, cur = {}, {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                words = line[1:].split()
                cur = words[0] if words else None
                if cur is not None:
                    by_id[cur] = ""
            elif cur is not None:
                by_id[cur] += line.strip().upper()
    for rid, seq in by_id.items():
        by_seq.setdefault(seq, []).append(rid)
    return by_seq, by_id


def identity(a, b):
    """Global alignment identity: matches / alignment length (match +2, mismatch -1, gap -2)."""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = -2 * i
    for j in range(1, m + 1):
        score[0][j] = -2 * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = score[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1)
            score[i][j] = max(d, score[i - 1][j] - 2, score[i][j - 1] - 2)
    i, j, matches, cols = n, m, 0, 0
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1):
            matches += a[i - 1] == b[j - 1]
            i, j = i - 1, j - 1
        elif i > 0 and score[i][j] == score[i - 1][j] - 2:
            i -= 1
        else:
            j -= 1
        cols += 1
    return matches / cols if cols else 0.0


def maybe(path, label):
    if path and Path(path).exists():
        return path
    print(f"WARNING: {label} not found ({path}); stage skipped", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paper-list", default="data/papers/paper_tested_peptides.tsv")
    ap.add_argument("--clean-fasta", default="data/dbsmorf_hmp_gut_clean.fasta")
    ap.add_argument("--amp-fasta", default="data/dbsmorf_amps.fasta")
    ap.add_argument("--rep-fasta", default="data/mmseqs_out_rep_seq.fasta")
    ap.add_argument("--cluster-tsv", default="data/mmseqs_out_cluster.tsv")
    ap.add_argument("--all-hits", default="data/blast_out/all_hits.tsv")
    ap.add_argument("--annotated", default="data/blast_filtered/templates_annotated.tsv")
    ap.add_argument("--final", default="data/blast_filtered/templates_filtered.tsv")
    ap.add_argument("--min-identity", type=float, default=0.80,
                    help="identity (0-1) for family matching when no exact match exists (default 0.80)")
    ap.add_argument("--out", default="data/blast_filtered/paper_peptides_check.tsv")
    args = ap.parse_args()

    paper = pd.read_csv(args.paper_list, sep="\t", dtype=str).fillna("")
    if "sequence" not in paper.columns:
        sys.exit("ERROR: --paper-list needs a 'sequence' column (run 09_build_paper_peptide_table.py)")

    fasta = {}
    for key, path, label in [("clean", args.clean_fasta, "clean FASTA"),
                             ("amp", args.amp_fasta, "AMP FASTA"),
                             ("rep", args.rep_fasta, "representative FASTA")]:
        p = maybe(path, label)
        fasta[key] = read_fasta(p) if p else None

    member2rep = {}
    p = maybe(args.cluster_tsv, "cluster TSV")
    if p:
        cl = pd.read_csv(p, sep="\t", header=None, names=["rep", "member"], dtype=str)
        member2rep = dict(zip(cl.member, cl.rep))

    p = maybe(args.all_hits, "all BLAST hits")
    hit_ids = set(pd.read_csv(p, sep="\t", header=None, usecols=[0], dtype=str)[0]) if p else None
    p = maybe(args.annotated, "annotated templates")
    ann_ids = set(pd.read_csv(p, sep="\t", dtype=str).qaccver) if p else None
    p = maybe(args.final, "final templates")
    final = pd.read_csv(p, sep="\t") if p else None
    final_ids = set(final.qaccver) if final is not None else None

    # k-mer index of the cleaned sORFs, for family matching
    kindex = defaultdict(set)
    if fasta["clean"] is not None:
        for rid, s in fasta["clean"][1].items():
            for i in range(len(s) - KMER + 1):
                kindex[s[i:i + KMER]].add(rid)

    def find(kind, seq):
        f = fasta[kind]
        return None if f is None or not seq else (f[0].get(seq) or [None])[0]

    def known(flag, loaded):
        return flag if loaded else None  # None = file missing, stage cannot be judged

    def trace(seq):
        """Follow one sORF sequence through the stages -> (stage index, label, matched id, via, template)."""
        clean_id, amp_id, rep_id = find("clean", seq), find("amp", seq), find("rep", seq)
        matched = rep_id or amp_id or clean_id or ""
        proxy, via = rep_id, ""
        if not rep_id and amp_id:
            rep = member2rep.get(amp_id)
            if rep and rep != amp_id:
                proxy, via = rep, rep
        stage_ok = [
            known(clean_id is not None, fasta["clean"] is not None),
            known(amp_id is not None, fasta["amp"] is not None),
            known(proxy is not None, fasta["rep"] is not None),
            known(proxy in hit_ids if hit_ids is not None else False, hit_ids is not None),
            known(proxy in ann_ids if ann_ids is not None else False, ann_ids is not None),
            known(proxy in final_ids if final_ids is not None else False, final_ids is not None),
        ]
        idx = -1
        for k, ok in enumerate(stage_ok):
            if ok is False:
                break
            if ok:
                idx = k
        label = STAGES[idx] if idx >= 0 else "not found in any file"
        if via and idx >= 2:
            label += f" (via representative {via})"
        tmpl = ""
        if final is not None and proxy:
            best = final[final.qaccver == proxy].head(1)
            if len(best):
                b = best.iloc[0]
                res = "NA" if pd.isna(b.get("resolution")) else f"{b['resolution']:.2f} A"
                tmpl = f"{b.pdb_id}_{b.chain} ({b.method}, {res}, {b.pident:.0f}% id, {b.qcovs:.0f}% cov)"
        return idx, label, matched, via, tmpl

    rows = []
    for _, r in paper.iterrows():
        seq = r.sequence.strip().upper()
        base = {"name": r["name"], "peptide_id": r.peptide_id, "active_in_paper": r.active_in_paper,
                "sequence_source": r.get("sequence_source", ""), "length": len(seq) if seq else "",
                "sequence": seq, "match_type": "", "match_identity": "", "family_size": "",
                "matched_id": "", "represented_by": "", "best_template": ""}
        if not seq:
            rows.append({**base, "last_stage_reached": "no sequence available"})
            continue

        if find("clean", seq) is not None or fasta["clean"] is None:
            idx, label, matched, via, tmpl = trace(seq)
            rows.append({**base, "match_type": "exact", "match_identity": 1.0, "family_size": 1,
                         "last_stage_reached": label, "matched_id": matched, "represented_by": via,
                         "best_template": tmpl})
            continue

        # family matching: candidate sORFs sharing a k-mer, then global identity
        cand = set()
        for i in range(len(seq) - KMER + 1):
            cand |= kindex.get(seq[i:i + KMER], set())
        members = {}  # distinct sequences >= threshold
        for rid in cand:
            s = fasta["clean"][1][rid]
            if s in members:
                continue
            ident = identity(seq, s)
            if ident >= args.min_identity:
                members[s] = ident
        if not members:
            rows.append({**base, "match_type": "none",
                         "last_stage_reached": f"not in sORF set (no sORF >= {args.min_identity:.0%} identical)"})
            continue
        best = None
        for s, ident in members.items():
            result = trace(s)
            key = (result[0], ident)
            if best is None or key > best[0]:
                best = (key, s, ident, result)
        _, s, ident, (idx, label, matched, via, tmpl) = best
        n_orfs = sum(len(fasta["clean"][0][x]) for x in members)
        rows.append({**base, "match_type": "family", "match_identity": round(ident, 3),
                     "family_size": f"{len(members)} sequences / {n_orfs} sORFs",
                     "last_stage_reached": label, "matched_id": matched, "represented_by": via,
                     "best_template": tmpl})

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, sep="\t", index=False)
    show = res[["name", "active_in_paper", "match_type", "match_identity", "matched_id",
                "last_stage_reached", "best_template"]]
    print(show.to_string(index=False))
    print("\nLast stage reached:")
    print(res.last_stage_reached.str.replace(r" \(via.*\)", "", regex=True).value_counts().to_string())
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()