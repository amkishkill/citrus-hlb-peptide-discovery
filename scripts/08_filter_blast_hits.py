#!/usr/bin/env python3
"""Step 5c: filter BLAST hits into usable homology-modelling templates.

Reads the combined outfmt-6 table produced by 07_run_blast_templates.sh and
  1. keeps hits passing identity / query-coverage cutoffs;
  2. annotates each template with its PDB experimental method and resolution
     (looked up from the RCSB Data API);
  3. applies the experimental filter (--method, --max-resolution);
  4. ranks the remaining templates per peptide.

Outputs (in --out-dir):
  templates_annotated.tsv  hits passing step 1, with method/resolution columns
                           (BEFORE the experimental filter, to inspect what was dropped)
  templates_filtered.tsv   hits passing all filters
  templates_best.tsv       top --top-n templates per peptide, best first
  peptides_with_template.txt / peptides_no_template.txt
  summary.txt              counts and the thresholds used

Usage:
  python scripts/08_filter_blast_hits.py --hits data/blast_out/all_hits.tsv \\
      --queries data/mmseqs_out_rep_seq.fasta --out-dir data/blast_filtered \\
      --min-identity 40 --min-qcov 70 --max-resolution 2.5

Experimental filter:
  --method          any (default), or a comma-separated list from
                    xray, nmr, em, neutron, multiple   (e.g. --method xray,nmr)
  --max-resolution  Angstrom cutoff (default 2.5); 'none' disables it.
                    The cutoff is applied ONLY to entries that have a numeric
                    resolution (X-ray, EM, ...). Entries without one (e.g. NMR)
                    are kept, unless --drop-no-resolution is given.
  Entries that RCSB does not know (obsolete IDs) are always dropped.

"""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

# Keep in sync with OUTFMT in 07_run_blast_templates.sh
COLS = ("qaccver saccver pident length mismatch positive ppos gapopen qstart qend "
        "sstart send evalue bitscore qlen slen qcovs qcovhsp stitle").split()

# Matches 1ABC_A, 1abc_A, pdb|1ABC|A
PDB_RE = re.compile(r"(?:pdb\|)?([0-9][A-Za-z0-9]{3})[_|]([A-Za-z0-9]+)")

RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
GQL_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info { experimental_method resolution_combined }
  }
}
"""
BATCH = 100
META_COLS = ["pdb_id", "method", "method_detail", "resolution"]
# --method choice -> value of rcsb_entry_info.experimental_method (lowercased)
METHODS = {"xray": "x-ray", "nmr": "nmr", "em": "em", "neutron": "neutron",
           "multiple": "multiple methods"}


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def split_pdb(acc):
    m = PDB_RE.match(str(acc))
    return (m.group(1).upper(), m.group(2)) if m else (None, None)


# ---------------------------------------------------------------- metadata ---
def fetch_batch(ids, retries=3):
    """Query RCSB Data API for a list of PDB IDs. Returns the 'entries' list."""
    payload = json.dumps({"query": GQL_QUERY, "variables": {"ids": ids}}).encode()
    req = urllib.request.Request(RCSB_GRAPHQL, data=payload,
                                 headers={"Content-Type": "application/json"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.load(resp)
            if not data.get("data"):
                raise RuntimeError(f"no data in response: {str(data)[:200]}")
            return data["data"]["entries"] or []
        except (urllib.error.URLError, TimeoutError, RuntimeError, json.JSONDecodeError) as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"RCSB request failed after {retries} attempts: {last}")


def parse_entry(entry):
    info = entry.get("rcsb_entry_info") or {}
    res = [r for r in (info.get("resolution_combined") or []) if r is not None]
    exptl = entry.get("exptl") or []
    return {
        "pdb_id": entry["rcsb_id"].upper(),
        "method": info.get("experimental_method") or "NA",
        "method_detail": "; ".join(e["method"] for e in exptl if e.get("method")) or "NA",
        "resolution": min(res) if res else float("nan"),
    }


def load_metadata(pdb_ids, cache_path, offline=False):
    """Return a DataFrame (META_COLS) for pdb_ids, fetching and caching any missing."""
    cache_path = Path(cache_path)
    if cache_path.exists():
        cache = pd.read_csv(cache_path, sep="\t")
    else:
        cache = pd.DataFrame(columns=META_COLS)
    missing = sorted(set(pdb_ids) - set(cache.pdb_id))
    if missing:
        if offline:
            sys.exit(f"ERROR: --offline set but {len(missing)} PDB IDs are missing from {cache_path}")
        print(f"Fetching experimental metadata for {len(missing)} PDB entries from RCSB ...")
        rows = []
        for i in range(0, len(missing), BATCH):
            chunk = missing[i:i + BATCH]
            got = {}
            for entry in fetch_batch(chunk):
                if entry:
                    row = parse_entry(entry)
                    got[row["pdb_id"]] = row
            for pid in chunk:  # unknown/obsolete IDs are cached as NOT_FOUND so we don't refetch
                rows.append(got.get(pid, {"pdb_id": pid, "method": "NOT_FOUND",
                                          "method_detail": "NOT_FOUND", "resolution": float("nan")}))
        cache = pd.concat([cache, pd.DataFrame(rows)], ignore_index=True)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(cache_path, sep="\t", index=False)
    return cache[cache.pdb_id.isin(pdb_ids)]


def parse_max_resolution(value):
    if str(value).strip().lower() in ("none", "off", "any"):
        return None
    return float(value)


def parse_methods(value):
    """'any' -> None; 'xray,nmr' -> {'x-ray', 'nmr'}"""
    names = [v.strip().lower().replace("-", "") for v in str(value).split(",") if v.strip()]
    if not names or "any" in names:
        return None
    bad = [n for n in names if n not in METHODS]
    if bad:
        raise argparse.ArgumentTypeError(f"unknown method(s) {bad}; choose from any, {', '.join(METHODS)}")
    return {METHODS[n] for n in names}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="combined BLAST outfmt-6 table")
    ap.add_argument("--queries", required=True, help="query FASTA (to report peptides with no hits)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-identity", type=float, default=30.0, help="min %% identity (default 30)")
    ap.add_argument("--min-qcov", type=float, default=70.0, help="min %% query coverage (default 70)")
    ap.add_argument("--max-evalue", type=float, default=None, help="optional E-value ceiling")
    ap.add_argument("--method", default="any", type=parse_methods,
                    help="'any' (default) or comma-separated: " + ", ".join(METHODS))
    ap.add_argument("--max-resolution", default="2.5", type=parse_max_resolution,
                    help="max resolution in Angstrom, applied only to entries that have one "
                         "(default 2.5); 'none' disables")
    ap.add_argument("--drop-no-resolution", action="store_true",
                    help="also drop entries without a numeric resolution (e.g. NMR) "
                         "when a --max-resolution cutoff is set")
    ap.add_argument("--pdb-meta", default="data/pdb_metadata.tsv",
                    help="cache file for PDB method/resolution (default data/pdb_metadata.tsv)")
    ap.add_argument("--offline", action="store_true", help="never access the network; use the cache only")
    ap.add_argument("--top-n", type=int, default=3, help="templates kept per peptide (default 3)")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    all_ids = read_fasta_ids(args.queries)
    if not all_ids:
        sys.exit("ERROR: no sequences found in --queries")

    try:
        df = pd.read_csv(args.hits, sep="\t", names=COLS, header=None, quoting=3)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=COLS)
    n_raw = len(df)

    # --- step 1: sequence-based filters ---
    keep = (df.pident >= args.min_identity) & (df.qcovs >= args.min_qcov)
    if args.max_evalue is not None:
        keep &= df.evalue <= args.max_evalue
    df = df[keep].copy()
    n_seq = len(df)
    with_seq = set(df.qaccver)

    pdb = df.saccver.map(split_pdb)
    df["pdb_id"] = [p[0] for p in pdb]
    df["chain"] = [p[1] for p in pdb]

    # --- step 2: experimental metadata ---
    ids = sorted(set(df.pdb_id.dropna()))
    if ids:
        meta = load_metadata(ids, args.pdb_meta, offline=args.offline)
    else:
        meta = pd.DataFrame(columns=META_COLS)
    df = df.merge(meta, on="pdb_id", how="left")
    df["method"] = df["method"].fillna("NA")
    df = df.sort_values(["qaccver", "bitscore", "pident"], ascending=[True, False, False])
    df.to_csv(out / "templates_annotated.tsv", sep="\t", index=False)

    # --- step 3: experimental filters ---
    n_unknown = int((df.method == "NOT_FOUND").sum())
    df = df[df.method != "NOT_FOUND"]
    if args.method is not None:
        df = df[df.method.str.lower().isin(args.method)]
    n_method = len(df)
    if args.max_resolution is not None:
        has_res = df.resolution.notna()
        too_low = has_res & (df.resolution > args.max_resolution)
        no_res = ~has_res
        df = df[~too_low & ~(no_res & args.drop_no_resolution)]
    df = df.sort_values(["qaccver", "bitscore", "pident", "resolution"],
                        ascending=[True, False, False, True])
    best = df.groupby("qaccver", sort=False).head(args.top_n)

    with_t = set(df.qaccver)
    no_t = [i for i in all_ids if i not in with_t]

    df.to_csv(out / "templates_filtered.tsv", sep="\t", index=False)
    best.to_csv(out / "templates_best.tsv", sep="\t", index=False)
    (out / "peptides_with_template.txt").write_text("\n".join(i for i in all_ids if i in with_t) + "\n")
    (out / "peptides_no_template.txt").write_text("\n".join(no_t) + "\n")

    meth_txt = "any" if args.method is None else ", ".join(sorted(args.method))
    if args.max_resolution is None:
        res_txt = "none"
    else:
        res_txt = (f"<= {args.max_resolution} A where a resolution exists; entries without one "
                   + ("dropped" if args.drop_no_resolution else "kept"))
    summary = (
        f"Query peptides:                         {len(all_ids)}\n"
        f"Raw hits:                               {n_raw}\n"
        f"Sequence thresholds: identity >= {args.min_identity}%, query coverage >= {args.min_qcov}%"
        + (f", E-value <= {args.max_evalue}" if args.max_evalue is not None else "") + "\n"
        f"  Hits passing:                         {n_seq}\n"
        f"  Peptides with >=1 hit:                {len(with_seq)} ({len(with_seq) / len(all_ids):.1%})\n"
        f"Experimental filter: method = {meth_txt}; resolution {res_txt}\n"
        f"  Hits dropped (PDB ID unknown to RCSB):{n_unknown:>6}\n"
        f"  Hits after method filter:             {n_method}\n"
        f"  Hits passing all filters:             {len(df)}\n"
        f"  Peptides with >=1 template:           {len(with_t)} ({len(with_t) / len(all_ids):.1%})\n"
        f"  Peptides with no template:            {len(no_t)}\n"
    )
    (out / "summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()