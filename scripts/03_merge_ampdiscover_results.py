"""
Step 3b: Merge the per-batch AMPDiscover result CSVs and filter down to
predicted antimicrobial peptides (AMPs).

After submitting each data/ampdiscover_batches/batch_NNN.fasta file to
AMPDiscover (https://biocom-ampdiscover.cicese.mx) using the general
ProtDCal-AMP_RF model, download each result CSV into a single folder
(e.g. data/ampdiscover_results/) and run this script.

Because AMPDiscover's exact label strings for "active"/"non-active" are
not documented in its Help page, this script first PRINTS the unique
values found in the result column so you can confirm the correct
--active-label before trusting the filtered output.

Usage:
    # first pass -- just inspect the column values
    python scripts/03_merge_ampdiscover_results.py \
        --results-dir data/ampdiscover_results \
        --out data/ampdiscover_merged.csv

    # second pass, once you know the correct label string
    python scripts/03_merge_ampdiscover_results.py \
        --results-dir data/ampdiscover_results \
        --out data/ampdiscover_merged.csv \
        --active-label Active \
        --require-ad-reliable
"""

import argparse
from pathlib import Path

import pandas as pd

RESULT_COL = "qsar_Result_ProtDCal-AMP_RF"
SCORE_COL = "qsar_Score_ProtDCal-AMP_RF"
AD_COL = "AD-R_Consensus_ProtDCal-AMP_RF"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Directory containing downloaded AMPDiscover CSV files (one per batch)",
    )
    parser.add_argument("--out", required=True, help="Path to write merged CSV")
    parser.add_argument(
        "--active-label",
        default=None,
        help="Exact string AMPDiscover uses to mark a peptide as an active AMP "
        "(e.g. 'Active'). If omitted, no filtering is applied -- the script "
        "just merges and reports the unique values found so you can pick one.",
    )
    parser.add_argument(
        "--require-ad-reliable",
        action="store_true",
        help="Additionally require AD-R_Consensus == 'IN' (prediction inside "
        "the applicability domain) before counting a peptide as active",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    csv_files = sorted(results_dir.glob("*.csv"))
    if not csv_files:
        raise SystemExit(f"No CSV files found in {results_dir}")

    print(f"Found {len(csv_files)} result files:")
    frames = []
    for f in csv_files:
        df = pd.read_csv(f)
        df["_source_batch"] = f.name
        frames.append(df)
        print(f"  {f.name}: {len(df)} rows")

    merged = pd.concat(frames, ignore_index=True)
    print(f"\nTotal merged rows: {len(merged)}")

    if RESULT_COL not in merged.columns:
        print(
            f"\nWARNING: expected column '{RESULT_COL}' not found. "
            f"Available columns: {list(merged.columns)}"
        )
        merged.to_csv(args.out, index=False)
        print(f"Wrote unfiltered merge to {args.out} -- inspect columns and re-run.")
        return

    print(f"\nUnique values in '{RESULT_COL}': "
          f"{merged[RESULT_COL].unique().tolist()}")
    if AD_COL in merged.columns:
        print(f"Unique values in '{AD_COL}': "
              f"{merged[AD_COL].unique().tolist()}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.active_label is None:
        merged.to_csv(out_path, index=False)
        print(
            f"\nWrote unfiltered merged CSV to {out_path}\n"
            f"Re-run with --active-label '<value from list above>' to filter."
        )
        return

    mask = merged[RESULT_COL].astype(str).str.strip().str.lower() == (
        args.active_label.strip().lower()
    )
    n_active_raw = mask.sum()
    print(f"\nPeptides matching --active-label '{args.active_label}': {n_active_raw}")

    if args.require_ad_reliable and AD_COL in merged.columns:
        ad_mask = merged[AD_COL].astype(str).str.strip().str.upper() == "IN"
        combined_mask = mask & ad_mask
        print(
            f"Of those, also AD-reliable ('IN'): {combined_mask.sum()}"
        )
        final_mask = combined_mask
    else:
        final_mask = mask

    filtered = merged[final_mask].copy()
    filtered.to_csv(out_path, index=False)
    print(f"\nWrote {len(filtered)} predicted-AMP rows to {out_path}")
    print(
        "\nFor reference, the paper reports 40,837 AMPs identified from "
        "122,826 sORFs at this stage."
    )


if __name__ == "__main__":
    main()