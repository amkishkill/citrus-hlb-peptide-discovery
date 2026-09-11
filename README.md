# Citrus HLB antiproteolysis peptide discovery (reproduction)

Reproducing the computational pipeline from:

> Zhao P. et al. "Targeted MYC2 stabilization confers citrus Huanglongbing
> resistance." *Science* 388, 191-198 (2025). DOI: 10.1126/science.adq7203

Specifically this repo reproduces the pipeline in Fig. 5A / Methods that
identifies antiproteolysis peptides (APPs) which bind and inhibit the
citrus E3 ubiquitin ligase **PUB21**, indirectly stabilizing the MYC2
transcription factor and conferring HLB resistance.

## Pipeline stages

| # | Stage | Tool | Status |
|---|-------|------|--------|
| 1 | Collect + annotate human gut microbiome sORFs | [DBsmORF](http://104.154.134.205:3838/DBsmORF/) ([SmORFinder](https://github.com/bhattlab/SmORFinder) backend, default significance filters) | done (manual download) |
| 2 | Remove duplicates / non-standard AA sequences | `scripts/01_filter_dbsmorf.py` | done — 122,826 sequences, exact match to paper |
| 3 | Classify antimicrobial peptides | [AMPDiscover](https://biocom-ampdiscover.cicese.mx) (web only, `ProtDCal-AMP_RF` general model) | done — 33,491 AMPs (paper: 40,837; see note below) |
| 4 | Deduplicate at 50% identity / 80% coverage | MMseqs2 (`scripts/05_run_mmseqs_cluster.sh`) | done -  12,951 AMPs (paper: 14,950)|
| 5 | Homology search | BLAST | TODO |
| 6 | Structure prediction (peptides + PUB21) | Modeller | TODO |
| 7 | Peptide-protein docking vs PUB21 | AutoDock Vina (paper) / peptide-aware alternatives | TODO |


## Step 1-2: DBsmORF download + cleaning

Downloaded from the DBsmORF **Download** tab with:

- Source Database: `HMP`
- Filter by Body Site: `Gut`
- Filter by taxa / NCBI UID / SRA ID / Gender: `Keep All`
- Individual model significance filters: SmORFinder defaults
  (pHMM E-value ≤ 1e-6, DSN1/DSN2 P(smORF) ≥ 0.9999)
- Overlapping model significance filters: SmORFinder defaults
  (pHMM E-value ≤ 1, DSN1/DSN2 P(smORF) ≥ 0.5)
- Download Fields: defaults + Protein Sequence
- File format: TSV

This produced 416,158 raw entries. After deduplication and removal of
non-standard-amino-acid sequences (`scripts/01_filter_dbsmorf.py`), **122,826
sequences remained — an exact match to the paper.**

```bash
pip install -r requirements.txt
python scripts/01_filter_dbsmorf.py \
    --input data/dbsmorf_hmp_gut_raw.tsv \
    --out-prefix data/dbsmorf_hmp_gut_clean
```

## Step 3: AMPDiscover classification

AMPDiscover (https://biocom-ampdiscover.cicese.mx) is a **web-only** tool —
no downloadable package or API exists. It also enforces a 1 MB per-submission
file size limit, so the cleaned FASTA (~7.1 MB) had to be split into batches,
each submitted manually through the web form, then merged.

Model used: **AMPDiscover (ProtDCal)**, the general `ProtDCal-AMP_RF` model.

```bash
# split into <1MB chunks
python scripts/02_split_fasta_for_ampdiscover.py \
    --input data/dbsmorf_hmp_gut_clean.fasta \
    --out-dir data/ampdiscover_batches

# --- manually submit each batch_NNN.fasta to AMPDiscover, download results
#     into data/ampdiscover_results/ ---

# merge + filter to predicted AMPs
python scripts/03_merge_ampdiscover_results.py \
    --results-dir data/ampdiscover_results \
    --out data/ampdiscover_active_amps.csv \
    --active-label AMP
```

### Discrepancy vs. the paper

This step yielded **33,491 predicted AMPs**, versus the paper's reported
**40,837** (~82% of the target). Before accepting this, we ruled out:

- **Data loss during batching/submission** — confirmed the sum of rows across
  all downloaded batch result CSVs equals 122,826, the full input count.
  Nothing was dropped.

Two plausible, non-diagnosable explanations remain:

1. The paper may have pooled results across multiple AMPDiscover models
   (e.g. combining the general model with activity-specific ones) rather
   than using the general `ProtDCal-AMP_RF` model alone — the Methods text
   doesn't specify precisely enough to rule this out.
2. AMPDiscover is a live web service with no version pinning; the deployed
   model/preprocessing may have changed since the paper's analysis was run.

Given both explanations point to information the paper doesn't disclose and
which can't be independently verified, **33,491 AMPs is treated as an
acceptable reproduction** and used as the working set for all downstream
steps.

## Step 3c: Extract AMP sequences

AMPDiscover's output CSV contains only IDs and prediction columns, not the
actual sequences. Join back against the cleaned FASTA:

```bash
python scripts/04_extract_amp_fasta.py \
    --clean-fasta data/dbsmorf_hmp_gut_clean.fasta \
    --amp-csv data/ampdiscover_active_amps.csv \
    --id-col IDs \
    --out data/dbsmorf_amps.fasta
```

## Step 4: MMseqs2 clustering

Per the Methods: *"Representative AMPs were chosen by using MMseqs2, with an
amino acid identity of 50% and coverage of 80%."*

MMseqs2 is a compiled binary, separate installation is required :

```bash

# Linux (conda)
conda install -c conda-forge -c bioconda mmseqs2


Then run:

```bash
bash scripts/05_run_mmseqs_cluster.sh \
    data/dbsmorf_amps.fasta \
    data/mmseqs_out
```

This produces `data/mmseqs_out_rep_seq.fasta` (one representative sequence
per cluster) and `data/mmseqs_out_cluster.tsv` (cluster membership). For
reference, the paper reports **14,950 representative AMPs** at this stage.

## Notes / gaps vs. the original paper

The paper's Methods and SI do not disclose several implementation
details needed for exact reproduction:

- The Modeller template(s) used to build the PUB21 homology model
  (no PDB entry exists for PUB21).
- AutoDock Vina docking box coordinates, exhaustiveness, and other
  run parameters.
- The exact BLAST reference database/parameters used in the homology
  search step.
- Whether AMPDiscover results were pooled across multiple models (see
  discrepancy note above).
