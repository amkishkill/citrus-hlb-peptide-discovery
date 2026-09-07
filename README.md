# Citrus HLB antiproteolysis peptide discovery

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
| 1 | Collect and annotate human gut microbiome sORFs | [DBsmORF](http://104.154.134.205:3838/DBsmORF/) | done (manual download) |
| 2 | Remove duplicates / non-standard AA sequences | `scripts/01_filter_dbsmorf.py` | done |
| 3 | Classify antimicrobial peptides | [AMPDiscover](https://biocom-ampdiscover.cicese.mx/) | TODO |
| 4 | Deduplicate at 50% identity / 80% coverage | [MMseqs2](https://github.com/soedinglab/mmseqs2) | TODO |
| 5 | Homology search | BLAST | TODO |
| 6 | Structure prediction (peptides + PUB21) | Modeller / AlphaFold | TODO |
| 7 | Peptide-protein docking vs PUB21 | AutoDock Vina (paper) / peptide-aware alternatives | TODO |

## Step 1-2: DBsmORF download

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

This produced 416,158 raw entries (before dedup / non-standard-AA
removal). The paper reports 122,826 sORFs remaining after that cleanup.

Place the raw TSV at `data/dbsmorf_hmp_gut_raw.tsv` and run:

```bash
pip install -r requirements.txt
python scripts/01_filter_dbsmorf.py \
    --input data/dbsmorf_hmp_gut_raw.tsv \
    --out-prefix data/dbsmorf_hmp_gut_clean
```

This writes `data/dbsmorf_hmp_gut_clean.fasta` and `.tsv`.

## Notes / gaps vs. the original paper

The paper's Methods and SI do not disclose several implementation
details needed for exact reproduction:

- The Modeller template(s) used to build the PUB21 homology model
  (no PDB entry exists for PUB21).
- AutoDock Vina docking box coordinates, exhaustiveness, and other
  run parameters.
- The exact BLAST reference database/parameters used in the homology
  search step.

These will be documented as design decisions in later scripts as the
pipeline is built out.
