# C. difficile clindamycin AMR reanalysis

[![DOI](https://zenodo.org/badge/1378243404.svg)](https://doi.org/10.5281/zenodo.22857761)

**Archived release v1.0.0:** https://doi.org/10.5281/zenodo.22857762

[![ORCID](https://img.shields.io/badge/ORCID-0009--0000--1843--3507-A6CE39)](https://orcid.org/0009-0000-1843-3507)

A reproducible bioinformatics reanalysis of genomic antimicrobial-resistance determinants associated with clindamycin resistance in *Clostridioides difficile*, using the EMBL-EBI Antimicrobial Resistance Portal 2026-07 phenotype-genotype release.

**Author:** Aditya Agrawal  
**Affiliation:** School of Health Sciences and Technology, UPES, Dehradun, Uttarakhand, India  
**Program:** B.Tech. Bioinformatics  
**ORCID:** [0009-0000-1843-3507](https://orcid.org/0009-0000-1843-3507)

## Research question

Does an interpretable machine-learning classifier add useful binary classification information beyond established AMR-marker rules for clindamycin resistance in the analyzed *C. difficile* cohort?

## Main finding

The analyzed cohort contained **102 isolates** (35 resistant, 67 susceptible). Logistic-regression out-of-fold classification achieved **96.1% accuracy**, **94.3% balanced accuracy**, **88.6% sensitivity**, and **100% specificity**. However, an **`erm(B)`-only rule produced exactly the same 102 binary classifications**. In this cohort, `cfr(C)` occurred only in `erm(B)`-positive isolates and did not rescue any `erm(B)`-negative resistant isolate.

This repository therefore presents a **reproducibility and baseline-audit study**, not a claim of a novel clinical predictor. The packaged model is **not for clinical use**.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE.md
├── requirements.txt
├── environment.txt
├── MODEL_CARD.md
├── DATA_DICTIONARY.md
├── src/
│   ├── run_pipeline.py
│   ├── followup_analysis.py
│   └── predict.py
├── data/
│   └── README.md
├── models/
│   └── cdiff_clindamycin_model.joblib
├── results/
├── figures/
├── docs/
└── manuscript/
```

## Data

The raw third-party EMBL-EBI dataset is intentionally **not committed to this repository**. See [`data/README.md`](data/README.md) for the exact release, expected filename, checksum, attribution, and placement instructions.

Expected local path after download:

```text
data/raw/phenotype_genotype_merged.csv.gz
```

The exact input used for the reported run had SHA-256:

```text
eb4111f5d3c58fcb7e0b29078182929847614c847cc0aa19c23cb4a3e5dbc92b
```

## Reproduce the analysis

Use Python 3.13.x and create an isolated environment if possible.

```bash
python -m pip install -r requirements.txt
python src/run_pipeline.py
python src/followup_analysis.py
```

The first script regenerates the feature matrix, cross-validation metrics, out-of-fold predictions, model artifact, feature importance, and principal figures. The second regenerates the marker-only comparisons, source-stratified audit, false-negative audit, and baseline figure.

### Windows PowerShell

```powershell
py -m pip install -r requirements.txt
py src/run_pipeline.py
py src/followup_analysis.py
```

## Demonstrate the fitted model

The repository includes the fitted proof-of-concept model used for demonstration:

```bash
python src/predict.py --features "erm(B),cfr(C)"
```

This accepts presence/absence of the packaged AMR features. It is a research artifact and must not be used to guide patient treatment.

## Key outputs

- `results/cv_metrics.csv` - repeated 5-fold x 5-repeat cross-validation results.
- `results/oof_predictions.csv` - per-isolate out-of-fold predictions.
- `results/baseline_comparison.csv` - marker rules versus ML classification.
- `results/false_negative_audit.csv` - four resistant isolates missed by both ML and the `erm(B)` rule.
- `results/source_stratified_analysis.csv` - source-specific cohort composition and marker performance.
- `figures/roc_curve.png` and `figures/pr_curve.png` - out-of-fold discrimination plots.
- `figures/baseline_comparison.png` - direct marker-rule versus ML comparison.

## Important limitations

1. The cohort is small (102 isolates).
2. Source composition is strongly imbalanced: the NCBI-antibiogram subset contains only resistant isolates in this filtered cohort.
3. No independent external validation cohort is included.
4. The retained feature space is very small and includes an invariant feature (`cplR`).
5. High apparent accuracy must not be interpreted as clinical validity.
6. The central `erm(B)` association is established biology and is not claimed as a discovery of this project.

## Manuscript

The repository includes the corresponding manuscript under [`manuscript/`](manuscript/). Its framing deliberately distinguishes reproducibility/baseline auditing from methodological novelty.

## AI-assistance disclosure

OpenAI ChatGPT (GPT-5.6 Sol) assisted with literature searching, code development, statistical checking, documentation, and manuscript drafting/organization. The human author is responsible for reviewing the repository, rerunning the final analysis, verifying references and results, and taking responsibility for any public or journal submission. AI is not an author.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). A DOI field should be added after a Zenodo DOI is reserved/issued for the archival release.

## Licensing and third-party material

See [`LICENSE.md`](LICENSE.md). The EMBL-EBI source data are third-party material and are not covered by the repository's permissions for original project material.
