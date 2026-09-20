# Data acquisition and provenance

The raw input is **not stored in this GitHub repository**.

## Source

- Resource: EMBL-EBI Antimicrobial Resistance (AMR) Portal
- Release used: **2026-07**
- Input: combined phenotype-genotype data
- Expected filename: `phenotype_genotype_merged.csv.gz`
- Source portal: https://www.ebi.ac.uk/amr/
- Release directory used during the project: https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2026-07/

## Integrity check

The exact compressed file used for the reported analysis had SHA-256:

```text
eb4111f5d3c58fcb7e0b29078182929847614c847cc0aa19c23cb4a3e5dbc92b
```

After downloading the matching release file, place it at:

```text
data/raw/phenotype_genotype_merged.csv.gz
```

Create `data/raw/` if it does not exist. The `.gitignore` file prevents this raw file from being committed accidentally.

### Check SHA-256 on Windows PowerShell

```powershell
Get-FileHash .\data\raw\phenotype_genotype_merged.csv.gz -Algorithm SHA256
```

### Check SHA-256 on Linux/macOS

```bash
sha256sum data/raw/phenotype_genotype_merged.csv.gz
```

Only proceed if the checksum matches the value above when exact reproduction of the reported run is required.

## Attribution and rights

The source dataset is third-party data from EMBL-EBI. This project does not claim ownership of it. The AMR Portal licensing information should be checked at the time of reuse and the source/release should be cited as required by the portal. The analysis filters and transforms the source records to derive isolate-level features, predictions, metrics, and figures.

## Generated data

Running `src/run_pipeline.py` creates `data/processed/feature_matrix.csv`, `data/processed/labels.csv`, and `data/processed/feature_schema.json` locally. They are excluded from this GitHub package because they are directly derived from the third-party input; the reported aggregate/per-isolate research outputs needed to audit the published run are retained under `results/`.
