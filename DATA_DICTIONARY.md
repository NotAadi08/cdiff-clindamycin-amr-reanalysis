# Data dictionary

- `data/raw/phenotype_genotype_merged.csv.gz`: source merged phenotype-genotype file supplied for the project.
- `data/processed/feature_matrix.csv`: isolate-by-feature binary matrix produced by the pipeline.
- `data/processed/labels.csv`: one phenotype label per retained BioSample.
- `data/processed/feature_schema.json`: ordered model feature names.
- `results/cv_metrics.csv`: repeated cross-validation metrics for all candidate models.
- `results/oof_predictions.csv`: one out-of-fold probability and prediction per isolate.
- `results/confusion_matrix.csv`: rows=true class, columns=predicted class; class order 0=susceptible, 1=resistant.
- `results/feature_importance.csv`: absolute coefficient magnitude for logistic regression or built-in tree importance for tree models.
- `results/top_feature_prevalence.csv`: feature prevalence by phenotype group.
- `results/summary.json`: compact machine-readable summary and source-file SHA-256 checksum.
