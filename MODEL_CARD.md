# Model card

**Task:** proof-of-concept prediction of clindamycin resistance in *Clostridioides difficile* from AMRFinderPlus resistance determinants.

**Input:** binary presence/absence of retained AMR determinants. The packaged run retained three features after requiring presence in at least two isolates.

**Training cohort:** 102 isolates (35 resistant, 67 susceptible) from the EMBL-EBI AMR Portal merged phenotype-genotype release used in the study.

**Model selection:** logistic regression, random forest, and Extra Trees were compared by repeated stratified 5-fold cross-validation (5 repeats). The model with highest mean AUROC was selected. A separate shuffled stratified 5-fold out-of-fold run generated the reported per-isolate predictions and plots.

**Intended use:** research reproducibility and methodological demonstration only.

**Not intended for:** clinical diagnosis, treatment selection, or patient care. The cohort is small and database-source imbalance may inflate apparent performance. External validation is required.
