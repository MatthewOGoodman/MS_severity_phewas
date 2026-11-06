# MS_severity_phewas

Analysis code for:

> Zhang T\*, Goodman M\*, Zhu F, Healy B, Carruthers R, Chitnis T, et al. (\*equal contribution)
> **Phenome-wide examination of comorbidity burden and multiple sclerosis disease severity.**
> *Neurology: Neuroimmunology & Neuroinflammation* 2020;7(6):e864. doi:[10.1212/NXI.0000000000000864](https://doi.org/10.1212/NXI.0000000000000864)

A phenome-wide association study (PheWAS) of electronic-health-record comorbidities against multiple sclerosis severity (MSSS), in a discovery cohort, with replication and meta-analysis in a second cohort.

## Files

| File | What it does |
|---|---|
| `run_PheWAS` | Runs the PheWAS of EHR phecodes on MSSS with the `PheWAS` package, adjusting for sex, race, smoking, age at onset, onset type, disease duration and healthcare utilization. |
| `replication_analysis.R` | Aligns discovery and replication results by phecode and applies the case-count thresholds (primary and secondary analyses). |
| `meta_analysis.R` | Fixed-effect meta-analysis of the two cohorts with `metafor`. |
| `plot_OR_vs_prev.R` | Odds ratio against comorbidity prevalence. |
| `plots_correlation_manhattan.R` | Clustering of p-values, correlation plots and Manhattan plots. |

These scripts are the analysis record. They need the cohorts' individual-level EHR data, which are not included and cannot be shared, and they read and write files in the working directory.

## Licence

MIT ([LICENSE](LICENSE)).
