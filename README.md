# PsychAD gene-set AUROC tester

Run from this folder:

```sh
pip install -r requirements.txt
streamlit run app.py --server.port 8503
```

Uses the supplied `data/Donghoon_Lee...Panos_Roussos/PsychAD_SupplementaryTable6.csv` from Lee et al., *Single-cell atlas of transcriptomic vulnerability across brain disorders*, https://doi.org/10.1038/s41586-025-09573-z. The table has 1,016,196 rows, 27 cell types and four phenotypes (dx_AD, CERAD, Braak, Dementia). The paper describes donor-level pseudobulk Dreamlet differential expression with mixed models and age, sex, PMI and expression-quality covariates. Figure 6d reports meta-analysis across brain banks.

Paste human gene symbols. Results update automatically on input changes. Leave the optional background box blank to use all source genes. Background is the eligible universe including targets; defaults to all source genes. Matching is case-sensitive and deduplicated. Each cell-type/phenotype comparison uses only available finite scores, excludes targets from the non-target comparison, and reports both group sizes. No gene-level significance filtering is applied.

Default ranking: signed t-statistic. Optional rankings: absolute t-statistic and signed effect estimate. AUROC = Mann–Whitney U / (n_target × n_comparison), with half credit for ties. Signed AUROC >0.5 means more positive associations, <0.5 more negative associations. Two-sided asymptotic Mann–Whitney tests include tie and continuity corrections. BH correction covers all valid tests, across all 27 cell types and four phenotypes. Stars use q<0.05/0.01/0.001. Empty groups remain untestable. Gene correlation is not modelled; p-values are exploratory competitive gene-set tests, not patient classification accuracy or a polygenic risk score.

The app defaults to a light theme. The heatmap shows all cell types and appears before the results table (sorted by descending sign(AUROC − 0.5) × −log10(p-value)), with downloadable PNG and full-precision CSV (including per-comparison target membership). The app's About panel describes phenotype coding and methods. Source `p.value` and `FDR` are not reused as gene-set significance.

Tests: `python -m unittest discover -s tests -v`.

## Default gene set

The default is the user-supplied list of 15 positively tau-correlated genes shown in Figure 2b of Huang W, Bartosch AM, Xiao H, et al. (2021), *An immune response characterizes early Alzheimer’s disease pathology and subjective cognitive impairment in hydrocephalus biopsies*, Nature Communications 12:5659, https://doi.org/10.1038/s41467-021-25902-y. This exact list matches that figure; the initially supplied Huang Y-N et al. (2025) plasma-proteomics citation (PMID 40810263) is a different study.

C4B, PIRT and SERPINA3 remain current human symbols. All three and RP11-203I2 are absent from the supplied PsychAD Table 6. RP11-203I2 is a clone identifier with no verified unique replacement here. The default retains all 15 submitted labels and reports unmatched genes rather than substituting other genes.

## Aging dataset

Select **Aging** to load `data/Hui_Yang...Panos_Roussos/Supplementary_Data_4.rds.filtered.csv.gz` from Hui Yang et al., *Lifespan single-cell transcriptomic atlas of the human prefrontal cortex*, https://www.nature.com/articles/s41586-026-10271-7. The supplied file contains 26 cell types × four age periods (104 tests). `group`, `t`, and `logFC` map to the analysis phenotype, statistic, and estimate columns. Positive signed scores indicate increasing expression with age within the indicated period, not a contrast between age groups. The eligible universe is restricted to the supplied filtered file; no additional significance filtering is applied. BH correction covers valid tests in the selected dataset only. The existing AD dataset remains the default.

## Whole-lifespan analysis

`Rscript scripts/prepare_lifespan.R` filters Supplementary Data 7 by excluding IDs matching `^ENSG`, keeps ID, assay, F, coef_1, coef_2, and writes `Supplementary_Data_7.filtered.csv`. No significance filter is applied. Select **Aging — whole lifespan** for the full-age polynomial model; **Aging — within age bins** uses Supplementary Data 4. The whole-lifespan default ranks overall F-statistics (unsigned association strength); optional signed coefficients represent the linear/quadratic orthogonal polynomial terms of log2(age+1), not a single monotonic aging slope. FDR is across the 26 valid cell-type tests.

Whole-lifespan mode now shows F-statistic and coarse AveExpr slope side by side. The R preparation joins the union of S7 genes and non-ENSG S4 genes, retaining unavailable measures as NA. Slopes require four finite bin means and use OLS against x=1,2,3,4 (equal bin spacing), not donor ages. Each column uses its own available gene universe; BH correction spans both columns. The slope is descriptive and unadjusted, not a reconstructed donor-level age effect.

Whole-lifespan heatmap columns now begin with four AveExpr z-score AUROCs (development, young, middle, late), followed by slope and F-statistic last. Z-scores standardize each gene-cell pair across its four bin means using sample SD (ddof=1). Missing or zero-variance profiles are excluded from z-score tests; no imputation. BH correction covers all valid tests across six columns.
