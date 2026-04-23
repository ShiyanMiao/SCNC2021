# Simulation + Analysis

QTL mapping of cotton fibre quality traits using a MAGIC population.  
Compares three marker selection methods: single-marker scan, PC-adjusted scan, and lasso regression.

---

## Files

| File | Description |
|---|---|
| `QTL_Analysis_Final.qmd` | Main Quarto report |
| `references.bib` | BibTeX references |
| `Genotype_data.txt` | Genotype matrix (256 lines × 6,523 markers) |
| `Genetic_map.txt` | Genetic map (SNP ID, chromosome, position in cM) |
| `Phenotype_data.txt` | Phenotype records stacked across three seasons |

---

## Requirements

- R ≥ 4.3.0
- Quarto ≥ 1.4
- TinyTeX (for PDF output): `tinytex::install_tinytex()`

R packages:

```r
install.packages(c(
  "tidyverse", "data.table", "glmnet", "knitr", "kableExtra",
  "scales", "irlba", "gt", "gtsummary", "GGally", "purrr", "patchwork"
))
```

---

## How to Reproduce

1. Place the three data files in this folder.
2. Run the following in a terminal:

```bash
quarto render QTL_Analysis_Final_Sent.qmd
```

Or open the file in RStudio and click **Render**.

The power analysis simulation takes 30–60 minutes on first render. Results are cached to `Result/power_results.rds` and reused on subsequent renders.

---

## Author

Shiyan Miao — SCNC2021, Australian National University, 2026
