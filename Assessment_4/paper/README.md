## Data

The genotype and phenotype data used in this study are **not included**
in this repository due to copyright restrictions. They are publicly
available at the CSIRO Data Access Portal:

> Li, Z., Zhu, Q.-H., Moncuquet, P., & Wilson, I. (2022).
> *Cotton Variety Phenotype and Genotype Data.*
> https://doi.org/10.25919/k18n-nk98

Download the following three files and place them in the `paper/`
folder before rendering:

| File | Description |
|------|-------------|
| `Genotype_data.txt` | SNP marker scores for 256 lines × 6,523 markers |
| `Genetic_map.txt` | Chromosome assignment and physical position (bp) |
| `Phenotype_data.txt` | Four fibre quality traits across three field seasons |

---

## How to Reproduce

### Requirements

- R ≥ 4.1
- Quarto ≥ 1.3

### R Packages

Install all required packages by running in R:

```r
install.packages(c(
  "tidyverse", "glmnet", "irlba", "lme4breeding",
  "gtsummary", "kableExtra", "GGally", "patchwork",
  "gt", "broom", "scales", "data.table"
))
```

### Render the Paper

1. Download the data files from the link above and place them in
   `Assessment4/paper/`
2. Open `QTL_Analysis_restructured_VFinal.qmd` in RStudio
3. Run the following in the terminal:
