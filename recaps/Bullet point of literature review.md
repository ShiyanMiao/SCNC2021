# Literature Review Coverage Points

## Background and Significance of Marker-Assisted Selection (MAS)
- **Role and advantages** of MAS in modern plant and animal breeding:
  - Increasing efficiency
  - Reducing costs
- **Central role** of QTL analysis in the MAS workflow

## Basic Workflow of QTL Analysis
- **Experimental design**: population type, population size, marker density, etc.
- **Sample collection and genotyping technologies**
- **Phenotypic data collection and quality control**
- **Choice of statistical models** and association analysis methods

## Key Factors Affecting QTL Detection Power
- Relationship between sample size and statistical power
- Marker effect size
- Influence of heritability
- Population structure and genetic diversity
- Impact of statistical model choice on detection outcomes

## Commonly Used Statistical Models in QTL Analysis
- Single marker models (SMM)
- Interval Mapping (IM)
- Composite Interval Mapping (CIM)
- Multiple marker models and mixed linear models (MLM, GWAS framework)
- Bayesian methods and whole-genome prediction models

## Progress and Challenges in Optimizing Experimental Design
- Consequences of over-sampling and under-sampling
- Effect of genetic homogeneity in populations on detection power
- Balancing marker density and genome coverage
- Resource allocation and cost–benefit considerations

## Applications of Simulation in QTL Research
- Simulation-based approaches for evaluating experimental design
- Advantages and limitations of existing QTL power analysis and simulation tools:
  - QTL Cartographer
  - R/qtl
  - GAPIT
- How simulation frameworks can help researchers understand the impact of design parameters

## Research Gaps and Needs
- Lack of user-friendly interactive tools for experimental design targeted at breeders without statistical backgrounds
- Lack of customizable QTL simulation and power evaluation frameworks for specific crops or traits

## Positioning and Contribution of This Project
- Develop a simulation-based framework for evaluating QTL experimental designs
- Build an interactive R Shiny application
- Provide practical recommendations for agronomists and statisticians regarding:
  - Sample size
  - Heritability
  - Population structure
- Potential application to ongoing QTL studies on blackleg severity in canola

---

## Suggested References and Links
1. **Broman, K. W., & Sen, Ś. (2009)**. *A Guide to QTL Mapping with R/qtl*. Springer.  
   [https://rqtl.org/book/](https://rqtl.org/book/)

2. **Lander, E. S., & Botstein, D. (1989)**. *Mapping Mendelian factors underlying quantitative traits using RFLP linkage maps*. Genetics, 121(1), 185–199.  
   [https://pubmed.ncbi.nlm.nih.gov/2563713/](https://pubmed.ncbi.nlm.nih.gov/2563713/)

3. **Churchill, G. A., & Doerge, R. W. (1994)**. *Empirical threshold values for quantitative trait mapping*. Genetics, 138(3), 963–971.  
   [https://www.genetics.org/content/138/3/963](https://www.genetics.org/content/138/3/963)

4. **Zeng, Z.-B. (1994)**. *Precision mapping of quantitative trait loci*. Genetics, 136(4), 1457–1468.  
   [https://www.genetics.org/content/136/4/1457](https://www.genetics.org/content/136/4/1457)

5. **Yu, J., Pressoir, G., Briggs, W. H., et al. (2006)**. *A unified mixed-model method for association mapping*. Nature Genetics, 38(2), 203–208.  
   [https://doi.org/10.1038/ng1702](https://doi.org/10.1038/ng1702)

6. **Xu, S. (2003)**. *Theoretical basis of the Beavis effect*. Genetics, 165(4), 2259–2268.  
   [https://www.genetics.org/content/165/4/2259](https://www.genetics.org/content/165/4/2259)
