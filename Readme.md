# Structural Breaks Detection: 'Simulation code'

---

The R code for illustrations, simulations and real data analyses in the paper "Multiple Testing of Local Extrema for Detection of Multiple Change Points in Linear Models" by Zhibing He, Dan Cheng, and Yunpeng Zhao (2021) available via (Arxiv link will be available soon)

## Dependencies

---

The following R packages are required for plots and tables presented in the paper.

- 'dSTEM' implementing our multiple testing for change points detection method (installed via  [r code](devtools::install_github('zhibinghe/ChangePoint') ))
- 'xtable' for exporting tables in latex format
- 'latex2exp' for mathematical symbols in plots
- 'foreach' for parallel computing

The following packages are additionally required for comparations with other methods

- 'not' for the Narrowest-Over-Threshold method
- 'strucchange' for B&P method
- 'nsp' for Narrowest-Significance-Pursuit method (installed from [github](https://github.com/pfryz/nsp) )

## Main Files

1. Illustration
2. Simulation 
3. Comparation with other methods
4. Real Data Analysis

