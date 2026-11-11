# FarmTest report

To reproduce [FinalReport.pdf](FinalReport.pdf), open [FinalReport.Rmd](FinalReport.Rmd) in RStudio and click **Knit to PDF**. This uses FarmTest 2.2.0, ggplot2, knitr, rmarkdown and LaTeX.

- `00_requirements.R` loads the packages.
- `01_GenerateData.R` simulates the gene readings with t₃ noise.
- `02_TTests.R` runs the t-tests and adjusts their p-values.
- `03_Evaluate.R` scores both methods and calculates the gene 1 and PCA examples.
- `04_Plots.R` draws the report figures.
- `run_analysis.R` runs the simulation and fits both methods each time you knit the report.
- `data/GSE7390_series_matrix.txt.gz` is the original GEO expression matrix for the planned breast tumor comparison.
- `results/illustration_scores.csv` contains the results for FarmTest and the t-tests shown in the report.
