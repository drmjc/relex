# Normalization and Statistical Analysis of Quantitative Proteomics Data Generated by Metabolic Labeling

This is an R Package, developed 2007-9 to support statistical analysis of quantitative proteomics data.

At the time, the field of microarray analysis had developed a suite of tools for data normalisation and empirical Bayes analysis of differential expression. This was the first study to apply these methods to quantgiative proteomics data, which had many parallels: radioactive-label swaps were equivalent to dye-swaps; there were biases linked to abundance; there were many peptides from a whole protein.

This package essentially helps to bridge the two fields, and allow *limma*, a popular R/Bioconductor package for microarray analysis to be compatible with proteomics data. We developed this to handle 'RelEx' data [PMID 14670053], which no longer appears to be active (http://fields.scripps.edu/relex/).

This code has not been developed since 2009, and is no longer actively supported. As of March 2019, `R CMD BUILD` does seem to still work, so the package at least appears valid. YMMV.

# Installation
## Dependencies
1. Install R
2. install limma: `R -e 'install.packages("limma")'`

## Build and install the 'relex' R package
1. obtain a copy of this repo, and unzip it
2. `R CMD BUILD relex`
3. `R CMD install relex_1.0.1.tar.gz`

# Usage
* start `R`, and load the relex package: `library(relex)`
* check the help for these key functions for more info:

```
?relex
?import.relex.experiment
?normalize.relex.experiment
```

# Citation
* Normalization and Statistical Analysis of Quantitative Proteomics Data Generated by Metabolic Labeling. Lily Ting, Mark J. Cowley, Seah Lay Hoon, Michael Guilhaus, Mark J. Raftery, Ricardo Cavicchioli. *Molecular & Cellular Proteomics*, October 1, 2009, First published on July 14, 2009, 8 (10) 2227-2242; DOI: 10.1074/mcp.M800462-MCP200
* https://www.ncbi.nlm.nih.gov/pubmed/19605365,

