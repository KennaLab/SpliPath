#  Splicing Pathology (SpliPath)

## Introduction

SpliPath is developed to nominate ultra-rare splicing quantitative trait loci (ur-sQTL) in genomewide analyses of paried DNAseq and RNA splicing datasets. It serves to mitigate key challenges of distinguishing rare bona fide aberrant splicing events from abundant artefects, and to identify causal genetic factors responsible for these events. SpliPath provides genomewide splice junction annotation, ur-sQTL nomination and an interactive Shiny-based data visualization and analysis tool. It enables users to look beyond consensus splice sites for causal genetic variants in deeper intronic regions, and can be applied to analyse any dataset that includes matched DNA and RNA sequencing data for donors believed to harbour rare pathogenic splice-altering variants. 

## Installation

### Clone repo
Clone repo
```{sh}
git clone https://github.com/KennaLab/SpliPath.git
```
Open an R session and install using devtools 
```{r}
devtools::install("SpliPath") # if the SpliPath directory is not in your current directory,  provide the full path to the cloned repo
```

### Install dependencies

* [RVAT](https://github.com/KennaLab/rvat)
```{r}
remotes::install_github("kennalab/rvat")
```

## Tutorials

* [SpliPath basics](https://github.com/KennaLab/SpliPath/blob/main/SpliPath_Tutorial.md)
