#  Splicing Pathology (SpliPath)

## Introduction

SpliPath is designed to identify and functionally cluster splice-altering variants in WGS that have similar effects on RNA splicing. We refer to these variants as collapsed splicing quantitative trait loci (csQTLs). SpliPath aims to address two main difficulties in explaining missing heritability in rare disorders: effective method to functionally interpret genetic variants and increase statistical power to establish associations between rare variants and phenotypes. First, SpliPath links the prediction of SpliceAI with reference transcriptomics data to identify genetic variants that induce splice changes actually occuring in disease-relevant transcriptomics profiles. Second, SpliPath aggregates variants with similar functional consequences into csQTLs for more powerful genetic association analyses.

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
