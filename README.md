<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# multinichenetr

<!-- badges: start -->

[![R build
status](https://github.com/browaeysrobin/multinichenetr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/browaeysrobin/multinichenetr/actions)
[![Coverage
Status](https://codecov.io/gh/browaeysrobin/multinichenetr/branch/master/graph/badge.svg?token=NKZBMJJDYA)](https://codecov.io/gh/browaeysrobin/multinichenetr)
<!-- badges: end -->

**multinichenetr: the R package containing multiple functionalities to
computationally study cell-cell communication from single-cell
transcriptomics data with complex multi-sample, multi-group design.**
The goal of this toolbox is to study differences in intercellular
communication between groups of samples of interest (eg patients of
different disease subtypes).

## Main functionalities of multinichenetr

Specific functionalities of this package include:

-   Finding differential expressed and active ligand-receptor
    interactions from different sender-receiver pairs between different
    sample groups

## Installation of multinichenetr

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc.

You can install multinichenetr (and required dependencies) from github
with:

    # install.packages("devtools")
    devtools::install_github("saeyslab/nichenetr")
    devtools::install_github("browaeysrobin/multinichenetr")

multinichenetr was tested on both Windows and Linux (most recently
tested R version: R 4.0.3)

## Learning to use multinichenetr

In the following vignettes, you can find how to do a multi-sample,
multi-group NicheNet analysis:

-   [Multi-Group Multi-Sample Cell-Cell Communication Analysis via
    NicheNet: HNSCC application](vignettes/basic_analysis.md):
    `vignette("basic_analysis", package="multinichenetr")`

## References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>