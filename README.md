<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# multinichenetr

<!-- badges: start -->

[![R build
status](https://github.com/browaeysrobin/multinichenetr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/browaeysrobin/multinichenetr/actions)
[![Coverage
Status](https://codecov.io/gh/browaeysrobin/multinichenetr/branch/master/graph/badge.svg?token=0X627I4TM7)](https://codecov.io/gh/browaeysrobin/multinichenetr)
<!-- badges: end -->

**multinichenetr: the R package containing multiple functionalities to
computationally study cell-cell communication from single-cell
transcriptomics data with complex multi-sample, multi-condition
designs.** The goal of this toolbox is to study differences in
intercellular communication between groups of samples of interest (eg
patients of different disease states).

The main goal of the MultiNicheNet package is to find which
ligand-receptor interactions are differentially expressed and
differentially active between conditions of interest, such as patient
groups. Compared to the normal NicheNet workflow, MultiNicheNet is more
suited to tackle complex experimental designs such as those with
multiple samples and conditions, and multiple receiver cell types of
interest.

In the MultiNicheNet approach, we allow the user to prioritize
differential cell-cell communication events (ligand-receptor
interactions and downstream signaling to target genes) based on the
following criteria:

-   Upregulation of the ligand in a sender cell type and/or upregulation
    of the receptor in a receiver cell type - in the condition of
    interest.
-   Sufficiently high expression levels of ligand and receptor in many
    samples of the same group (to mitigate the influence of outlier
    samples).
-   Cell-type and condition specific expression of the ligand in the
    sender cell type and receptor in the receiver cell type (to mitigate
    the influence of upregulated but still relatively weakly expressed
    ligands/receptors)
-   High NicheNet ligand activity, to further prioritize ligand-receptor
    pairs based on their predicted effect of the ligand-receptor
    interaction on the gene expression in the receiver cell type

MultiNicheNet combines all these criteria in a single prioritization
score, which is also comparable between all sender-receiver pairs. This
way, MultiNicheNet extends on the prioritization done by NicheNet, which
is only based on the ligand activity score.

Users can customize the weights of these different factors to prioritize
some of these criteria stronger, or neglect them altogether. If wanted,
the user can also take into account relative abundance of the sender and
receiver cell type in the condition of interest.

At the basis of MultiNicheNet for defining differentially expressed
ligands, receptors and target genes, is the the differential state
analysis as discussed by muscat, which provides a framework for
cell-level mixed models or methods based on aggregated “pseudobulk” data
(<https://doi.org/10.1038/s41467-020-19894-4>,
<https://bioconductor.org/packages/release/bioc/html/muscat.html>). We
use this muscat framework to make inferences on the sample-level (as
wanted in a multi-sample, multi-condition setting) and not the classic
cell-level differential expression analysis of Seurat
(Seurat::FindMarkers), because muscat allows us to overcome some of the
limitations of cell-level analyses for differential state analyses. Some
of these limitations include: a bias towards samples with more cells of
cell type, a lack of flexibility to work with complex study designs, and
a too optimistic estimation of the statistical power since the analysis
is done at the cell-level and not at the sample level.

In the future, we might extend the differential expression analyses
options to include other frameworks than muscat.

## Main functionalities of multinichenetr

-   Prioritizing the most important ligand-receptor interactions from
    different sender-receiver pairs between different sample groups,
    according to criteria such as condition specificity, cell-type
    specificity, signaling activity, and more.
-   Finding differential expressed ligand-receptor interactions from
    different sender-receiver pairs between different sample groups
    (Differential Ligand-Receptor network inference).
-   Predicting the most active ligand-receptor interactions in different
    sample groups based on predicted signaling effects (NicheNet ligand
    activity analysis).
-   Predicting specific downstream affected target genes of
    ligand-receptor links of interest (NicheNet ligand-target
    inference).

## Installation of multinichenetr

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc.

You can install multinichenetr (and required dependencies) from github
with:

    # install.packages("devtools")
    devtools::install_github("saeyslab/nichenetr")
    devtools::install_github("browaeysrobin/multinichenetr")

multinichenetr is tested via Github Actions version control on Windows,
Linux (Ubuntu) and Mac (most recently tested R version: R 4.0.4).

## Learning to use multinichenetr

In the following vignettes, you can find how to do a MultiNicheNet
analysis:

Performing an all-vs-all analysis with sender cell types and receiver
cell types in the same Seurat object:

-   [Multi-sample Multi-condition Cell-Cell Communication Analysis via
    NicheNet: HNSCC application; All-vs-All: step-by-step
    (recommended)](vignettes/basic_analysis_steps.md):
    `vignette("basic_analysis_steps", package="multinichenetr")`

-   [Multi-sample Multi-condition Cell-Cell Communication Analysis via
    NicheNet: HNSCC application;
    All-vs-All](vignettes/basic_analysis.md):
    `vignette("basic_analysis", package="multinichenetr")`

Performing an analysis with predefined sender and receiver cell types in
different Seurat objects:

-   [Multi-sample Multi-condition Cell-Cell Communication Analysis via
    NicheNet: HNSCC application; Separate: step-by-step
    (recommended)](vignettes/basic_analysis_steps_separate.md):
    `vignette("basic_analysis_separate_steps", package="multinichenetr")`

-   [Multi-sample Multi-condition Cell-Cell Communication Analysis via
    NicheNet: HNSCC application;
    Separate](vignettes/basic_analysis_separate.md):
    `vignette("basic_analysis_separate", package="multinichenetr")`

## References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>

Crowell, H.L., Soneson, C., Germain, PL. et al. muscat detects
subpopulation-specific state transitions from multi-sample
multi-condition single-cell transcriptomics data. Nat Commun 11, 6077
(2020). <https://doi.org/10.1038/s41467-020-19894-4>
