<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# multinichenetr

<!-- badges: start -->

[![R build
status](https://github.com/saeyslab/multinichenetr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saeyslab/multinichenetr/actions)
[![Coverage
Status](https://codecov.io/gh/saeyslab/multinichenetr/branch/master/graph/badge.svg?token=0X627I4TM7)](https://codecov.io/gh/saeyslab/multinichenetr)
<!-- badges: end -->

**multinichenetr: the R package for differential cell-cell communication
analysis from single-cell transcriptomics data with complex
multi-sample, multi-condition designs.** The goal of this toolbox is to
study differences in intercellular communication between groups of
samples of interest (eg patients of different disease states).

You can read all about MultiNicheNet in the following preprint:
<https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1>

## Aims of MultiNicheNet

The main goal of the MultiNicheNet package is to find which
ligand-receptor interactions are differentially expressed and
differentially active between conditions of interest, such as patient
groups. Compared to the normal NicheNet workflow, MultiNicheNet 1)
considers more information to prioritize interactions (differential and
cell-type specific expression in addition to ligand activity); 2) is
more suited to tackle complex experimental designs such as those with
multiple samples and conditions, and multiple receiver cell types of
interest; and 3) enables incorporation of complementary data types for
prioritization of cell-cell interactions (e.g., proteomics).

<br><br> ![](vignettes/overview_figure.png) <br><br>

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
    pairs based on the enrichment of their predicted target genes in the
    receiver cell type

MultiNicheNet combines all these criteria in a single prioritization
score, which is also comparable between all sender-receiver pairs. This
way, MultiNicheNet extends on the prioritization done by NicheNet, which
is only based on the ligand activity score. Importantly, this
prioritization framework is flexible and extendable: it is possible to
design an additional criterion based on complementary data such as serum
proteomics.

To summarize the main functionalities of multinichenetr:

-   Prioritizing the most differential ligand-receptor interactions from
    different sender-receiver pairs between different sample groups,
    according to criteria such as condition specificity, cell-type
    specificity, ligand activity (= downstream signaling activity /
    target gene enrichement), and more.
-   Predicting the specific downstream affected target genes of
    prioritized differential ligand-receptor links of interest.
-   Predicting intercellular gene regulatory networks, connecting
    ligands to ligand- or receptor-encoding target genes in other cell
    types, enabling predictions concerning intercellular cascade and
    feedback mechanisms.

Note: at the basis of MultiNicheNet for defining differentially
expressed ligands, receptors and target genes, is the the differential
state analysis as discussed by **muscat**, which provides a framework
for cell-level mixed models or methods based on aggregated “pseudobulk”
data (<https://doi.org/10.1038/s41467-020-19894-4>,
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

# IMPORTANT UPDATE MAY 2024: multinichenetr 2.0.0.

Since May 13 2024, we have put out multinichenetr 2.0.0. on the main
branch. Please install this latest version of multinichenetr, which has
been improved at several levels compared to the first version:

### Novelties

Here is a concise list of methodological novelties compared to the
previous version. See the tutorials for more in-depth explanation on
each of these aspects. These novelties will also be discussed in the
revised version of the manuscript.

-   It is now possible to add new prioritization criteria based on
    additional data modalities. This option can helps users prioritize
    interactions better if they have additional complementary data, such
    as proteomics data. This option can help overcome major limitations
    of cell-cell communication inference from transcriptomics data.
-   An additional alternative workflow has been added to prioritize
    interactions involving condition-specific cell types. This option
    was added to not ignore biological signal from potentially relevant
    condition-specifc cell types.
-   More emphasis on the usability of the “Intercellular regulatory
    network” to prune prioritized interactions.
-   Integration with the Omnipath L-R database to check database source,
    quality and curation effort behind prioritized ligand-receptor
    interactions.

### Software changes

Here is a concise list of software changes compared to the previous
version. See the tutorials for more in-depth explanation on each of
these aspects.

-   Gene filtering has been modified slightly, resulting in fewer and
    easier to interpret parameters and stronger consistency with the
    rest of the MultiNicheNet methodology
-   Prioritization can now be done by only using the up-regulatory
    ligand activities instead of both up- and downregulatory activities.
    This option was added because of the difficulty of interpreting
    downregulatory activities.
-   The option to flexibly change prioritization weights has been
    replaced by the option to select biological scenario’s. This reduces
    the number of parameters for the end users and limits unwanted
    tunability of end results.  
-   The standard interpretable bubble plot visualization has been
    extended and provides now information about cell-type specificity,
    fraction of expression, and curation effort of the ligand-receptor
    pairs according to Omnipath. As a result of this, this plot now
    summarizes all the criteria used for prioritization and gives an
    indication about the curation effort of the ligand-receptor prior
    knowledge. This can help users now to get better insights in their
    results and define better candidate interactions for follow-up
    validation.

### Notes

-   Some visualization options of multinichenetr version &gt;= 2.0.0.
    are only compatible with MultiNicheNet output object generated by
    version multinichenetr version &gt;= 2.0.0. Please use
    multinichenetr version &gt;= 2.0.0. for buth running and
    interpreting the analysis.

### Call for feedback

multinichenetr is in ongoing development. We always appreciate it if you
let us know how we can improve the
software/documentation/algorithm/output visualizations further. The best
way to do this is through the issues page:
<https://github.com/saeyslab/multinichenetr/issues>

# How to install multinichenetr?

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc.

You can install multinichenetr (and required dependencies) from github
with:

    # install.packages("devtools")
    devtools::install_github("saeyslab/nichenetr")
    devtools::install_github("saeyslab/multinichenetr")

It is possible that during installation the following warning is thrown:

“glmmTMB was built with TMB version 1.9.4” “Current TMB version is 1.9.5

This warning can be safely ignored since this does not affect
multinichenetr’s installation and functionalities.

multinichenetr is tested via Github Actions version control on Windows,
Linux (Ubuntu) and Mac (most recently tested R version: R 4.3.1.).

# How to use multinichenetr?

We provide several vignettes demonstrating the different types of
analysis that can be performed with MultiNicheNet, and the several types
of downstream visualizations that can be created.

## Tutorials

We recommend users to start with the following vignette, which
demonstrates the different steps in the analysis and exploration of the
output. This is the recommended vignette to learn MultiNicheNet.

-   [**MultiNicheNet - comprehensive tutorial** - Condition A vs
    Condition B vs Condition
    C](vignettes/basic_analysis_steps_MISC.knit.md) | [*R Markdown
    version*](vignettes/basic_analysis_steps_MISC.Rmd) | [*HTML
    version*](vignettes/basic_analysis_steps_MISC.html)

That vignette provides an example of a comparison between 3 groups. The
following vignettes demonstrate how to analyze cell-cell communication
differences in other settings. These vignettes are the best vignettes to
learn how to apply MultiNicheNet to different datastes for addressing
different questions. To reduce the length of these vignettes, the
sections on downstream analysis has been reduced strongly and a wrapper
function is sometimes used to perform the core analysis. So we strongly
recommend to read these vignettes to learn how to perform the analysis
in different settings, but still perform all additional analyses and
checks as demonstrated in the comprehensive tutorial vignette.

-   [Condition A vs Condition B - without repeated
    subjects](vignettes/pairwise_analysis_MISC.knit.md) | [*R Markdown
    version*](vignettes/pairwise_analysis_MISC.Rmd) | [*HTML
    version*](vignettes/pairwise_analysis_MISC.html)
-   [Condition A vs Condition B - with **repeated subjects**: paired
    analysis with subject\_id as
    covariate](vignettes/paired_analysis_SCC.knit.md) | [*R Markdown
    version*](vignettes/paired_analysis_SCC.Rmd) | [*HTML
    version*](vignettes/paired_analysis_SCC.html)
-   [Condition A vs Condition B - from integrated atlas data: correcting
    for **batch
    effects**](vignettes/batch_correction_analysis_LungAtlas.knit.md) |
    [*R Markdown
    version*](vignettes/batch_correction_analysis_LungAtlas.Rmd) |
    [*HTML version*](vignettes/batch_correction_analysis_LungAtlas.html)
-   [**Multifactorial analysis**: condition-differences in cell-cell
    communication dynamics: (Time2\_ConditionA - Time1\_ConditionA) vs
    (Time2\_ConditionB -
    Time1\_ConditionB)](vignettes/multifactorial_analysis_BreastCancer.knit.md)
    | [*R Markdown
    version*](vignettes/multifactorial_analysis_BreastCancer.Rmd) |
    [*HTML
    version*](vignettes/multifactorial_analysis_BreastCancer.html)

In addition to these vignettes, we also provide 2 other vignettes
showcasing the flexibility and extendibility of the MultiNicheNet
framework: how to include additional complementary data modalities and
assess condition-specific cell types:

-   [Incorporating **additional prioritization criteria** based on
    **complementary data** types: serum
    **proteomics**](vignettes/add_proteomics_MISC.knit.md) | [*R
    Markdown version*](vignettes/add_proteomics_MISC.Rmd) | [*HTML
    version*](vignettes/add_proteomics_MISC.html)
-   [Altnernative workflow to assess interactions involving
    **condition-specific cell
    types**](vignettes/condition_specific_celltype_MISC.knit.md) | [*R
    Markdown version*](vignettes/condition_specific_celltype_MISC.Rmd) |
    [*HTML version*](vignettes/condition_specific_celltype_MISC.html)

When applying MultiNicheNet on datasets with many samples and cell
types, it is often needed to run the analysis on HPC infrastructure. In
those cases, we recommend to first set up your pipeline locally on a
subset of the data (eg subset of 2/3 cell types). Then we recommend to
run the core analysis on HPC infrastructure and save the output, and
finally interpret this output locally. In the following scripts you can
see an example of how we split up the analysis in two parts: 1) running
MultiNicheNet (with qsub on gridengine cluster) and saving necessary
output and plots; and 2) interpreting the results and generating
visualizations. These scripts are illustrative and will not work
directly when you would just run them.

-   [MultiNicheNet core analysis on gridengine cluster via
    HPC](vignettes/multinichenet_qsub.Rmd)
-   [Interpreting MultiNicheNet results
    locally](vignettes/multinichenet_intepretation.Rmd)

## Guidelines for parameter changes and interpretation of output figures

To help users in interpreting parameter values and output figures, we
provide the following two files:

-   [Parameter interpretation](parameter_interpretation.pdf): provides
    an explanation of different parameter choices - can help users in
    deciding when the default parameter values would not be optimal for
    their own dataset
-   [Output figure interpretation](output_interpretation.pdf): provides
    an explanation of different output figures - can help users in
    drawing hypotheses based on the output figures

## Frequently recurring questions and issues

-   The input data needed for MultiNicheNet should be raw counts, and
    metadata of cells giving information about the sample, condition and
    cell type. In all vignettes, we assume that the data has been
    preprocessed adequately (proper cell filtering, doublet removal,
    ambient RNA correction,…).
-   We strongly recommend having at least 4 samples in each of the
    groups/conditions you want to compare. With less samples, the
    benefits of performing a pseudobulk-based DE analysis are less clear
    and non-multi-sample tools for differential cell-cell communication
    might be better alternatives. If you want to perform differential
    cell-cell communication with a MultiNicheNet-like prioritization
    framework, you can have a look at this vignette: [Differential
    cell-cell Communication for datasets with limited samples:
    “sample-agnostic/cell-level”
    MultiNicheNet](vignettes/basic_analysis_steps_MISC_SACL.knit.md).
    Just realize that the analysis is based on a limited number of
    samples, and it will be hard to draw strong conclusions. This may
    often be the best you can get out of your data, but it is not a
    practice we would recommend.
-   Visualization functions of multinichenetr v.2.0.0 require output
    objects created by multinichenetr v.2.0.0

# References

Browaeys, R. et al. MultiNicheNet: a flexible framework for differential
cell-cell communication analysis from multi-sample multi-condition
single-cell transcriptomics data. (preprint).
<https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1>

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<a href="https://doi:10.1038/s41592-019-0667-5"
class="uri">https://doi:10.1038/s41592-019-0667-5</a>

Crowell, H.L., Soneson, C., Germain, PL. et al. muscat detects
subpopulation-specific state transitions from multi-sample
multi-condition single-cell transcriptomics data. Nat Commun 11, 6077
(2020). <https://doi.org/10.1038/s41467-020-19894-4>
