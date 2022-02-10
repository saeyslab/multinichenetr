## description of data

#' SingleCellExperiment object containing scRNAseq data (subsampled)
#'
#' SingleCellExperiment object containing scRNAseq data (subsampled). Source of the data: Puram et al., Cell 2017: “Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer.”. This example data was downsampled (features and cells). Interesting metadata/colData columns: tumor, pEMT, pEMT_fine, celltype, batch
#'
#' @format An object of class SingleCellExperiment
#'
"sce"
#' Gene annotation information: version 2 - january 2022 - suited for alias conversion
#'
#' A data.frame/tibble describing HGNC human gene symbols, their entrez ids and potential aliases.
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{symbol}{human gene symbol}
#'   \item{entrez}{human gene entrez}
#'   \item{alias}{human gene alias}
#'   }
#'
"geneinfo_alias_human"
#' Gene annotation information: version 2 - january 2022 - suited for alias conversion
#'
#' A data.frame/tibble describing MGI mouse gene symbols, their entrez ids and potential aliases.
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{symbol}{mouse gene symbol}
#'   \item{entrez}{mouse gene entrez}
#'   \item{alias}{mouse gene alias}
#'   }
#'
"geneinfo_alias_mouse"