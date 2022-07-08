#' @title multi_nichenet_analysis
#'
#' @description \code{multi_nichenet_analysis}  Perform a MultiNicheNet analysis. See `multi_nichenet_analysis_separate` and `multi_nichenet_analysis_combined` for more information.
#' @usage multi_nichenet_analysis(sender_receiver_separate = TRUE, ...)
#'
#' @param sender_receiver_separate Indicates whether the user gives as input one separate SingleCellExperiment object with sender cell types and one with receiver cell types (TRUE) or whether only one SingleCellExperiment object with both sender and receiver cell types of interest (FALSE).
#' TRUE calls the function `multi_nichenet_analysis_separate`, FALSE calls the function `multi_nichenet_analysis_combined`. Default: TRUE.
#' @param ... Arguments to `multi_nichenet_analysis_separate` or `multi_nichenet_analysis_combined`. 
#'
#' @return  List containing different types of information and output of the MultiNicheNet analysis. 
#' See `multi_nichenet_analysis_separate` and `multi_nichenet_analysis_combined` for more information.
#'
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' batches = NA
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      covariates = covariates,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl, 
#'      sender_receiver_separate = FALSE
#'      )
#' }
#'
#' @export
#'
multi_nichenet_analysis = function(sender_receiver_separate = TRUE, ...){

  if(!is.logical(sender_receiver_separate)) {
    stop("The sender_receiver_separate argument should be TRUE or FALSE")
  }

  if(sender_receiver_separate == TRUE){
    output = multi_nichenet_analysis_separate(...)
  } else {
    output = multi_nichenet_analysis_combined(...)
  }
  return(output)
}

#' @title multi_nichenet_analysis_separate
#'
#' @description \code{multi_nichenet_analysis_separate}  Perform a MultiNicheNet analysis between sender cell types and receiver cell types of interest.
#' @usage multi_nichenet_analysis_separate(
#' sce_receiver, sce_sender,celltype_id_receiver,celltype_id_sender,sample_id,group_id, batches, covariates, lr_network,ligand_target_matrix,contrasts_oi,contrast_tbl, fraction_cutoff = 0.05,
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0),
#' assay_oi_pb ="counts",fun_oi_pb = "sum",de_method_oi = "edgeR",min_cells = 10,logFC_threshold = 0.25,p_val_threshold = 0.05, p_val_adj = FALSE, empirical_pval = TRUE, top_n_target = 250, verbose = FALSE, n.cores = 1, return_lr_prod_matrix = FALSE, findMarkers = FALSE, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7, top_n_LR = 2500)
#'
#' @param sce_receiver SingleCellExperiment object containing the receiver cell types of interest
#' @param sce_sender SingleCellExperiment object containing the sender cell types of interest
#' @param celltype_id_receiver Name of the meta data column that indicates the cell type of a cell (in sce_receiver).
#' @param celltype_id_sender Name of the meta data column that indicates the cell type of a cell (in sce_sender). 
#' @param sample_id Name of the meta data column that indicates from which sample/patient a cell comes from (in both sce_receiver and sce_sender)
#' @param group_id Name of the meta data column that indicates from which group/condition a cell comes from (in both sce_receiver and sce_sender)
#' @param batches NA if no batches should be corrected for. If there should be corrected for batches during DE analysis and pseudobulk expression calculation, this argument should be the name(s) of the columns in the meta data that indicate the batch(s). Should be categorical. Pseudobulk expression values will be corrected for the first element of this vector. 
#' @param covariates NA if no covariates should be corrected for. If there should be corrected for covariates uring DE analysis, this argument should be the name(s) of the columns in the meta data that indicate the covariate(s). Can both be categorical and continuous. Pseudobulk expression values will not be corrected for the first element of this vector. 
#' @param lr_network Prior knowledge Ligand-Receptor network (columns: ligand, receptor)
#' @param ligand_target_matrix Prior knowledge model of ligand-target regulatory potential (matrix with ligands in columns and targets in rows). See https://github.com/saeyslab/nichenetr.
#' @param contrasts_oi String indicating the contrasts of interest (= which groups/conditions will be compared) for the differential expression and MultiNicheNet analysis. 
#' We will demonstrate here a few examples to indicate how to write this. Check the limma package manuals for more information about defining design matrices and contrasts for differential expression analysis. \cr
#' If wanting to compare group A vs B: `contrasts_oi = c("'A-B'")` \cr
#' If wanting to compare group A vs B & B vs A: `contrasts_oi = c("'A-B','B-A'")` \cr
#' If wanting to compare group A vs B & A vs C & A vs D: `contrasts_oi = c("'A-B','A-C', 'A-D'")` \cr
#' If wanting to compare group A vs B and C: `contrasts_oi = c("'A-(B+C)/2'")` \cr
#' If wanting to compare group A vs B, C and D: `contrasts_oi = c("'A-(B+C+D)/3'")` \cr
#' If wanting to compare group A vs B, C and D & B vs A,C,D: `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")` \cr
#' Note that the groups A, B, ... should be present in the meta data column 'group_id'.
#' @param contrast_tbl Data frame providing names for each of the contrasts in contrasts_oi in the 'contrast' column, and the corresponding group of interest in the 'group' column. Entries in the 'group' column should thus be present in the group_id column in the metadata. 
#' Example for `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")`:
#' `contrast_tbl = tibble(contrast = c("A-(B+C+D)/3","B-(A+C+D)/3"), group = c("A","B"))`
#' @param fraction_cutoff Cutoff indicating the minimum fraction of cells of a cell type in a specific sample that are necessary to consider the gene as expressed. 
#' @param prioritizing_weights Named vector indicating the relative weights of each prioritization criterion included in MultiNicheNet.
#' Default: prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0)
#' Details about the meaning of this naming: \cr
#' de_ligand: importance of DE of the ligand in a certain sender (indicates upregulation of the ligand in condition of interest compared to other conditions): based on a scaling of the following calculation: -log10(p-value) and logFC. \cr
#' de_receptor: importance of DE of the receptor in a certain receiver (indicates upregulation of the receptor in condition of interest compared to other conditions): : based on a scaling of the following calculation: -log10(p-value) and logFC.\cr
#' activity_scaled: importance of the scaled ligand activity of the ligand in a certain receiver-condition combination (indicates active signaling of a ligand in a receiver cell type in a certain condition compared to other conditions). Scaled activity: indicates the ranking of ligands in a certain receiver-condition combination. \cr
#' exprs_ligand: importance of the condition-and-sender specific expression of a ligand (taking into account average expression value and fraction of cells expressing a ligand).  \cr
#' exprs_receptor: importance of the condition-and-receiver specific expression of a receptor (taking into account average expression value and fraction of cells expressing a receptor). \cr
#' frac_exprs_ligand_receptor: importance of the fraction of samples in a group that show high enough expression of the specific ligand-receptor pair. \cr
#' abund_sender: importance of relative cell type abundance of the sender cell type. \cr
#' abund_receiver: importance of relative cell type abundance of the receiver cell type. \cr
#' @param assay_oi_pb Indicates which information of the assay of interest should be used (counts, scaled data,...). Default: "counts". See `muscat::aggregateData`.
#' @param fun_oi_pb Indicates way of doing the pseudobulking. Default: "sum". See `muscat::aggregateData`.
#' @param de_method_oi Indicates the DE method that will be used after pseudobulking. Default: "edgeR". See `muscat::pbDS`.
#' @param min_cells Indicates the minimal number of cells that a sample should have to be considered in the DE analysis. Default: 10. See `muscat::pbDS`.
#' @param logFC_threshold For defining the gene set of interest for NicheNet ligand activity: what is the minimum logFC a gene should have to belong to this gene set? Default: 0.25/
#' @param p_val_threshold For defining the gene set of interest for NicheNet ligand activity: what is the maximam p-value a gene should have to belong to this gene set? Default: 0.05.
#' @param p_val_adj For defining the gene set of interest for NicheNet ligand activity: should we look at the p-value corrected for multiple testing? Default: FALSE.
#' @param empirical_pval For defining the gene set of interest for NicheNet ligand activity - and for ranking DE ligands and receptors: should we use the normal p-values, or the p-values that are corrected by the empirical null procedure. The latter could be beneficial if p-value distribution histograms indicate potential problems in the model definition (eg not all relevant batches are selected, etc). Default: TRUE.
#' @param top_n_target For defining NicheNet ligand-target links: which top N predicted target genes. See `nichenetr::get_weighted_ligand_target_links()`.
#' @param verbose Indicate which different steps of the pipeline are running or not. Default: FALSE.
#' @param n.cores The number of cores used for parallel computation of the ligand activities per receiver cell type. Default: 1 - no parallel computation.
#' @param return_lr_prod_matrix Indicate whether to calculate a senderLigand-receiverReceptor matrix, which could be used for unsupervised analysis of the cell-cell communication. Default FALSE. Setting to FALSE might be beneficial to avoid memory issues.
#' @param findMarkers Indicate whether we should also calculate DE results with the classic scran::findMarkers approach. Default (recommended): FALSE.
#' @param filterByExpr.min.count check edgeR::filterByExpr documentation - Default = 7. Increase this if you want more stringent filtering in terms of required counts per gene per sample.
#' @param filterByExpr.min.total.count check edgeR::filterByExpr documentation -Default = 15. Increase this if you want more stringent filtering in terms of required counts per gene per sample.
#' @param filterByExpr.large.n check edgeR::filterByExpr documentation -Default = 4. Increase this if you want more stringent filtering in terms of required samples with expression of the gene.
#' @param filterByExpr.min.prop check edgeR::filterByExpr documentation - Default = 0.7. Increase this if you want more stringent filtering in terms of required samples with expression of the gene.
#' @param  top_n_LR top nr of LR pairs for which correlation with target genes will be calculated. Is 2500 by default. If you want to calculate correlation for all expressed LR pairs, set this argument to NA.


#' @return List containing information and output of the MultiNicheNet analysis.\cr
#' celltype_info: contains average expression value and fraction of each cell type - sample combination, \cr
#' celltype_de: contains output of the differential expression analysis,  \cr
#' sender_receiver_info: links the expression information of the ligand in the sender cell types to the expression of the receptor in the receiver cell types,  \cr
#' sender_receiver_de: links the differential information of the ligand in the sender cell types to the expression of the receptor in the receiver cell types \cr
#' ligand_activities_targets_DEgenes: contains the output of the NicheNet ligand activity analysis, and the NicheNet ligand-target inference \cr
#' prioritization_tables: contains the tables with the final prioritization scores \cr
#' lr_prod_mat: matrix of the ligand-receptor expression product of the expressed senderLigand-receiverReceptor pairs, \cr
#' grouping_tbl: data frame showing the group per sample  \cr
#' lr_target_prior_cor: data frame showing the expression correlation between ligand-receptor pairs and DE genes + NicheNet regulatory potential scores indicating the amount of prior knowledge supporting a LR-target regulatory link
#' @import dplyr
#' @import ggplot2
#' @importFrom generics setdiff intersect union
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidyr spread
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sce_receiver = sce %>% subset(subset = celltype == "Malignant")
#' sce_sender = sce %>% subset(subset = celltype == "CAF")
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id_receiver = "celltype"
#' celltype_id_sender = "celltype"
#' batches = NA
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis_separate(
#'      sce_receiver = sce_receiver, 
#'      sce_sender = sce_sender,
#'      celltype_id_receiver = celltype_id_receiver, 
#'      celltype_id_sender = celltype_id_sender,
#'      sample_id = sample_id, 
#'      group_id = group_id, 
#'      batches = batches,
#'      covariates = covariates,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#'}
#'
#' @export
#'
multi_nichenet_analysis_separate = function(sce_receiver,
                                            sce_sender,
                                            celltype_id_receiver,
                                            celltype_id_sender,
                                            sample_id,
                                            group_id,
                                            batches,
                                            covariates,
                                            lr_network,
                                            ligand_target_matrix,
                                            contrasts_oi,
                                            contrast_tbl,
                                            fraction_cutoff = 0.05,
                                            prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0),
                                            assay_oi_pb ="counts",
                                            fun_oi_pb = "sum",
                                            de_method_oi = "edgeR",
                                            min_cells = 10,
                                            logFC_threshold = 0.25,
                                            p_val_threshold = 0.05,
                                            p_val_adj = FALSE,
                                            empirical_pval = TRUE,
                                            top_n_target = 250, verbose = FALSE, n.cores = 1, return_lr_prod_matrix = FALSE, findMarkers = FALSE, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7, top_n_LR = 2500){


  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # input checks
  if (class(sce_receiver) != "SingleCellExperiment") {
    stop("sce_receiver should be a SingleCellExperiment object")
  }
  if (class(sce_sender) != "SingleCellExperiment") {
    stop("sce_sender should be a SingleCellExperiment object")
  }
  if (!celltype_id_receiver %in% colnames(SummarizedExperiment::colData(sce_receiver))) {
    stop("celltype_id_receiver should be a column name in the metadata dataframe of sce_receiver")
  }
  if (celltype_id_receiver != make.names(celltype_id_receiver)) {
    stop("celltype_id_receiver should be a syntactically valid R name - check make.names")
  }
  if (!celltype_id_sender %in% colnames(SummarizedExperiment::colData(sce_sender))) {
    stop("celltype_id_sender should be a column name in the metadata dataframe of sce_sender")
  }
  if (celltype_id_sender != make.names(celltype_id_sender)) {
    stop("celltype_id_sender should be a syntactically valid R name - check make.names")
  }
  if (!sample_id %in% colnames(SummarizedExperiment::colData(sce_receiver))) {
    stop("sample_id should be a column name in the metadata dataframe of sce_receiver")
  }
  if (!sample_id %in% colnames(SummarizedExperiment::colData(sce_sender))) {
    stop("sample_id should be a column name in the metadata dataframe of sce_sender")
  }
  if (sample_id != make.names(sample_id)) {
    stop("sample_id should be a syntactically valid R name - check make.names")
  }
  if (!group_id %in% colnames(SummarizedExperiment::colData(sce_receiver))) {
    stop("group_id should be a column name in the metadata dataframe of sce_receiver")
  }
  if (!group_id %in% colnames(SummarizedExperiment::colData(sce_sender))) {
    stop("group_id should be a column name in the metadata dataframe of sce_sender")
  }
  if (group_id != make.names(group_id)) {
    stop("group_id should be a syntactically valid R name - check make.names")
  }
  
  # group, sample and cell type id must not be numeric
  if(is.double(SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver])){
    stop("SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce_receiver)[,group_id])){
    stop("SummarizedExperiment::colData(sce_receiver)[,group_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce_receiver)[,sample_id])){
    stop("SummarizedExperiment::colData(sce_receiver)[,sample_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce_sender)[,celltype_id_sender])){
    stop("SummarizedExperiment::colData(sce_sender)[,celltype_id_sender] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce_sender)[,group_id])){
    stop("SummarizedExperiment::colData(sce_sender)[,group_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce_sender)[,sample_id])){
    stop("SummarizedExperiment::colData(sce_sender)[,sample_id] should be a character vector or a factor")
  }
  
  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce_sender)[,celltype_id_sender])){
    is_make_names = levels(SummarizedExperiment::colData(sce_sender)[,celltype_id_sender]) == make.names(levels(SummarizedExperiment::colData(sce_sender)[,celltype_id_sender]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_sender)[,celltype_id_sender]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_sender)[,celltype_id_sender] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce_sender)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce_sender)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce_sender)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_sender)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_sender)[,group_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce_sender)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce_sender)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce_sender)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_sender)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_sender)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver])){
    is_make_names = levels(SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver]) == make.names(levels(SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce_receiver)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce_receiver)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce_receiver)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_receiver)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_receiver)[,group_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce_receiver)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce_receiver)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce_receiver)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce_receiver)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce_receiver)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  }
  
 
  if(!is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector")
  }
  if(!is.data.frame(contrast_tbl)){
    stop("contrast_tbl should be a data frame / tibble")
  }
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = SummarizedExperiment::colData(sce_receiver)[,group_id] %>% unique()

  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    # stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",","," ,", ", ")) %>% unlist() %>% unique()
  conditions_oi = conditions_oi[is.na(suppressWarnings(as.numeric(conditions_oi)))]
  
  if(length(contrasts_oi) != 1 | !is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }
  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_oi_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts_oi should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  if (sum(contrasts_oi_simplified %in% unique(contrast_tbl$contrast)) != length(contrasts_oi_simplified)) {
    stop("conditions written in contrasts_oi should be in the contrast column of contrast_tbl column! This is not the case, which can lead to errors downstream.")
  }

  
  # 
  groups_oi_contrast_tbl = contrast_tbl$group %>% unique()
  if(sum(groups_oi_contrast_tbl %in% groups_oi) != length(groups_oi_contrast_tbl)){
    stop("You have defined some groups in contrast_tbl$group that are not present SummarizedExperiment::colData(sce)[,group_id]. This will result in lack of information downstream. We recommend to change your metadata or this contrast_tbl appropriately.")
  }
  if(length(groups_oi_contrast_tbl) != length(contrast_tbl$group)){
    warning("According to your contrast_tbl, some of your contrasts will be assigned to the same group. This should not be a problem if this was intended, but be aware not to make mistakes in the further interpretation and plotting of the results.")
  }
  
  
  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce_receiver))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce_receiver")
    }
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce_sender))) != length(batches) ) {
      stop("batches should be NA or all present column name(s) in the metadata dataframe of sce_sender")
    }
  }
  
  # Check concordance ligand-receptor network and ligand-target network - and concordance with sce object features
  if (!is.matrix(ligand_target_matrix)){
    stop("ligand_target_matrix should be a matrix")
  }
  if (!is.data.frame(lr_network)){
    stop("lr_network should be a data frame / tibble")
  }

  if(! "ligand" %in% colnames(lr_network)){
    if("from" %in% colnames(lr_network)){
      lr_network = lr_network %>% dplyr::rename(ligand = from)
    } else {
      stop("The ligand-receptor network should have the columns ligand and receptor (or: from and to)")
    }
  }
  if(! "receptor" %in% colnames(lr_network)){
    if("to" %in% colnames(lr_network)){
      lr_network = lr_network %>% dplyr::rename(receptor = to)
    } else {
      stop("The ligand-receptor network should have the columns ligand and receptor (or: from and to)")
    }
  }
  
  ligands_lrnetwork = lr_network$ligand %>% unique()
  receptors_lrnetwork = lr_network$receptor %>% unique()
  ligands_ligand_target_matrix = colnames(ligand_target_matrix)

  if(length(ligands_ligand_target_matrix) < length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }
  if(length(ligands_lrnetwork) < length(ligands_ligand_target_matrix)){
    warning("Not all Ligands from your ligand-target matrix are in the ligand-receptor network")
  }

  if(length(rownames(sce_sender) %>% generics::intersect(ligands_lrnetwork)) < 25 ){
    warning("Less than 25 ligands from your ligand-receptor network are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(rownames(sce_sender) %>% generics::intersect(ligands_ligand_target_matrix)) < 25 ){
    warning("Less than 25 ligands from your ligand-target matrix are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(rownames(sce_receiver) %>% generics::intersect(receptors_lrnetwork)) < 25 ){
    warning("Less than 25 receptors from your ligand-receptor network are in your expression matrix of the receiver cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }

  if(length(prioritizing_weights) != 8 | !is.double(prioritizing_weights)) {
    stop("prioritizing_weights should be a numeric vector with length 8")
  }
  names_prioritizing_weights = c("de_ligand",
                                 "de_receptor",
                                 "activity_scaled",
                                 "exprs_ligand",
                                 "exprs_receptor",
                                 "frac_exprs_ligand_receptor",
                                 "abund_sender",
                                 "abund_receiver")
  if(sum(names_prioritizing_weights %in% names(prioritizing_weights)) != length(names_prioritizing_weights)) {
    stop("prioritizing_weights should be have the correct names. Check the vignettes and code documentation")
  }
  
  if(!is.character(assay_oi_pb)){
    stop("assay_oi_pb should be a character vector")
  } else {
    if(assay_oi_pb != "counts"){
      warning("are you sure you don't want to use the counts assay?")
    }
  }
  if(!is.character(fun_oi_pb)){
    stop("fun_oi_pb should be a character vector")
  }
  if(!is.character(de_method_oi)){
    stop("de_method_oi should be a character vector")
  }

  if(!is.double(min_cells)){
    stop("min_cells should be numeric")
  } else {
    if(min_cells <= 0) {
      warning("min_cells is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(!is.double(logFC_threshold)){
    stop("logFC_threshold should be numeric")
  } else {
    if(logFC_threshold <= 0) {
      warning("logFC_threshold is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(!is.double(p_val_threshold)){
    stop("p_val_threshold should be numeric")
  } else {
    if(p_val_threshold <= 0 | p_val_threshold > 1) {
      warning("p_val_threshold is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.10, 0 excluded.")
    }
  }
  if(!is.double(fraction_cutoff)){
    stop("fraction_cutoff should be numeric")
  } else {
    if(fraction_cutoff <= 0 | fraction_cutoff > 1) {
      stop("fraction_cutoff is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.25, 0 excluded.")
    }
  }
  if(!is.double(top_n_target)){
    stop("top_n_target should be numeric")
  } else {
    if(top_n_target <= 0 ) {
      warning("top_n_target is now 0 or smaller. We recommend having a positive, non-zero value for this parameter.")
    }
  }
  if(!is.logical(p_val_adj)){
    stop("p_val_adj should be TRUE or FALSE")
  }
  if(!is.logical(verbose)){
    stop("verbose should be TRUE or FALSE")
  }
  if(!is.logical(empirical_pval)){
    stop("empirical_pval should be TRUE or FALSE")
  }
  if(!is.logical(return_lr_prod_matrix)){
    stop("return_lr_prod_matrix should be TRUE or FALSE")
  }
  
  if(!is.double(n.cores)){
    stop("n.cores should be numeric")
  } else {
    if(n.cores <= 0 ) {
      warning("n.cores is now 0 or smaller. We recommend having a positive, non-zero value for this parameter.")
    }
  }
  
  if(verbose == TRUE){
    print("Make diagnostic cell type abundance plots")
  }
  
  
  ### Perform the DE analysis ----------------------------------------------------------------
  
  if(verbose == TRUE){
    print("Calculate differential expression for all cell types")
  }
  if(findMarkers == FALSE){
    DE_info_receiver = get_DE_info(sce = sce_receiver, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id_receiver, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells, 
                                   assay_oi_pb = assay_oi_pb,
                                   fun_oi_pb = fun_oi_pb,
                                   de_method_oi = de_method_oi, 
                                   findMarkers = findMarkers, 
                                   filterByExpr.min.count = filterByExpr.min.count, 
                                   filterByExpr.min.total.count = filterByExpr.min.total.count, 
                                   filterByExpr.large.n = filterByExpr.large.n, 
                                   filterByExpr.min.prop = filterByExpr.min.prop)
    
    DE_info_sender = get_DE_info(sce = sce_sender, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id_sender, batches = batches, covariates = covariates , contrasts_oi = contrasts_oi, min_cells = min_cells, 
                                 assay_oi_pb = assay_oi_pb,
                                 fun_oi_pb = fun_oi_pb,
                                 de_method_oi = de_method_oi,
                                 findMarkers = findMarkers,
                                 filterByExpr.min.count = filterByExpr.min.count, 
                                 filterByExpr.min.total.count = filterByExpr.min.total.count, 
                                 filterByExpr.large.n = filterByExpr.large.n, 
                                 filterByExpr.min.prop = filterByExpr.min.prop)
  } else {
    DE_info_receiver = get_DE_info(sce = sce_receiver, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id_receiver, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells, 
                                   assay_oi_pb = assay_oi_pb,
                                   fun_oi_pb = fun_oi_pb,
                                   de_method_oi = de_method_oi, 
                                   findMarkers = findMarkers,
                                   contrast_tbl = contrast_tbl)
    
    DE_info_sender = get_DE_info(sce = sce_sender, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id_sender, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells, 
                                 assay_oi_pb = assay_oi_pb,
                                 fun_oi_pb = fun_oi_pb,
                                 de_method_oi = de_method_oi,
                                 findMarkers = findMarkers,
                                 contrast_tbl = contrast_tbl)

  }

  empirical_pval_receiver = empirical_pval
  empirical_pval_sender = empirical_pval
    
  if(empirical_pval_receiver == TRUE){
    if(findMarkers == TRUE){
      DE_info_emp_receiver = get_empirical_pvals(DE_info_receiver$celltype_de_findmarkers)
    } else {
      DE_info_emp_receiver = get_empirical_pvals(DE_info_receiver$celltype_de$de_output_tidy)
    }
  } 
  if(empirical_pval_sender == TRUE){
    if(findMarkers == TRUE){
      DE_info_emp_sender = get_empirical_pvals(DE_info_sender$celltype_de_findmarkers)
    } else {
      DE_info_emp_sender = get_empirical_pvals(DE_info_sender$celltype_de$de_output_tidy)
    }
  } 
  
  if(empirical_pval_receiver == FALSE){
    if(findMarkers == TRUE){
      celltype_de_receiver = DE_info_receiver$celltype_de_findmarkers
    } else {
      celltype_de_receiver = DE_info_receiver$celltype_de$de_output_tidy
    }
  } else {
    celltype_de_receiver = DE_info_emp_receiver$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
  }
  if(empirical_pval_sender == FALSE){
    if(findMarkers == TRUE){
      celltype_de_sender = DE_info_sender$celltype_de_findmarkers
    } else {
      celltype_de_sender = DE_info_sender$celltype_de$de_output_tidy
    }
  } else {
    celltype_de_sender = DE_info_emp_sender$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
  }
  
  
  ### Define Senders-Receivers
  senders_oi = celltype_de_sender$cluster_id %>% unique()
  receivers_oi = celltype_de_receiver$cluster_id %>% unique()
  genes_oi_sender = celltype_de_sender$gene %>% unique()
  genes_oi_receiver = celltype_de_receiver$gene %>% unique()
  
  sce_receiver = sce_receiver[genes_oi_receiver, SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver] %in% c(receivers_oi)]
  sce_sender = sce_sender[genes_oi_sender, SummarizedExperiment::colData(sce_sender)[,celltype_id_sender] %in% c(senders_oi)]
  
  sender_receiver_de = combine_sender_receiver_de(
    sender_de = celltype_de_sender,
    receiver_de = celltype_de_receiver,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network
  )
  
  ### Receiver abundance plots + Calculate expression information
  if(verbose == TRUE){
    print("Make diagnostic abundance plots + Calculate expression information")
  }
  
  abundance_expression_info = get_abundance_expression_info_separate(sce_receiver = sce_receiver, sce_sender = sce_sender, sample_id = sample_id, group_id = group_id, celltype_id_receiver = celltype_id_receiver, celltype_id_sender = celltype_id_sender, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches, min_cells = min_cells)

  
  ### Use the DE analysis for defining the DE genes in the receiver cell type and perform NicheNet ligand activity and ligand-target inference ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Calculate NicheNet ligand activities and ligand-target links")
  }
  ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_receiver,
    receivers_oi = receivers_oi,
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )))
  
  ### Combine the three types of information calculated above to prioritize ligand-receptor interactions ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Combine all the information in prioritization tables")
  }
  ### Remove types of information that we don't need anymore:
  
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  metadata_combined = SummarizedExperiment::colData(sce_receiver) %>% tibble::as_tibble()
  
  if(!is.na(batches)){
    grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group",batches)
  } else {
    grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group")
  }
  rm(sce_sender)
  rm(sce_receiver)
  prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    prioritizing_weights = prioritizing_weights,
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender
  ))
  
  # Prepare Unsupervised analysis of samples! ------------------------------------------------------------------------------------------------------------
  
  if(return_lr_prod_matrix == TRUE){
    
    ids_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% dplyr::pull(id) %>% unique()
    
    lr_prod_df = abundance_expression_info$sender_receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sample, id, ligand_receptor_pb_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(id, ligand_receptor_pb_prod)
    lr_prod_mat = lr_prod_df %>% dplyr::select(-sample) %>% data.frame() %>% as.matrix()
    rownames(lr_prod_mat) = lr_prod_df$sample
    
    col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
    row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
    
    lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(col_remove),colnames(.) %>% generics::setdiff(col_remove)]
  } else {
    lr_prod_mat = NULL
  }
  
  
  # Add information on prior knowledge and expression correlation between LR and target expression ------------------------------------------------------------------------------------------------------------
  lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de_receiver, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj, top_n_LR = top_n_LR)
  
  
  multinichenet_output = list(
    receiver_info = abundance_expression_info$receiver_info,
    receiver_de = celltype_de_receiver,
    sender_info = abundance_expression_info$sender_info,
    sender_de = celltype_de_sender,
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver, 
    abundance_data_sender = abundance_expression_info$abundance_data_sender, 
    # sender_receiver_info = abundance_expression_info$sender_receiver_info,
    # sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    lr_prod_mat = lr_prod_mat,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  )
  
  multinichenet_output = multinichenet_output %>% make_lite_output(top_n_LR = top_n_LR)
  
  return(multinichenet_output)
  
}
#' @title multi_nichenet_analysis_combined
#'
#' @description \code{multi_nichenet_analysis_combined}  Perform a MultiNicheNet analysis in an all-vs-all setting: all cell types in the data will be considered both as sender and receiver.
#' @usage multi_nichenet_analysis_combined(
#' sce, celltype_id, sample_id,group_id, batches, covariates, lr_network,ligand_target_matrix,contrasts_oi,contrast_tbl, senders_oi = NULL,receivers_oi = NULL, fraction_cutoff = 0.05,
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0),
#' assay_oi_pb ="counts",fun_oi_pb = "sum",de_method_oi = "edgeR",min_cells = 10,logFC_threshold = 0.25,p_val_threshold = 0.05,p_val_adj = FALSE, empirical_pval = TRUE, top_n_target = 250, verbose = FALSE, n.cores = 1, return_lr_prod_matrix = FALSE, findMarkers = FALSE, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7, top_n_LR = 2500)
#'
#' @param sce SingleCellExperiment object of the scRNAseq data of interest. Contains both sender and receiver cell types.
#' @param celltype_id Name of the column in the meta data of sce that indicates the cell type of a cell.
#' @param senders_oi Default NULL: all celltypes will be considered as senders If you want to select specific senders_oi: you can add this here as character vector.
#' @param receivers_oi Default NULL: all celltypes will be considered as receivers. If you want to select specific receivers_oi: you can add this here as character vector.
#' @inheritParams multi_nichenet_analysis_separate
#'
#'
#' @return List containing information and output of the MultiNicheNet analysis.
#' celltype_info: contains average expression value and fraction of each cell type - sample combination,
#' celltype_de: contains output of the differential expression analysis, 
#' sender_receiver_info: links the expression information of the ligand in the sender cell types to the expression of the receptor in the receiver cell types, 
#' sender_receiver_de: links the differential information of the ligand in the sender cell types to the expression of the receptor in the receiver cell types
#' ligand_activities_targets_DEgenes: contains the output of the NicheNet ligand activity analysis, and the NicheNet ligand-target inference
#' prioritization_tables: contains the tables with the final prioritization scores
#' lr_prod_mat: matrix of the ligand-receptor expression product of the expressed senderLigand-receiverReceptor pairs,
#' grouping_tbl: data frame showing the group per sample 
#' lr_target_prior_cor: data frame showing the expression correlation between ligand-receptor pairs and DE genes + NicheNet regulatory potential scores indicating the amount of prior knowledge supporting a LR-target regulatory link
#' 
#' @import dplyr
#' @import ggplot2
#' @importFrom generics setdiff intersect union
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidyr spread
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' batches = NA
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis_combined(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id, 
#'      batches = batches,
#'      covariates = covariates,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl)
#' }
#'
#' @export
#'
multi_nichenet_analysis_combined = function(sce,
                                            celltype_id,
                                            sample_id,
                                            group_id,
                                            batches,
                                            covariates,
                                            lr_network,
                                            ligand_target_matrix,
                                            contrasts_oi,
                                            contrast_tbl,
                                            senders_oi = NULL,
                                            receivers_oi = NULL,
                                            fraction_cutoff = 0.05,
                                            prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0),
                                            assay_oi_pb ="counts",
                                            fun_oi_pb = "sum",
                                            de_method_oi = "edgeR",
                                            min_cells = 10,
                                            logFC_threshold = 0.25,
                                            p_val_threshold = 0.05,
                                            p_val_adj = FALSE,
                                            empirical_pval = TRUE,
                                            top_n_target = 250, verbose = FALSE, n.cores = 1, return_lr_prod_matrix = FALSE, findMarkers = FALSE, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7, top_n_LR = 2500){


  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # input checks
  if (class(sce) != "SingleCellExperiment") {
    stop("sce should be a SingleCellExperiment object")
  }
  
  if (!celltype_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("celltype_id should be a column name in the metadata dataframe of sce")
  }
  if (celltype_id != make.names(celltype_id)) {
    stop("celltype_id should be a syntactically valid R name - check make.names")
  }
  if (!sample_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("sample_id should be a column name in the metadata dataframe of sce")
  }
  if (sample_id != make.names(sample_id)) {
    stop("sample_id should be a syntactically valid R name - check make.names")
  }
  if (!group_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("group_id should be a column name in the metadata dataframe of sce")
  }
  if (group_id != make.names(group_id)) {
    stop("group_id should be a syntactically valid R name - check make.names")
  }
  
  if(is.double(SummarizedExperiment::colData(sce)[,celltype_id])){
    stop("SummarizedExperiment::colData(sce)[,celltype_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,group_id])){
    stop("SummarizedExperiment::colData(sce)[,group_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,sample_id])){
    stop("SummarizedExperiment::colData(sce)[,sample_id] should be a character vector or a factor")
  }

  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce)[,celltype_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,celltype_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,celltype_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,celltype_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  }
  
  
  if(!is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector")
  }
  if(!is.data.frame(contrast_tbl)){
    stop("contrast_tbl should be a data frame / tibble")
  }
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = SummarizedExperiment::colData(sce)[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    # stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",","," ,", ", ")) %>% unlist() %>% unique()
  conditions_oi = conditions_oi[is.na(suppressWarnings(as.numeric(conditions_oi)))]

  if(length(contrasts_oi) != 1 | !is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }
  
  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_oi_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts_oi should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  if (sum(contrasts_oi_simplified %in% unique(contrast_tbl$contrast)) != length(contrasts_oi_simplified)) {
    stop("conditions written in contrasts_oi should be in the contrast column of contrast_tbl column! This is not the case, which can lead to errors downstream.")
  }
  
  #
  groups_oi_contrast_tbl = contrast_tbl$group %>% unique()
  if(sum(groups_oi_contrast_tbl %in% groups_oi) != length(groups_oi_contrast_tbl)){
    stop("You have defined some groups in contrast_tbl$group that are not present SummarizedExperiment::colData(sce)[,group_id]. This will result in lack of information downstream. We recommend to change your metadata or this contrast_tbl appropriately.")
  }
  
  if(length(groups_oi_contrast_tbl) != length(contrast_tbl$group)){
    warning("According to your contrast_tbl, some of your contrasts will be assigned to the same group. This should not be a problem if this was intended, but be aware not to make mistakes in the further interpretation and plotting of the results.")
  }
  
  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce_receiver")
    }
  }
  
  # Check concordance ligand-receptor network and ligand-target network - and concordance with sce object features
  if (!is.matrix(ligand_target_matrix)){
    stop("ligand_target_matrix should be a matrix")
  }
  if (!is.data.frame(lr_network)){
    stop("lr_network should be a data frame / tibble")
  }

  if(! "ligand" %in% colnames(lr_network)){
    if("from" %in% colnames(lr_network)){
      lr_network = lr_network %>% dplyr::rename(ligand = from)
    } else {
      stop("The ligand-receptor network should have the columns ligand and receptor (or: from and to)")
    }
  }
  if(! "receptor" %in% colnames(lr_network)){
    if("to" %in% colnames(lr_network)){
      lr_network = lr_network %>% dplyr::rename(receptor = to)
    } else {
      stop("The ligand-receptor network should have the columns ligand and receptor (or: from and to)")
    }
  }
  
  ligands_lrnetwork = lr_network$ligand %>% unique()
  receptors_lrnetwork = lr_network$receptor %>% unique()
  ligands_ligand_target_matrix = colnames(ligand_target_matrix)

  if(length(ligands_ligand_target_matrix) < length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }
  if(length(ligands_lrnetwork) < length(ligands_ligand_target_matrix)){
    warning("Not all Ligands from your ligand-target matrix are in the ligand-receptor network")
  }
  

  if(length(rownames(sce) %>% generics::intersect(ligands_lrnetwork)) < 25 ){
    warning("Less than 25 ligands from your ligand-receptor network are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(rownames(sce) %>% generics::intersect(ligands_ligand_target_matrix)) < 25 ){
    warning("Less than 25 ligands from your ligand-target matrix are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(rownames(sce) %>% generics::intersect(receptors_lrnetwork)) < 25 ){
    warning("Less than 25 receptors from your ligand-receptor network are in your expression matrix of the receiver cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }

  if(length(prioritizing_weights) != 8 | !is.double(prioritizing_weights)) {
    stop("prioritizing_weights should be a numeric vector with length 8")
  }
  names_prioritizing_weights = c("de_ligand",
                                 "de_receptor",
                                 "activity_scaled",
                                 "exprs_ligand",
                                 "exprs_receptor",
                                 "frac_exprs_ligand_receptor",
                                 "abund_sender",
                                 "abund_receiver")
  if(sum(names_prioritizing_weights %in% names(prioritizing_weights)) != length(names_prioritizing_weights)) {
    stop("prioritizing_weights should be have the correct names. Check the vignettes and code documentation")
  }


  if(!is.character(assay_oi_pb)){
    stop("assay_oi_pb should be a character vector")
  } else {
    if(assay_oi_pb != "counts"){
      warning("are you sure you don't want to use the counts assay?")
    }
  }
  if(!is.character(fun_oi_pb)){
    stop("fun_oi_pb should be a character vector")
  }
  if(!is.character(de_method_oi)){
    stop("de_method_oi should be a character vector")
  }

  if(!is.double(min_cells)){
    stop("min_cells should be numeric")
  } else {
    if(min_cells <= 0) {
      warning("min_cells is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(!is.double(logFC_threshold)){
    stop("logFC_threshold should be numeric")
  } else {
    if(logFC_threshold <= 0) {
      warning("logFC_threshold is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(!is.double(p_val_threshold)){
    stop("p_val_threshold should be numeric")
  } else {
    if(p_val_threshold <= 0 | p_val_threshold > 1) {
      warning("p_val_threshold is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.10, 0 excluded.")
    }
  }
  if(!is.double(fraction_cutoff)){
    stop("fraction_cutoff should be numeric")
  } else {
    if(fraction_cutoff <= 0 | fraction_cutoff > 1) {
      stop("fraction_cutoff is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.25, 0 excluded.")
    }
  }
  if(!is.double(top_n_target)){
    stop("top_n_target should be numeric")
  } else {
    if(top_n_target <= 0 ) {
      warning("top_n_target is now 0 or smaller. We recommend having a positive, non-zero value for this parameter.")
    }
  }
  if(!is.logical(p_val_adj)){
    stop("p_val_adj should be TRUE or FALSE")
  }
  if(!is.logical(verbose)){
    stop("verbose should be TRUE or FALSE")
  }
  if(!is.logical(empirical_pval)){
    stop("empirical_pval should be TRUE or FALSE")
  }
  if(!is.logical(return_lr_prod_matrix)){
    stop("return_lr_prod_matrix should be TRUE or FALSE")
  }
  

  
  
  if(!is.double(n.cores)){
    stop("n.cores should be numeric")
  } else {
    if(n.cores <= 0 ) {
      warning("n.cores is now 0 or smaller. We recommend having a positive, non-zero value for this parameter.")
    }
  }
  if(is.null(senders_oi)){
    senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
  }
  if(is.null(receivers_oi)){
    receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
  }
  sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  
  ### Perform the DE analysis ----------------------------------------------------------------

  if(verbose == TRUE){
    print("Calculate differential expression for all cell types")
  }
  
  if(findMarkers == FALSE){
    DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells,
                          assay_oi_pb = assay_oi_pb,
                          fun_oi_pb = fun_oi_pb,
                          de_method_oi = de_method_oi,
                          findMarkers = findMarkers, 
                          filterByExpr.min.count = filterByExpr.min.count, 
                          filterByExpr.min.total.count = filterByExpr.min.total.count, 
                          filterByExpr.large.n = filterByExpr.large.n, 
                          filterByExpr.min.prop = filterByExpr.min.prop)
  } else {
    DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells,
                          assay_oi_pb = assay_oi_pb,
                          fun_oi_pb = fun_oi_pb,
                          de_method_oi = de_method_oi,
                          findMarkers = findMarkers,
                          contrast_tbl = contrast_tbl)
  }

  if(empirical_pval == TRUE){
    if(findMarkers == TRUE){
      DE_info_emp = get_empirical_pvals(DE_info$celltype_de_findmarkers)
    } else {
      DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
    }
  } 

  if(empirical_pval == FALSE){
    if(findMarkers == TRUE){
      celltype_de = DE_info$celltype_de_findmarkers
    } else {
      celltype_de = DE_info$celltype_de$de_output_tidy
    }
  } else {
    celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
  }
  
  senders_oi = celltype_de$cluster_id %>% unique()
  receivers_oi = celltype_de$cluster_id %>% unique()
  genes_oi = celltype_de$gene %>% unique()
  sce = sce[genes_oi, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  
  sender_receiver_de = suppressMessages(combine_sender_receiver_de(
    sender_de = celltype_de,
    receiver_de = celltype_de,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network
  ))

  ### Receiver abundance plots + Calculate expression information
  if(verbose == TRUE){
    print("Make diagnostic abundance plots + Calculate expression information")
  }
  
  abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)
  rm(sce)
  ### Use the DE analysis for defining the DE genes in the receiver cell type and perform NicheNet ligand activity and ligand-target inference ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Calculate NicheNet ligand activities and ligand-target links")
  }
  ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = receivers_oi,
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )))

  ### Combine the three types of information calculated above to prioritize ligand-receptor interactions ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Combine all the information in prioritization tables")
  }
  ### Remove types of information that we don't need anymore:

  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  
  if(!is.na(batches)){
    grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group",batches)
  } else {
    grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group")
  }

  prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    prioritizing_weights = prioritizing_weights,
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender
  ))

  # Prepare Unsupervised analysis of samples! ------------------------------------------------------------------------------------------------------------
  
  if(return_lr_prod_matrix == TRUE){
    
    ids_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% dplyr::pull(id) %>% unique()
    
    lr_prod_df = abundance_expression_info$sender_receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sample, id, ligand_receptor_pb_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(id, ligand_receptor_pb_prod)
    lr_prod_mat = lr_prod_df %>% dplyr::select(-sample) %>% data.frame() %>% as.matrix()
    rownames(lr_prod_mat) = lr_prod_df$sample
    
    col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
    row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
    
    lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(col_remove),colnames(.) %>% generics::setdiff(col_remove)]
  } else {
    lr_prod_mat = NULL
  }

  # Add information on prior knowledge and expression correlation between LR and target expression ------------------------------------------------------------------------------------------------------------
  lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj, top_n_LR = top_n_LR)
  
  multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver, 
    abundance_data_sender = abundance_expression_info$abundance_data_sender, 
    celltype_de = celltype_de,
    # sender_receiver_info = abundance_expression_info$sender_receiver_info,
    # sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    lr_prod_mat = lr_prod_mat,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
  
  multinichenet_output = multinichenet_output %>% make_lite_output(top_n_LR = top_n_LR)
  
  return(multinichenet_output)

}

