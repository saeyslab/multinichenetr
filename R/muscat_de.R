#' @title perform_muscat_de_analysis
#'
#' @description \code{perform_muscat_de_analysis} Perform differential expression analysis via Muscat - Pseudobulking approach.
#' @usage perform_muscat_de_analysis(seurat_obj, sample_id, celltype_id, group_id, covariates, contrasts, assay_oi_sce = "RNA", assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @param contrasts String indicating the contrasts of interest (= which groups/conditions will be compared) for the differential expression and MultiNicheNet analysis. 
#' We will demonstrate here a few examples to indicate how to write this. Check the limma package manuals for more information about defining design matrices and contrasts for differential expression analysis.
#' If wanting to compare group A vs B: `contrasts_oi = c("'A-B'")`
#' If wanting to compare group A vs B & B vs A: `contrasts_oi = c("'A-B','B-A'")`
#' If wanting to compare group A vs B & A vs C & A vs D: `contrasts_oi = c("'A-B','A-C', 'A-D'")`
#' If wanting to compare group A vs B and C: `contrasts_oi = c("'A-(B+C)/2'")`
#' If wanting to compare group A vs B, C and D: `contrasts_oi = c("'A-(B+C+D)/3'")`
#' If wanting to compare group A vs B, C and D & B vs A,C,D: `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")`
#' Note that the groups A, B, ... should be present in the meta data column 'group_id'.
#' @return List with a SingleCellExperiment object and the output of the differential expression analysis (`muscat::pbDS()`)
#'
#' @import Seurat
#' @import dplyr
#' @import muscat
#' @importFrom SummarizedExperiment assayNames 
#' @importFrom S4Vectors metadata
#' @importFrom limma makeContrasts
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' senders_oi = Idents(seurat_obj) %>% unique()
#' receivers_oi = Idents(seurat_obj) %>% unique()
#' celltype_de = perform_muscat_de_analysis(
#'    seurat_obj = seurat_obj,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#'}
#'
#' @export
#'
perform_muscat_de_analysis = function(seurat_obj, sample_id, celltype_id, group_id, covariates, contrasts, assay_oi_sce = "RNA", assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10){

  if (class(seurat_obj) != "Seurat") {
    stop("seurat_obj should be a Seurat object")
  }
  if (!celltype_id %in% colnames(seurat_obj@meta.data)) {
    stop("celltype_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if (celltype_id != make.names(celltype_id)) {
    stop("celltype_id should be a syntactically valid R name - check make.names")
  }
  if (!sample_id %in% colnames(seurat_obj@meta.data)) {
    stop("sample_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if (sample_id != make.names(sample_id)) {
    stop("sample_id should be a syntactically valid R name - check make.names")
  }
  if (!group_id %in% colnames(seurat_obj@meta.data)) {
    stop("group_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if (group_id != make.names(group_id)) {
    stop("group_id should be a syntactically valid R name - check make.names")
  }
  
  if(is.double(seurat_obj@meta.data[,celltype_id])){
    stop("seurat_obj@meta.data[,celltype_id] should be a character vector or a factor")
  }
  if(is.double(seurat_obj@meta.data[,group_id])){
    stop("seurat_obj@meta.data[,group_id] should be a character vector or a factor")
  }
  if(is.double(seurat_obj@meta.data[,sample_id])){
    stop("seurat_obj@meta.data[,sample_id] should be a character vector or a factor")
  }
  
  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(seurat_obj@meta.data[,celltype_id])){
    if(levels(seurat_obj@meta.data[,celltype_id]) != make.names(levels(seurat_obj@meta.data[,celltype_id])))
      stop("The levels of the factor seurat_obj@meta.data[,celltype_id] should be a syntactically valid R names - see make.names")
  }
  if(is.factor(seurat_obj@meta.data[,group_id])){
    if(levels(seurat_obj@meta.data[,group_id]) != make.names(levels(seurat_obj@meta.data[,group_id])))
      stop("The levels of the factor seurat_obj@meta.data[,group_id] should be a syntactically valid R names - see make.names")
  }
  if(is.factor(seurat_obj@meta.data[,sample_id])){
    if(levels(seurat_obj@meta.data[,sample_id]) != make.names(levels(seurat_obj@meta.data[,sample_id])))
      stop("The levels of the factor seurat_obj@meta.data[,sample_id] should be a syntactically valid R names - see make.names")
  }
  
  if(!is.character(contrasts)){
    stop("contrasts should be a character vector")
  }

  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = seurat_obj@meta.data[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts, "'") %>% unlist() %>% unique() %>%
    stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()
  
  if(length(contrasts) != 1 | !is.character(contrasts)){
    stop("contrasts should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }
  
  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_simplified = stringr::str_split(contrasts, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()
  
  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }

  if(!is.na(covariates)){
    if (sum(covariates %in% colnames(seurat_obj@meta.data)) != length(covariates) ) {
      stop("covariates should be NA or all present as column name(s) in the metadata dataframe of seurat_obj_receiver")
    }
  }
  
  if(!is.character(assay_oi_sce)){
    stop("assay_oi_sce should be a character vector")
  } else {
    if(assay_oi_sce != "RNA"){
      warning("are you sure you don't want to use the RNA assay?")
    }
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

  requireNamespace("Seurat")
  requireNamespace("dplyr")
  
  # convert seurat to SCE object
  sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = assay_oi_sce)

  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  # drop all other colData columns ----------------- change to false

  pb = muscat::aggregateData(sce,
                             assay = assay_oi_pb, fun = fun_oi_pb,
                             by = c("cluster_id", "sample_id"))

  # prepare the design and contrast matrix for the muscat DE analysis
  if(length(covariates) > 1){
    covariates_present = TRUE
  } else {
    if(!is.na(covariates)){
      covariates_present = TRUE
    } else {
      covariates_present = FALSE

    }
  }

  if(covariates_present){
    extra_metadata = seurat_obj@meta.data %>% dplyr::select(all_of(sample_id), all_of(covariates)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
  } else {
    extra_metadata = seurat_obj@meta.data %>% dplyr::select(all_of(sample_id)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
  }
  if('sample_id' != sample_id){
    extra_metadata$sample_id = extra_metadata[[sample_id]]
  }
  ei = S4Vectors::metadata(sce)$experiment_info

  ei = ei %>%  dplyr::inner_join(extra_metadata, by = "sample_id")

  if(covariates_present){
    covariates_string = paste0("ei$",covariates) %>% paste(collapse = " + ")
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id + ", covariates_string, " ) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id), make.names(colnames(design)[(length(levels(ei$group_id))+1):length(colnames(design))])) )

  } else {
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id)))

  }

  contrast = eval(parse(text=paste("limma::makeContrasts(", contrasts, ",levels=design)",sep="")))

  # check which cell types will be excluded - Inner Muscat code
  n_cells = S4Vectors::metadata(pb)$n_cells
  celltypes = SummarizedExperiment::assayNames(pb)
  names(celltypes) = celltypes
  excluded_celltypes = celltypes %>% lapply(function(k){
    rmv = n_cells[k, ] < min_cells
    d = design[colnames(y <- pb[, !rmv]), , drop = FALSE]
    if (any(tabulate(y$group_id) < 2) || qr(d)$rank == nrow(d) ||
        qr(d)$rank < ncol(d)) {
      return(k)
    }
  }) %>% unlist() %>% unique()

  if(length(excluded_celltypes) > 0){
    print("excluded cell types are:")
    print(excluded_celltypes)
    print("These celltypes are not considered in the analysis. After removing samples that contain less cells than the required minimal, some groups don't have 2 or more samples anymore. As a result the analysis cannot be run. To solve this: decrease the number of min_cells or change your group_id and pool all samples that belong to groups that are not of interest! ")
  }

  # run DS analysis
  res = muscat::pbDS(pb, method = de_method_oi , design = design, contrast = contrast, min_cells = min_cells, verbose = FALSE, filter = "none")

  return(list(sce = sce, de_output = res))
}
