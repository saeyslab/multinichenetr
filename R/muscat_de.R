#' @title perform_muscat_de_analysis
#'
#' @description \code{perform_muscat_de_analysis}  XXXX
#' @usage perform_muscat_de_analysis(seurat_obj, sample_id, celltype_id, group_id, covariates, contrasts, assay_oi_sce = "RNA", assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#'
#' @return XXXX
#'
#' @import Seurat
#' @import dplyr
#' @import muscat
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
perform_muscat_de_analysis = function(seurat_obj, sample_id, celltype_id, group_id, covariates, contrasts, assay_oi_sce = "RNA", assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10){

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
  ei = metadata(sce)$experiment_info

  ei = ei %>%  dplyr::inner_join(extra_metadata, by = "sample_id")

  if(covariates_present){
    covariates_string = paste0("ei$",covariates) %>% paste(collapse = " + ")
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id + ", covariates_string, " ) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id), make.names(colnames(design)[(length(levels(ei$group_id))+1):length(colnames(design))])) )

  } else {
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id)))

  }

  contrast = eval(parse(text=paste("makeContrasts(", contrasts, ",levels=design)",sep="")))

  # check which cell types will be excluded
  n_cells = metadata(pb)$n_cells
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
