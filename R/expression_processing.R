#' @title get_muscat_exprs_frac
#'
#' @description \code{get_muscat_exprs_frac} Calculate sample- and group-average of fraction of cells in a cell type expressing a gene.
#' @usage get_muscat_exprs_frac(sce, sample_id, celltype_id, group_id)
#'
#' @inheritParams multi_nichenet_analysis
#'
#' @return List with two dataframes: one with fraction of cells in a cell type expressing a gene, averaged per sample; and one averaged per group.
#'
#' @import dplyr
#' @import muscat
#' @import tibble
#' @import tidyr
#' @importFrom purrr map
#' @importFrom magrittr set_names
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' muscat_exprs_frac = get_muscat_exprs_frac(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_muscat_exprs_frac = function(sce, sample_id, celltype_id, group_id){

  requireNamespace("dplyr")
  
  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  #

  samples = sce$sample_id %>% unique() %>% as.character()
  groups = sce$group_id  %>% unique() %>% as.character()

  frq = muscat::calcExprFreqs(sce, assay = "counts", th = 0) # gives NaN sometimes...
  if(is.null(rownames(frq))){
    all_genes = rownames(sce) # eg in case the first cell type in frq is absent in one sample -- then the rownames of frq will be NULL
  } else {
    all_genes = rownames(frq)
  }  
  celltypes = sce$cluster_id %>% unique() %>% as.character()
  
  frq_lists = celltypes %>% lapply(function(celltype_oi, frq, all_genes){
    frq_celltype = frq@assays@data[[celltype_oi]]
    
    celltype_genes = rownames(frq_celltype)
    if(sum(celltype_genes == all_genes) != length(all_genes) ){
      warning(paste0("For Celltype ",celltype_oi, " gene names got lost while calculating the fraction of expression of each gene with muscat::calcExprFreqs (due to a bug in this function). We temporality fixed this ourselves for the moment."))
      rownames(frq_celltype) = all_genes
    }
    
    frq_celltype_samples = frq_celltype[,samples]
    frq_celltype_samples = frq_celltype_samples %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, fraction_sample, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)
    
    frq_celltype_groups = frq_celltype[,groups]
    frq_celltype_groups = frq_celltype_groups %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(group, fraction_group, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)
    
    return(list(frq_celltype_samples = frq_celltype_samples, frq_celltype_groups = frq_celltype_groups))
  },frq, all_genes) %>% magrittr::set_names(sce$cluster_id %>% unique())
  
  
  frq_celltype_samples = frq_lists %>% purrr::map("frq_celltype_samples") %>% dplyr::bind_rows()
  frq_celltype_groups = frq_lists %>% purrr::map("frq_celltype_groups") %>% dplyr::bind_rows()
  rm(frq_lists)

  return(list(frq_celltype_samples = frq_celltype_samples, frq_celltype_groups = frq_celltype_groups))
}
#' @title get_muscat_exprs_avg
#'
#' @description \code{get_muscat_exprs_avg}  Calculate sample- and group-average of gene expression per cell type.
#' @usage get_muscat_exprs_avg(sce, sample_id, celltype_id, group_id)
#'
#' @inheritParams multi_nichenet_analysis
#'
#' @return Data frame with average gene expression per sample and per group.
#'
#' @import dplyr
#' @import muscat
#' @import tibble
#' @import tidyr
#' @importFrom scater logNormCounts computeLibraryFactors
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' muscat_exprs_avg = get_muscat_exprs_avg(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_muscat_exprs_avg = function(sce, sample_id, celltype_id, group_id){

  requireNamespace("dplyr")
  
  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  #
  if(! "logcounts" %in% SummarizedExperiment::assayNames(sce)){
    sce = scater::computeLibraryFactors(sce)
    sce = scater::logNormCounts(sce)
    # if wanting to use scran: --- uncomment this --> add this to the description file again (via add_data.R + add this to this import file)
    # set.seed(1919) # seed for reproducibility of the scran::quickCluster function (see https://bioconductor.org/books/release/OSCA/normalization.html)
    # 
    # clusters = suppressWarnings(scran::quickCluster(sce)) # why suppressWarnings --> on the lite dataset input ":regularize.values(x, y, ties, missing(ties), na.rm = na.rm)" pops always up, not on the entire dataset. Avoid these types of warnings for checks + anyway: logcounts assay are not much used anyway in the prioritizaton and visualization.
    # sce = suppressWarnings(scran::computeSumFactors(sce, clusters=clusters))
    # sce = suppressWarnings(scater::logNormCounts(sce))
  }

  avg = muscat::aggregateData(sce, assay = "logcounts", fun = "mean", by = c("cluster_id", "sample_id"))

  avg_df = sce$cluster_id %>% unique() %>% lapply(function(celltype_oi, avg){
    avg_celltype = avg@assays@data[[celltype_oi]]

    avg_celltype_samples = avg_celltype[,sce$sample_id %>% unique()]
    avg_celltype_samples = avg_celltype_samples %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, average_sample, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)

  },avg) %>% dplyr::bind_rows()

  return(avg_df)
}
#' @title get_pseudobulk_logCPM_exprs
#'
#' @description \code{get_pseudobulk_logCPM_exprs}  Calculate the 'library-size' normalized pseudbulk counts per sample for each gene - returned values are similar to logCPM. 
#' @usage get_pseudobulk_logCPM_exprs(sce, sample_id, celltype_id, group_id, batches = NA, assay_oi_pb = "counts", fun_oi_pb = "sum")
#'
#' @inheritParams multi_nichenet_analysis
#'
#' @return Data frame with logCPM-like values of the library-size corrected pseudobulked counts (`pb_sample`) per gene per sample. pb_sample = log2( ((pb_raw/effective_library_size) \* 1000000) + 1). effective_library_size  = lib.size \* norm.factors (through edgeR::calcNormFactors).
#'
#' @import dplyr
#' @import muscat
#' @import tibble
#' @import tidyr
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom sva ComBat_seq
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' pseudobulk_logCPM_exprs = get_pseudobulk_logCPM_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_pseudobulk_logCPM_exprs = function(sce, sample_id, celltype_id, group_id, batches = NA, assay_oi_pb = "counts", fun_oi_pb = "sum"){
  
  requireNamespace("dplyr")
  
  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  #
  
  pb = muscat::aggregateData(sce, assay = assay_oi_pb, fun = fun_oi_pb, by = c("cluster_id", "sample_id"))
  
  # Prepare the design (group and batch effects) for the combat correction
  if(length(batches) > 1){
    batches_present = TRUE
    batches = batches[1] ## only keep the first batch for batch correction!
    warning("You used more than 1 batch/batch to correct for. This is OK for the DE analysis. But, for the Combat batch correction of pseudobulked counts (used in visualization), only the first batch will be considered as batch. If you want to use both as batch, redo the analysis after merging both in one variable.")
  } else {
    if(!is.na(batches)){
      batches_present = TRUE
    } else {
      batches_present = FALSE
    }
  }
  if(batches_present){ ## then perform the correction!
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batches)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batches")
    ei = S4Vectors::metadata(sce)$experiment_info
    ei = ei %>%  dplyr::inner_join(extra_metadata, by = "sample_id")
    
    # do a check: will we able to correct for the batches on a group-by-group basis?
    n_combinations_observed = ei %>% dplyr::select(group_id, batches) %>% dplyr::distinct() %>% nrow()
    n_combinations_possible = length(levels(ei$group_id)) * length(levels(ei$batches))
    if(n_combinations_observed != n_combinations_possible){
      warning("Not all possible group-batch/batch combinations are present in your data. This will result in errors during the batch effect correction process of Combat and/or Muscat DE analysis. Please reconsider the groups and batches you defined.")
    } 
    
    # do another necessary check: each batch-group combinations should have at least one or preferably two observations
    # print(ei %>% dplyr::select(sample_id, group_id, batches) %>% dplyr::distinct() %>% dplyr::group_by(group_id, batches) %>% dplyr::count() %>% arrange(n))
    # print(ei %>% dplyr::select(sample_id, batches) %>% dplyr::distinct() %>% dplyr::group_by(batches) %>% dplyr::count() %>% arrange(n))
    batch2 = as.factor(ei$batches)
    if (any(table(batch2) <= 1)) {
      warning("For some batch values, only one sample is present in your data. Combat correction of the expression for downstream visualization cannot handle this.")
      pb_df = sce$cluster_id %>% unique() %>% lapply(function(celltype_oi, pb){ # no correction of the pseudobulk counts
        
        pseudobulk_counts_celltype = edgeR::DGEList(pb@assays@data[[celltype_oi]])
        
        non_zero_samples = pseudobulk_counts_celltype %>% apply(2,sum) %>% .[. > 0] %>% names()
        pseudobulk_counts_celltype = pseudobulk_counts_celltype[,non_zero_samples]
        
        pseudobulk_counts_celltype = edgeR::calcNormFactors(pseudobulk_counts_celltype)
        
        pseudobulk_counts_celltype_df = dplyr::inner_join(
          pseudobulk_counts_celltype$sample %>% data.frame() %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(effective_library_size  = lib.size * norm.factors),
          pseudobulk_counts_celltype$counts %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, pb_raw, -gene)
        )
        
        pseudobulk_counts_celltype_df = pseudobulk_counts_celltype_df %>% dplyr::mutate(pb_norm = pb_raw / effective_library_size) %>% dplyr::mutate(pb_sample = log2( (pb_norm * 1000000) + 1)) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi) # +1: pseudocount - should be introduced after library size correction, otherwise: we think some samples have a higher expression even though it was zero!
        
      },pb) %>% dplyr::bind_rows() %>% dplyr::select(gene, sample, pb_sample, celltype) %>% dplyr::distinct()
      
    } else {
      pb_df = sce$cluster_id %>% unique() %>% lapply(function(celltype_oi, pb, ei){
        # get the count matrix
        count_matrix = pb@assays@data[[celltype_oi]] %>% .[,ei$sample_id]
        non_zero_samples = count_matrix %>% apply(2,sum) %>% .[. > 0] %>% names()
        count_matrix = count_matrix[,non_zero_samples]
        # adjust the count matrix
        ei = ei %>% dplyr::filter(sample_id %in% non_zero_samples)
        batch2 = as.factor(ei$batches) # still check whether we have enough samples for each celltype to correct for!!
        if (any(table(batch2) <= 1)) {
          warning(paste0("For some batch values for celltype ",celltype_oi, " only one sample is present in your data. Combat correction of the expression for downstream visualization cannot handle this. Therefore we continue with non-corrected pseudobulk expression values here."))
          adj_count_matrix = count_matrix 
        } else {
          quiet <- function(x) { # to suppress sva-combat output
            sink(tempfile()) 
            on.exit(sink()) 
            invisible(force(x)) 
          } 
          adj_count_matrix = quiet(sva::ComBat_seq(count_matrix, batch=ei$batches, group=ei$group_id, full_mod=TRUE)) # to suppress its output
        }
        # normalize the adjusted count matrix, just like we do for the non-adjusted one
        pseudobulk_counts_celltype = edgeR::DGEList(adj_count_matrix)
        pseudobulk_counts_celltype = edgeR::calcNormFactors(pseudobulk_counts_celltype)
        
        pseudobulk_counts_celltype_df = dplyr::inner_join(
          pseudobulk_counts_celltype$sample %>% data.frame() %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(effective_library_size  = lib.size * norm.factors),
          pseudobulk_counts_celltype$counts %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, pb_raw, -gene)
        )
        
        pseudobulk_counts_celltype_df = pseudobulk_counts_celltype_df %>% dplyr::mutate(pb_norm = pb_raw / effective_library_size) %>% dplyr::mutate(pb_sample = log2( (pb_norm * 1000000) + 1)) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi) # +1: pseudocount - should be introduced after library size correction, otherwise: we think some samples have a higher expression even though it was zero!
        
      },pb, ei) %>% dplyr::bind_rows() %>% dplyr::select(gene, sample, pb_sample, celltype) %>% dplyr::distinct()
    }
    
  } else { # no correction of the pseudobulk counts
    pb_df = sce$cluster_id %>% unique() %>% lapply(function(celltype_oi, pb){
      
      pseudobulk_counts_celltype = edgeR::DGEList(pb@assays@data[[celltype_oi]])
      
      non_zero_samples = pseudobulk_counts_celltype %>% apply(2,sum) %>% .[. > 0] %>% names()
      pseudobulk_counts_celltype = pseudobulk_counts_celltype[,non_zero_samples]
      
      pseudobulk_counts_celltype = edgeR::calcNormFactors(pseudobulk_counts_celltype)
      
      pseudobulk_counts_celltype_df = dplyr::inner_join(
        pseudobulk_counts_celltype$sample %>% data.frame() %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(effective_library_size  = lib.size * norm.factors),
        pseudobulk_counts_celltype$counts %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, pb_raw, -gene)
      )
      
      pseudobulk_counts_celltype_df = pseudobulk_counts_celltype_df %>% dplyr::mutate(pb_norm = pb_raw / effective_library_size) %>% dplyr::mutate(pb_sample = log2( (pb_norm * 1000000) + 1)) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi) # +1: pseudocount - should be introduced after library size correction, otherwise: we think some samples have a higher expression even though it was zero!
      
    },pb) %>% dplyr::bind_rows() %>% dplyr::select(gene, sample, pb_sample, celltype) %>% dplyr::distinct()
  }
  
  
  return(pb_df)
}
#' @title fix_frq_df
#'
#' @description \code{fix_frq_df}  Fix muscat-feature/bug in fraction calculation: in case a there are no cells of a cell type in a sample, that expression fraction will be NA / NaN. Change these NA/NaN to 0.
#' @usage fix_frq_df(sce, frq_celltype_samples)
#'
#' @inheritParams multi_nichenet_analysis
#' @param frq_celltype_samples Sample-average data frame output of `get_muscat_exprs_frac`
#'
#' @return Fixed data frame with fraction of cells expressing a gene.
#'
#' @import dplyr
#' @importFrom magrittr set_names
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' frq_df = get_muscat_exprs_frac(sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id) %>% .$frq_celltype_samples
#' if(nrow(frq_df %>% dplyr::filter(is.na(fraction_sample))) > 0 | nrow(frq_df %>% dplyr::filter(is.nan(fraction_sample))) > 0) {
#'   frq_df = fix_frq_df(sce, frq_df)
#'   }
#' }
#'
#' @export
#'
fix_frq_df = function(sce, frq_celltype_samples){
  
  requireNamespace("dplyr")
  
  genes = rownames(sce)
  gene_mapping = genes %>% magrittr::set_names(seq(length(genes)))

  frq_celltype_samples_OK = frq_celltype_samples %>% dplyr::filter(gene %in% genes)

  # Fix 0's incase gene names absent
  frq_celltype_samples_FIX = frq_celltype_samples %>% dplyr::filter(!gene %in% genes)

  frq_celltype_samples_FIX = frq_celltype_samples_FIX %>% dplyr::mutate(gene = gene_mapping[gene])

  frq_celltype_samples_FIX = frq_celltype_samples_FIX %>% dplyr::mutate(fraction_sample = 0)

  frq_celltype_samples = frq_celltype_samples_OK %>% dplyr::bind_rows(frq_celltype_samples_FIX)

  # Fix NA/NaNs to Os in case gene names present
  frq_celltype_samples_OK = frq_celltype_samples %>% dplyr::filter(!is.na(fraction_sample)) 
  
  frq_celltype_samples_FIX = frq_celltype_samples %>% dplyr::filter(is.na(fraction_sample)) 
  frq_celltype_samples_FIX = frq_celltype_samples_FIX %>% dplyr::mutate(fraction_sample = 0)
  
  frq_celltype_samples = frq_celltype_samples_OK %>% dplyr::bind_rows(frq_celltype_samples_FIX)
  
  return(frq_celltype_samples)

}
#' @title get_avg_pb_exprs
#'
#' @description \code{get_avg_pb_exprs}  Calculate the average and normalized pseudobulk expression of each gene per sample and per group.
#' @usage get_avg_pb_exprs(sce, sample_id, celltype_id, group_id, batches = NA, min_cells = 10)
#'
#' @inheritParams multi_nichenet_analysis
#' 
#' @return List containing data frames with average and normalized pseudobulk of expression per sample and per group.
#'
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_pb_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_avg_pb_exprs = function(sce, sample_id, celltype_id, group_id, batches = NA, min_cells = 10){
  
  requireNamespace("dplyr")
  
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

  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  }
  ## calculate averages, fractions, relative abundance of a cell type in a group

  # calculate average expression
  avg_df = get_muscat_exprs_avg(sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)

  # calculate pseudobulked counts
  pb_df = get_pseudobulk_logCPM_exprs(sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id, batches = batches, assay_oi_pb = "counts", fun_oi_pb = "sum") # should be these parameters
  
  # check whether something needs to be fixed
  if(nrow(avg_df %>% dplyr::filter(is.na(average_sample))) > 0 | nrow(avg_df %>% dplyr::filter(is.nan(average_sample))) > 0) {
    warning("There are some genes with NA average expression.")
  }

  # calculate these above metrics per group
  metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  if ("sample_id" != sample_id) {
    metadata$sample_id = metadata[[sample_id]]
  }
  if ("group_id" != sample_id) {
    metadata$group_id = metadata[[group_id]]
  }
  if ("celltype_id" != celltype_id) {
    metadata$celltype_id = metadata[[celltype_id]]
  }
  #metadata_abundance = metadata %>% dplyr::select(sample_id, group_id, celltype_id) %>% tibble::as_tibble()
  #colnames(metadata_abundance) =c("sample", "group", "celltype")
  #abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample , celltype) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample , group ), by = "sample")
  #abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, levels = c(TRUE,FALSE)))
  
  metadata_abundance = metadata %>% dplyr::select(sample_id, group_id, celltype_id) %>% tibble::as_tibble()
  #colnames(metadata_abundance) =c("sample", "group", "celltype")
  #abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample , celltype) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample , group ), by = "sample")
  #abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, levels = c(TRUE,FALSE)))
  
  # Ensure column names are consistent
  colnames(metadata_abundance) <- c("sample", "group", "celltype")
  
  # Get unique sample and celltype combinations, fill in missing with n = 0
  all_combinations <- metadata_abundance %>%
    dplyr::distinct(sample, celltype) %>%
    tidyr::expand(sample, celltype)
  
  # Calculate abundance
  abundance_data <- metadata_abundance %>%
    dplyr::group_by(sample, celltype) %>%
    dplyr::count() %>%
    dplyr::right_join(all_combinations, by = c("sample", "celltype")) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n))  # Set missing counts to 0
  
  # Add group information
  abundance_data <- abundance_data %>%
    dplyr::left_join(metadata_abundance %>% distinct(sample, group), by = "sample")
  
  abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, levels = c(TRUE,FALSE)))
  
  grouping_df = metadata %>% dplyr::select(sample_id, group_id) %>% 
    tibble::as_tibble() %>% dplyr::distinct() %>% dplyr::rename(sample = sample_id, 
                                                                group = group_id)
  
  grouping_df_filtered = grouping_df %>% inner_join(abundance_data) %>% filter(keep_sample == TRUE) # continue only with samples that have sufficient cells of a cell type in a sample
  
  avg_df_group = avg_df %>% dplyr::inner_join(grouping_df_filtered) %>% 
    dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(average_group = mean(average_sample))
  pb_df_group = pb_df %>% dplyr::inner_join(grouping_df_filtered) %>% 
    dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(pb_group = mean(pb_sample))
  
  return(list(avg_df = avg_df, pb_df = pb_df, avg_df_group = avg_df_group, pb_df_group = pb_df_group))

}
#' @title get_frac_exprs
#'
#' @description \code{get_frac_exprs}  Calculate the average fraction of expression of each gene per sample and per group.
#' @usage get_frac_exprs(sce, sample_id, celltype_id, group_id, batches = NA, min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 0.5)
#'
#' @inheritParams multi_nichenet_analysis
#' 
#' @return List containing data frames with the fraction of expression per sample and per group.
#'
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' frac_info = get_frac_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_frac_exprs = function(sce, sample_id, celltype_id, group_id, batches = NA, min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 0.5){
  
  requireNamespace("dplyr")
  
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
  
  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  }
  # calculate fraction of expression
  frq_df = get_muscat_exprs_frac(sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id) %>% .$frq_celltype_samples
  
  if(nrow(frq_df %>% dplyr::filter(is.na(fraction_sample))) > 0 | nrow(frq_df %>% dplyr::filter(is.nan(fraction_sample))) > 0) {
    warning("There are some genes with NA/NaN fraction of expression. This is the result of the muscat function `calcExprFreqs` which will give NA/NaN when there are no cells of a particular cell type in a particular group or no cells of a cell type in one sample. As a temporary fix, we give all these genes an expression fraction of 0 in that group for that cell type")
    frq_df = fix_frq_df(sce, frq_df)
  }
  
  # calculate these above metrics per group
  metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  if ("sample_id" != sample_id) {
    metadata$sample_id = metadata[[sample_id]]
  }
  if ("group_id" != sample_id) {
    metadata$group_id = metadata[[group_id]]
  }
  if ("celltype_id" != celltype_id) {
    metadata$celltype_id = metadata[[celltype_id]]
  }
  metadata_abundance = metadata %>% dplyr::select(sample_id, group_id, celltype_id) %>% tibble::as_tibble()
  #colnames(metadata_abundance) =c("sample", "group", "celltype")
  #abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample , celltype) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample , group ), by = "sample")
  #abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, levels = c(TRUE,FALSE)))
  
  # Ensure column names are consistent
  colnames(metadata_abundance) <- c("sample", "group", "celltype")
  
  # Get unique sample and celltype combinations, fill in missing with n = 0
  all_combinations <- metadata_abundance %>%
    dplyr::distinct(sample, celltype) %>%
    tidyr::expand(sample, celltype)
  
  # Calculate abundance
  abundance_data <- metadata_abundance %>%
    dplyr::group_by(sample, celltype) %>%
    dplyr::count() %>%
    dplyr::right_join(all_combinations, by = c("sample", "celltype")) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n))  # Set missing counts to 0
  
  # Add group information
  abundance_data <- abundance_data %>%
    dplyr::left_join(metadata_abundance %>% distinct(sample, group), by = "sample")
  
  abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, levels = c(TRUE,FALSE)))
  
  
  grouping_df = metadata %>% dplyr::select(sample_id, group_id) %>% 
    tibble::as_tibble() %>% dplyr::distinct() %>% dplyr::rename(sample = sample_id, 
                                                                group = group_id)
  
  grouping_df_filtered = grouping_df %>% inner_join(abundance_data) %>% filter(keep_sample == TRUE) # continue only with samples that have sufficient cells of a cell type in a sample
  
  print(paste0("Samples are considered if they have more than ", min_cells, " cells of the cell type of interest")) 
  
  frq_df_group = frq_df %>% dplyr::inner_join(grouping_df_filtered) %>% 
    dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(fraction_group = mean(fraction_sample))
  
  # define whether genes are expressed - inspired by edgeR::filterByExprs but more suited for pseudobulk data  
  # expressed in a cell type = for n samples: non-zero counts in >= fraction_cutoff% of cells in a sample
  # n = min_sample_prop fraction of samples of the smallest group
  # min_sample_prop = 0.5 by default
  # fraction_cutoff = 0.05 by default
  n_smallest_group_tbl = grouping_df_filtered %>% dplyr::group_by(group, celltype) %>% dplyr::count() %>% dplyr::group_by(celltype) %>% dplyr::summarize(n_smallest_group = min(n)) %>% dplyr::mutate(n_min = min_sample_prop * n_smallest_group) %>% dplyr::mutate(n_min = pmax(n_min, 2)) %>% distinct()
  
  print(paste0("Genes with non-zero counts in at least ",fraction_cutoff*100, "% of cells of a cell type of interest in a particular sample will be considered as expressed in that sample.")) 
  
  for(i in seq(length(unique(n_smallest_group_tbl$celltype)))){
    celltype_oi = unique(n_smallest_group_tbl$celltype)[i]
    n_min = n_smallest_group_tbl %>% filter(celltype == celltype_oi) %>% pull(n_min)
    print(paste0("Genes expressed in at least ",n_min, " samples will considered as expressed in the cell type: ",celltype_oi)) 
  }
  
  frq_df = frq_df %>% dplyr::inner_join(grouping_df) %>% dplyr::mutate(expressed_sample = fraction_sample >= fraction_cutoff)
  
  expressed_df = frq_df %>% inner_join(n_smallest_group_tbl) %>% inner_join(abundance_data) %>% dplyr::group_by(gene, celltype) %>% dplyr::summarise(n_expressed = sum(expressed_sample)) %>% dplyr::mutate(expressed = n_expressed >= n_min) %>% distinct(celltype, gene, expressed)
  for(i in seq(length(unique(expressed_df$celltype)))){
    celltype_oi = unique(expressed_df$celltype)[i]
    n_genes = expressed_df %>% filter(celltype == celltype_oi) %>% filter(expressed == TRUE) %>% pull(gene) %>% unique() %>% length()
    print(paste0(n_genes, " genes are considered as expressed in the cell type: ",celltype_oi)) 
  }
  return(list(frq_df = frq_df, frq_df_group = frq_df_group, expressed_df = expressed_df))
  
}
#' @title get_frac_exprs_sampleAgnostic
#'
#' @description \code{get_frac_exprs_sampleAgnostic}  Calculate the average fraction of expression of each gene per group. All cells from all samples will be pooled per group/condition.
#' @usage get_frac_exprs_sampleAgnostic(sce, sample_id, celltype_id, group_id, batches = NA, min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 0.5)
#'
#' @inheritParams get_frac_exprs
#' @param min_sample_prop Default, and only recommended value = 1. Hereby, the gene should be expressed in at least one group/condition.
#' 
#' @return List containing data frames with the fraction of expression per sample and per group.
#'
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' frac_info = get_frac_exprs_sampleAgnostic(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_frac_exprs_sampleAgnostic = function(sce, sample_id, celltype_id, group_id, batches = NA, 
                                         min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 1){
  
  sample_id = group_id
  frq_df = get_muscat_exprs_frac(sce, sample_id = sample_id, 
                                 celltype_id = celltype_id, group_id = group_id) %>% .$frq_celltype_samples
  metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  metadata_abundance = metadata %>% dplyr::select(group_id, celltype_id) 
  print(head(metadata_abundance))
  metadata_abundance = cbind(metadata[, group_id], metadata_abundance) 
  print(head(metadata_abundance))
  colnames(metadata_abundance) = c("sample", "group", "celltype")
  print(head(metadata_abundance))
  metadata_abundance = metadata_abundance %>% tibble::as_tibble()
  print(head(metadata_abundance))
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% 
    dplyr::group_by(sample, celltype) %>% dplyr::count() %>% 
    dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% 
                        dplyr::distinct(sample, group), by = "sample")
  abundance_data = abundance_data %>% dplyr::mutate(keep_sample = n >= 
                                                      min_cells) %>% dplyr::mutate(keep_sample = factor(keep_sample, 
                                                                                                        levels = c(TRUE, FALSE)))
  grouping_df = abundance_data[,c("sample","group")] %>% 
    tibble::as_tibble() %>% dplyr::distinct() 
  grouping_df_filtered = grouping_df %>% inner_join(abundance_data) %>% 
    filter(keep_sample == TRUE)
  print(paste0("Groups are considered if they have more than ", 
               min_cells, " cells of the cell type of interest"))
  frq_df_group = frq_df %>% dplyr::inner_join(grouping_df_filtered) %>% 
    dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(fraction_group = mean(fraction_sample))
  n_smallest_group_tbl = grouping_df_filtered %>% dplyr::group_by(group, 
                                                                  celltype) %>% dplyr::count() %>% dplyr::group_by(celltype) %>% 
    dplyr::summarize(n_smallest_group = min(n)) %>% dplyr::mutate(n_min = min_sample_prop * 
                                                                    n_smallest_group) %>% distinct()
  print(paste0("Genes with non-zero counts in at least ", fraction_cutoff * 
                 100, "% of cells of a cell type of interest in a particular group/condition will be considered as expressed in that group/condition"))
  for (i in seq(length(unique(n_smallest_group_tbl$celltype)))) {
    celltype_oi = unique(n_smallest_group_tbl$celltype)[i]
    n_min = n_smallest_group_tbl %>% filter(celltype == celltype_oi) %>% 
      pull(n_min)
    print(paste0("Genes expressed in at least ", n_min, " group will considered as expressed in the cell type: ", 
                 celltype_oi))
  }
  frq_df = frq_df %>% dplyr::inner_join(grouping_df) %>% dplyr::mutate(expressed_sample = fraction_sample >= 
                                                                         fraction_cutoff)
  expressed_df = frq_df %>% inner_join(n_smallest_group_tbl) %>% 
    inner_join(abundance_data) %>% dplyr::group_by(gene, 
                                                   celltype) %>% dplyr::summarise(n_expressed = sum(expressed_sample)) %>% 
    dplyr::mutate(expressed = n_expressed >= n_min) %>% distinct(celltype, 
                                                                 gene, expressed)
  for (i in seq(length(unique(expressed_df$celltype)))) {
    celltype_oi = unique(expressed_df$celltype)[i]
    n_genes = expressed_df %>% filter(celltype == celltype_oi) %>% 
      filter(expressed == TRUE) %>% pull(gene) %>% unique() %>% 
      length()
    print(paste0(n_genes, " genes are considered as expressed in the cell type: ", 
                 celltype_oi))
  }
  return(list(frq_df = frq_df, frq_df_group = frq_df_group, 
              expressed_df = expressed_df))
}
#' @title process_info_to_ic
#'
#' @description \code{process_info_to_ic}  Process cell type expression information into intercellular communication focused information. Only keep information of ligands for the sender cell type setting, and information of receptors for the receiver cell type.
#' @usage process_info_to_ic(info_object, ic_type = "sender", lr_network)
#'
#' @inheritParams multi_nichenet_analysis
#' @param info_object Output of `get_avg_frac_exprs_abund`
#' @param ic_type "sender" or "receiver": indicates whether we should keep ligands or receptors respectively.
#'
#' @return List with expression information of ligands (sender case) or receptors (receiver case) - similar to output of `get_avg_frac_exprs_abund`.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' }
#'
#' @export
#'
process_info_to_ic = function(info_object, ic_type = "sender", lr_network){
  
  requireNamespace("dplyr")

  ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
  receptors = lr_network %>% dplyr::pull(receptor) %>% unique()

  if(ic_type == "sender"){
    avg_df = info_object$avg_df %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, avg_ligand = average_sample)
    frq_df = info_object$frq_df %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, fraction_ligand = fraction_sample)
    pb_df = info_object$pb_df %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, pb_ligand = pb_sample)
      
    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, avg_ligand_group = average_group)
    frq_df_group = info_object$frq_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, fraction_ligand_group = fraction_group)
    pb_df_group = info_object$pb_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, pb_ligand_group = pb_group)
    
    rel_abundance_df = info_object$rel_abundance_df %>% dplyr::rename(sender = celltype, rel_abundance_scaled_sender = rel_abundance_scaled)
  }
  if(ic_type == "receiver"){
    avg_df = info_object$avg_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, avg_receptor = average_sample)
    frq_df = info_object$frq_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor = fraction_sample)
    pb_df = info_object$pb_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, pb_receptor = pb_sample)
    
    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, avg_receptor_group = average_group)
    frq_df_group = info_object$frq_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor_group = fraction_group)
    pb_df_group = info_object$pb_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, pb_receptor_group = pb_group)
    
    rel_abundance_df = info_object$rel_abundance_df %>% dplyr::rename(receiver = celltype, rel_abundance_scaled_receiver = rel_abundance_scaled)
  }

  return(list(avg_df = avg_df, frq_df = frq_df, pb_df = pb_df,  avg_df_group = avg_df_group, frq_df_group = frq_df_group, pb_df_group = pb_df_group, rel_abundance_df = rel_abundance_df))
}


#' @title process_abund_info
#'
#' @description \code{process_abund_info}  Process cell type abundance information into intercellular communication focused information. 
#' @usage process_abund_info(abund_data, ic_type = "sender")
#'
#' @param abund_data Data frame with number of cells per cell type - sample combination
#' @param ic_type "sender" or "receiver": indicates whether we should rename celltype to sender or receiver.
#'
#' @return Data frame with number of cells per cell type - sample combination; but now adapted to sender/receiver nomenclature
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' min_cells = 10
#' metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)]
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' receiver_abundance_data = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' sender_abundance_data = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#' }
#'
#' @export
#'
process_abund_info = function(abund_data, ic_type = "sender"){
  
  requireNamespace("dplyr")
  
  abund_data = abund_data %>% rename(sample = sample_id, group = group_id)
  
  if(ic_type == "sender"){
    
    if("celltype_id_sender" %in% colnames(abund_data)){
      abund_data = abund_data %>% rename(sender = celltype_id_sender)
    } else {
      abund_data = abund_data %>% rename(sender = celltype_id)
    }
    
    abund_data = abund_data %>% rename(keep_sender = keep, n_cells_sender = n) %>% mutate(keep_sender = as.double(as.logical(keep_sender)))
  }
  if(ic_type == "receiver"){
    if("celltype_id_receiver" %in% colnames(abund_data)){
      abund_data = abund_data %>% rename(receiver = celltype_id_receiver)
    } else {
      abund_data = abund_data %>% rename(receiver = celltype_id)
    }
    
    abund_data = abund_data %>% rename(keep_receiver = keep, n_cells_receiver = n) %>% mutate(keep_receiver = as.double(as.logical(keep_receiver)))
  }
  
  return(abund_data)
}

#' @title combine_sender_receiver_info_ic
#'
#' @description \code{combine_sender_receiver_info_ic}  Link the ligand-expression information of the Sender cell type to the receptor-expression information of the Receiver cell type. Linking via prior knowledge ligand-receptor network.
#' @usage combine_sender_receiver_info_ic(sender_info, receiver_info, senders_oi, receivers_oi, lr_network)
#'
#' @inheritParams multi_nichenet_analysis
#' @param sender_info Output of `process_info_to_ic` with `ic_type = "sender"`
#' @param receiver_info Output of `process_info_to_ic` with `ic_type = "receiver"`
#' @param senders_oi Character vector indicating the names of the sender cell types of interest
#' @param receivers_oi Character vector indicating the names of the receiver cell types of interest
#'
#' @return List with data frames containing ligand-receptor sender-receiver combined expression information (see output of `get_avg_frac_exprs_abund` and `process_info_to_ic`)
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#'   
#' }
#'
#' @export
#'
combine_sender_receiver_info_ic = function(sender_info, receiver_info, senders_oi, receivers_oi, lr_network){

  requireNamespace("dplyr")
  
  # combine avg_df
  avg_df_sender = sender_info$avg_df %>% dplyr::filter(sender %in% senders_oi)
  avg_df_receiver = receiver_info$avg_df %>% dplyr::filter(receiver %in% receivers_oi)

  avg_df_sender_receiver = avg_df_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(avg_df_receiver, by = c("receptor","sample"))
  avg_df_sender_receiver = avg_df_sender_receiver %>% dplyr::mutate(ligand_receptor_prod = avg_ligand * avg_receptor) %>% dplyr::arrange(-ligand_receptor_prod) %>% dplyr::select(sample, sender, receiver, ligand, receptor, avg_ligand, avg_receptor, ligand_receptor_prod) %>% dplyr::distinct()

  # combine avg_df_group
  avg_df_group_sender = sender_info$avg_df_group %>% dplyr::filter(sender %in% senders_oi)
  avg_df_group_receiver = receiver_info$avg_df_group %>% dplyr::filter(receiver %in% receivers_oi)

  avg_df_group_sender_receiver = avg_df_group_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(avg_df_group_receiver, by = c("receptor","group"))
  avg_df_group_sender_receiver = avg_df_group_sender_receiver %>% dplyr::mutate(ligand_receptor_prod_group = avg_ligand_group * avg_receptor_group) %>% dplyr::arrange(-ligand_receptor_prod_group) %>% dplyr::select(group, sender, receiver, ligand, receptor, avg_ligand_group, avg_receptor_group, ligand_receptor_prod_group) %>% dplyr::distinct()

  # combine frq_df

  frq_df_sender = sender_info$frq_df %>% dplyr::filter(sender %in% senders_oi)
  frq_df_receiver = receiver_info$frq_df %>% dplyr::filter(receiver %in% receivers_oi)

  frq_df_sender_receiver = frq_df_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(frq_df_receiver, by = c("receptor","sample"))
  frq_df_sender_receiver = frq_df_sender_receiver %>% dplyr::mutate(ligand_receptor_fraction_prod = fraction_ligand * fraction_receptor) %>% dplyr::arrange(-ligand_receptor_fraction_prod) %>% dplyr::select(sample, sender, receiver, ligand, receptor, fraction_ligand, fraction_receptor, ligand_receptor_fraction_prod) %>% dplyr::distinct()

  # combine frq_df_group

  frq_df_group_sender = sender_info$frq_df_group %>% dplyr::filter(sender %in% senders_oi)
  frq_df_group_receiver = receiver_info$frq_df_group %>% dplyr::filter(receiver %in% receivers_oi)

  frq_df_group_sender_receiver = frq_df_group_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(frq_df_group_receiver, by = c("receptor","group"))
  frq_df_group_sender_receiver = frq_df_group_sender_receiver %>% dplyr::mutate(ligand_receptor_fraction_prod_group = fraction_ligand_group * fraction_receptor_group) %>% dplyr::arrange(-ligand_receptor_fraction_prod_group) %>% dplyr::select(group, sender, receiver, ligand, receptor, fraction_ligand_group, fraction_receptor_group, ligand_receptor_fraction_prod_group) %>% dplyr::distinct()

  # combine pb_df
  pb_df_sender = sender_info$pb_df %>% dplyr::filter(sender %in% senders_oi)
  pb_df_receiver = receiver_info$pb_df %>% dplyr::filter(receiver %in% receivers_oi)
  
  pb_df_sender_receiver = pb_df_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(pb_df_receiver, by = c("receptor","sample"))
  pb_df_sender_receiver = pb_df_sender_receiver %>% dplyr::mutate(ligand_receptor_pb_prod = pb_ligand * pb_receptor) %>% dplyr::arrange(-ligand_receptor_pb_prod) %>% dplyr::select(sample, sender, receiver, ligand, receptor, pb_ligand, pb_receptor, ligand_receptor_pb_prod) %>% dplyr::distinct()
  
  # combine pb_df_group
  pb_df_group_sender = sender_info$pb_df_group %>% dplyr::filter(sender %in% senders_oi)
  pb_df_group_receiver = receiver_info$pb_df_group %>% dplyr::filter(receiver %in% receivers_oi)
  
  pb_df_group_sender_receiver = pb_df_group_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(pb_df_group_receiver, by = c("receptor","group"))
  pb_df_group_sender_receiver = pb_df_group_sender_receiver %>% dplyr::mutate(ligand_receptor_pb_prod_group = pb_ligand_group * pb_receptor_group) %>% dplyr::arrange(-ligand_receptor_pb_prod_group) %>% dplyr::select(group, sender, receiver, ligand, receptor, pb_ligand_group, pb_receptor_group, ligand_receptor_pb_prod_group) %>% dplyr::distinct()
  
  # combine relative abundances
  rel_abundance_df_sender = sender_info$rel_abundance_df %>% dplyr::filter(sender %in% senders_oi)
  rel_abundance_df_receiver = receiver_info$rel_abundance_df %>% dplyr::filter(receiver %in% receivers_oi)

  rel_abundance_df_sender_receiver = rel_abundance_df_sender %>%  dplyr::inner_join(rel_abundance_df_receiver, by = "group") %>% dplyr::mutate(sender_receiver_rel_abundance_avg = 0.5*(rel_abundance_scaled_sender  + rel_abundance_scaled_receiver))

  # return
  return(list(avg_df = avg_df_sender_receiver, frq_df = frq_df_sender_receiver, pb_df = pb_df_sender_receiver,  avg_df_group = avg_df_group_sender_receiver, frq_df_group = frq_df_group_sender_receiver, pb_df_group = pb_df_group_sender_receiver, rel_abundance_df = rel_abundance_df_sender_receiver))
}

#' @title combine_sender_receiver_de
#'
#' @description \code{combine_sender_receiver_de}  Combine Muscat differential expression output for senders and receivers by linkgin ligands to receptors based on the prior knowledge ligand-receptor network.
#' @usage combine_sender_receiver_de(sender_de, receiver_de, senders_oi, receivers_oi, lr_network)
#'
#' @inheritParams multi_nichenet_analysis
#' @inheritParams combine_sender_receiver_info_ic
#' @param sender_de Differential expression analysis output for the sender cell types. `de_output_tidy` slot of the output of `perform_muscat_de_analysis`.
#' @param receiver_de Differential expression analysis output for the receiver cell types. `de_output_tidy` slot of the output of `perform_muscat_de_analysis`.
#'
#' @return Data frame combining Muscat DE output for sender and receiver linked to each other through joining by the ligand-receptor network.
#'
#' @import dplyr
#' @import muscat
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' batches = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    contrasts = contrasts_oi)
#'sender_receiver_de = combine_sender_receiver_de(
#'  sender_de = celltype_de$de_output_tidy,
#'  receiver_de = celltype_de$de_output_tidy,
#'  senders_oi = senders_oi,
#'  receivers_oi = receivers_oi,
#'  lr_network = lr_network)
#' }
#'
#' @export
#'
combine_sender_receiver_de = function(sender_de, receiver_de, senders_oi, receivers_oi, lr_network){

  requireNamespace("dplyr")
  
  de_output_tidy_sender = sender_de %>% dplyr::filter(cluster_id %in% senders_oi)
  de_output_tidy_receiver = receiver_de %>% dplyr::filter(cluster_id %in% receivers_oi)
  
  de_output_tidy_sender = de_output_tidy_sender %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast) %>% dplyr::filter(cluster_id %in% senders_oi) %>% dplyr::rename(ligand = gene, lfc_ligand = logFC, p_val_ligand = p_val,  p_adj_ligand = p_adj, sender = cluster_id)
  de_output_tidy_receiver =  de_output_tidy_receiver %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast) %>% dplyr::filter(cluster_id %in% receivers_oi) %>% dplyr::rename(receptor = gene, lfc_receptor = logFC, p_val_receptor = p_val,  p_adj_receptor = p_adj, receiver = cluster_id)

  de_tbl_sender_receiver = de_output_tidy_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(de_output_tidy_receiver, by = c("receptor","contrast"))
  de_tbl_sender_receiver = de_tbl_sender_receiver %>% dplyr::mutate(ligand_receptor_lfc_avg = (lfc_receptor + lfc_ligand)/2) %>% dplyr::arrange(-ligand_receptor_lfc_avg) %>% dplyr::select(contrast, sender, receiver, ligand, receptor, lfc_ligand, lfc_receptor, ligand_receptor_lfc_avg, p_val_ligand, p_adj_ligand, p_val_receptor, p_adj_receptor) %>% dplyr::distinct()

  return(de_tbl_sender_receiver)
}
