#' @title perform_muscat_de_analysis
#'
#' @description \code{perform_muscat_de_analysis} Perform differential expression analysis via Muscat - Pseudobulking approach - FOR ONE CELLTYPE.
#' @usage perform_muscat_de_analysis(sce, sample_id, celltype_id, group_id, batches, covariates, contrasts, expressed_df, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10)
#'
#' @inheritParams multi_nichenet_analysis
#' @param expressed_df tibble with three columns: gene, celltype, expressed; this data frame indicates which genes can be considered as expressed in each cell type. 
#' @param contrasts String indicating the contrasts of interest (= which groups/conditions will be compared) for the differential expression and MultiNicheNet analysis. 
#' We will demonstrate here a few examples to indicate how to write this. Check the limma package manuals for more information about defining design matrices and contrasts for differential expression analysis.
#' If wanting to compare group A vs B: `contrasts_oi = c("'A-B'")`
#' If wanting to compare group A vs B & B vs A: `contrasts_oi = c("'A-B','B-A'")`
#' If wanting to compare group A vs B & A vs C & A vs D: `contrasts_oi = c("'A-B','A-C', 'A-D'")`
#' If wanting to compare group A vs B and C: `contrasts_oi = c("'A-(B+C)/2'")`
#' If wanting to compare group A vs B, C and D: `contrasts_oi = c("'A-(B+C+D)/3'")`
#' If wanting to compare group A vs B, C and D & B vs A,C,D: `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")`
#' Note that the groups A, B, ... should be present in the meta data column 'group_id'.
#' @return List with output of the differential expression analysis in 1) default format(`muscat::pbDS()`), and 2) in a tidy table format (`muscat::resDS()`).
#'
#' @import dplyr
#' @import muscat
#' @importFrom SummarizedExperiment assayNames colData
#' @importFrom S4Vectors metadata
#' @importFrom limma makeContrasts
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "CAF"
#' batches = NA
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' frq_list = get_frac_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    covariates = covariates,
#'    contrasts = contrasts_oi,
#'    expressed_df = frq_list$expressed_df)
#'}
#'
#' @export
#'
#'
perform_muscat_de_analysis = function(sce, sample_id, celltype_id, group_id, batches, covariates, contrasts, expressed_df, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10,...){
  requireNamespace("dplyr")

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
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,celltype_id]))))){
      stop("All the cell type labels in SummarizedExperiment::colData(sce)[,celltype_id] should be syntactically valid R names - see make.names")
    }
  }

  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,group_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,group_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,group_id]))))){
      stop("All the group/condition labels in SummarizedExperiment::colData(sce)[,group_id] should be syntactically valid R names - see make.names")
    }
  }
  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,sample_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,sample_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,sample_id]))))){
      stop("All the sample_id labels in SummarizedExperiment::colData(sce)[,sample_id] should be syntactically valid R names - see make.names")
    }
  }
  
  
  if(!is.character(contrasts)){
    stop("contrasts should be a character vector")
  }

  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = SummarizedExperiment::colData(sce)[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts, "'") %>% unlist() %>% unique() %>%
    # stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",","," ,", ", ")) %>% unlist() %>% unique()
  conditions_oi = conditions_oi[is.na(suppressWarnings(as.numeric(conditions_oi)))]
  
  if(length(contrasts) != 1 | !is.character(contrasts)){
    stop("contrasts should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }
  
  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_simplified = stringr::str_split(contrasts, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()
  
  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }

  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  }
  if(length(covariates) > 1){
    covariates_present = TRUE
    if (sum(covariates %in% colnames(SummarizedExperiment::colData(sce))) != length(covariates) ) {
      stop("covariates should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  } else {
    if(!is.na(covariates)){
      covariates_present = TRUE
      if (sum(covariates %in% colnames(SummarizedExperiment::colData(sce))) != length(covariates) ) {
        stop("covariates should be NA or all present as column name(s) in the metadata dataframe of sce")
      }
    } else {
      covariates_present = FALSE
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

  requireNamespace("dplyr")
  

  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  # drop all other SummarizedExperiment::colData columns ----------------- change to false

  # test to see whether sample_ids are unique
  if (sum(table(sce$sample_id, sce$group_id) %>% apply(1, function(row_oi){sum(row_oi > 0)}) > 1) > 0){
    stop("One or more of your sample_ids belongs to more than one group/condition of interest. Please make sure that all sample_ids are uniquely divided over your groups/conditions.")
  } 
  
  pb = muscat::aggregateData(sce,
                             assay = assay_oi_pb, fun = fun_oi_pb,
                             by = c("cluster_id", "sample_id"),...)
  
  if(assay_oi_pb == "counts"){
    libsizes = colSums(SummarizedExperiment::assay(pb))
    if (!isTRUE(all(libsizes == floor(libsizes)))) {
      warning("non-integer library sizes: are you sure you are working with raw counts?")
    }
  }

  # prepare the experiment info (ei) table if batches present
  if(length(batches) > 1){
    batches_present = TRUE
  } else {
    if(!is.na(batches)){
      batches_present = TRUE
    } else {
      batches_present = FALSE

    }
  }

  if(batches_present){
    extra_metadata = SummarizedExperiment::colData(sce)  %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batches)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
  } else {
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor) 
  } 
  if('sample_id' != sample_id){
    extra_metadata$sample_id = extra_metadata[[sample_id]]
  }
  ei = S4Vectors::metadata(sce)$experiment_info

  ei = ei %>%  dplyr::inner_join(extra_metadata, by = "sample_id")

  
  # prepare the experiment info (ei) table if covariates present
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
    extra_metadata = SummarizedExperiment::colData(sce)  %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(covariates)) %>% dplyr::distinct() ## no factorization of the continuous covariates!
  } else {
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor) 
  } 
  if('sample_id' != sample_id){
    extra_metadata$sample_id = extra_metadata[[sample_id]]
  }
  if(is.null(ei)){
    ei = S4Vectors::metadata(sce)$experiment_info
  }
  
  ei = ei %>%  dplyr::inner_join(extra_metadata, by = "sample_id")
  
  # prepare the design and contrast matrix for the muscat DE analysis
  
  if(batches_present == TRUE & covariates_present == FALSE){
    
    batches_string = paste0("ei$",batches) %>% paste(collapse = " + ")
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id + ", batches_string, " ) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id), make.names(colnames(design)[(length(levels(ei$group_id))+1):length(colnames(design))])) )

  } else if (batches_present == FALSE & covariates_present == TRUE) {
    
    covariates_string = paste0("ei$",covariates) %>% paste(collapse = " + ")
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id + ", covariates_string, " ) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id), make.names(colnames(design)[(length(levels(ei$group_id))+1):length(colnames(design))])) )

  } else if (batches_present == TRUE & covariates_present == TRUE) {
    
    batches_string = paste0("ei$",batches) %>% paste(collapse = " + ")
    covariates_string = paste0("ei$",covariates) %>% paste(collapse = " + ")
    
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id + ", batches_string, " + ", covariates_string, " ) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id), make.names(colnames(design)[(length(levels(ei$group_id))+1):length(colnames(design))])) )

  } else if (batches_present == FALSE & covariates_present == FALSE) {
    
    design = eval(parse(text=paste("model.matrix(~ 0 + ei$group_id) ",sep="")))
    dimnames(design) = list(ei$sample_id, c(levels(ei$group_id)))

  }

  contrast = eval(parse(text=paste("limma::makeContrasts(", contrasts, ",levels=design)",sep="")))

  # filter genes
  celltype_oi = SummarizedExperiment::assayNames(pb) %>% .[1]
  genes_to_keep = expressed_df %>% filter(celltype == celltype_oi & expressed == TRUE) %>% pull(gene) %>% unique()
  pb = pb[genes_to_keep, ]
  
  # run DS analysis
  res = muscat::pbDS(pb, method = de_method_oi , design = design, contrast = contrast, min_cells = min_cells, verbose = FALSE, filter = "none",...)
  de_output_tidy = muscat::resDS(sce, res, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble() %>% dplyr::rename(p_adj = p_adj.glb)
  
  # # check which cell types were excluded 
  celltypes = SummarizedExperiment::assayNames(pb)
  names(celltypes) = celltypes
  
  excluded_celltypes = celltypes %>% generics::setdiff(de_output_tidy$cluster_id) %>% unique()
  
  if(length(excluded_celltypes) > 0){
    print("excluded cell types are:")
    print(excluded_celltypes)
    print("These celltypes are not considered in the analysis. After removing samples that contain less cells than the required minimal, some groups/conditions don't have 2 or more samples anymore or there are not all batch/categorical.covariate levels are present in all groups/conditions. As a result the analysis cannot be run. To solve this: decrease the number of min_cells or change your group_id and pool all samples that belong to groups that are not of interest or do not consider the batch/categorical.covariate! ")
  }
  if(length(excluded_celltypes) == length(celltypes)){
    print("None of the cell types passed the check. This might be due to 2 reasons. 1) no cell type has enough cells in >=2 samples per group. 2) problem in batch or categorical.covariate definition: not all levels of your batch/categorical.covariate are in each group - Also for groups not included in your contrasts! To solve this:pool all samples that belong to groups that are not of interest or do not consider the batch/categorical.covariate!")
  }
  
  return(list(de_output = res, de_output_tidy = de_output_tidy))
}
#' @title Get empirical p-values and adjusted p-values.
#'
#' @description \code{p.adjust_empirical} Get emperical p-values and adjusted p-values. Credits to Jeroen Gillis (cf satuRn package)
#' @usage p.adjust_empirical(pvalues, tvalues, plot = FALSE, celltype = NULL, contrast = NULL)
#'
#' @param pvalues Vector of original p-values
#' @param tvalues Vector of original t-values: in case of DE analysis: log fold changes
#' @param plot TRUE or FALSE (default): should we plot the z-score distribution?
#' @param celltype NULL, or name of the cell type of interest - this will be added to the plot title if plot = TRUE
#' @param contrast NULL, or name of the contrast of interest - thhis will be added to the plot title if plot = TRUE
#' @return List with empirical p-values and adjusted p-values + ggplot output and estimated delta and sigma.
#'
#' @importFrom stats p.adjust dnorm glm median pnorm poisson poly qnorm quantile
#' @importFrom graphics hist lines
#' @import locfdr
#' @import ggplot2
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
#' de_output_tidy = muscat::resDS(celltype_de$sce, celltype_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()
#' emp_res = p.adjust_empirical(de_output_tidy %>% pull(p_val), de_output_tidy  %>% pull(p_val), plot = T)
#'}
#'
#' @export
#'
p.adjust_empirical <- function (pvalues, tvalues, plot = FALSE, celltype = NULL, contrast = NULL) 
{
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  zvalues <- qnorm(pvalues/2) * sign(tvalues) # sign(t) for you will be sign(LFC)
  
  zvalues_mid <- zvalues[abs(zvalues) < 10] 
  zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)]
  # avoid numeric issues by removing extreme z-scores for estimating the empirical null
  # the good thing is that these are only removed for estimating the null;
  # you will still have results for them in the end
  
  # compute empricial null distribution -> from locfdr
  #### start
  N <- length(zvalues_mid)
  b <- 4.3 * exp(-0.26 * log(N, 10))
  med <- median(zvalues_mid)
  sc <- diff(quantile(zvalues_mid)[c(2, 4)])/(2 * qnorm(0.75))
  mlests <- locfdr:::locmle(zvalues_mid, xlim = c(med, b * 
                                                    sc))
  lo <- min(zvalues_mid)
  up <- max(zvalues_mid)
  bre <- 120
  breaks <- seq(lo, up, length = bre)
  zzz <- pmax(pmin(zvalues_mid, up), lo)
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  x <- (breaks[-1] + breaks[-length(breaks)])/2
  sw <- 0
  X <- cbind(1, poly(x, df = 7))
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  y <- zh$counts
  f <- glm(y ~ poly(x, df = 7), poisson)$fit
  Cov.in <- list(x = x, X = X, f = f, sw = sw)
  ml.out <- locfdr:::locmle(zvalues_mid, xlim = c(mlests[1], 
                                                  b * mlests[2]), d = mlests[1], s = mlests[2], Cov.in = Cov.in)
  mlests <- ml.out$mle # MLEs for the empirical null; same interpretation as
  # the MLE I discussed in the locfdr plot
  #### end
  
  # "Correct" the z-scores; based on the empirical null MLE, we force the bulk
  # of the test statistics to follow a z-distribution. This is done by taking
  # the original z-scores, substractin with the MLE for delta and dividing by
  # sigma
  zval_empirical <- (zvalues - mlests[1])/mlests[2]
  
  # based on the new z-scores, compute new p-values
  pval_empirical <- 2 * pnorm(-abs(zval_empirical), mean = 0, 
                              sd = 1)
  
  # make a plot similar to the one you saw above from locfdr
  if (plot) {
    zval_empirical <- zval_empirical[!is.na(zval_empirical)]
    lo <- min(zval_empirical)
    up <- max(zval_empirical)
    lo <- min(lo, -1 * up)
    up <- max(up, -1 * lo) # to make plot symmetrical
    bre <- 120
    breaks <- seq(lo, up, length = bre)
    zzz <- pmax(pmin(zval_empirical, up), lo) # do not plot extreme zvals
    
    gg_data <- as.data.frame(zzz)
    gg_plot <- ggplot(data = gg_data, aes(x=zzz)) +
      geom_histogram(aes(y=..density..),
                     breaks=breaks,
                     color="black",
                     fill="grey84") +
      theme_bw() +
      ggtitle("Empirical Z-scores") +
      theme(plot.title = element_text(size = 14, face = "bold")) +
      xlab("Empirical Z-scores")
    
    # add standard-normal density
    xfit <- seq(min(zzz), max(zzz), length = 4000)
    yfit <- dnorm(xfit/mlests[3], mean = 0, sd = 1)
    dens_theor <- as.data.frame(yfit)
    colnames(dens_theor) <- "standardNormal"
    dens_theor$range <- xfit
    gg_plot <- gg_plot +
      geom_line(data = dens_theor, 
                aes(x=range, y=standardNormal),
                color="darkgreen",
                lwd=1) + ggtitle(paste0("Empirical distribution of z-scores\n",celltype, " : ", contrast))
  } else{
    gg_plot <- NA
  }
  
  # return FDR values; these are normal FDR, not local FDR!
  FDR <- p.adjust(pval_empirical, method = "BH")
  # newList <- list(pval = pval_empirical, FDR = FDR, plot_output = plot_output)
  newList <- list(plot_output = gg_plot,
                  delta = mlests[1], # add delta estimate to output
                  sigma = mlests[2], # add sigma estimate to output
                  pval = pval_empirical, 
                  FDR = FDR)
  return(newList)
}
#' @title Add empirical p-values and adjusted p-values to a subset of the DE output table.
#' @description \code{get_FDR_empirical}  Add empirical p-values and adjusted p-values to a subset of the DE output table. This is the function that works under  the hood of `add_empirical_pval_fdr`. Credits to Jeroen Gillis (cf satuRn package)
#' @usage get_FDR_empirical(de_output_tidy, cluster_id_oi, contrast_oi, plot = FALSE)
#'
#' @param de_output_tidy Data frame of DE results, containing at least the following columns: cluster_id, contrast, p_val, logFC.
#' @param cluster_id_oi Indicate which celltype DE results should be filtered for.
#' @param contrast_oi Indicate which contrast DE results should be filtered for.
#' @param plot TRUE or FALSE (default): should we plot the z-score distribution?
#' @return de_output_tidy dataframe (for the celltype and contrast of interest) with two columns added: p_emp and p_adj_emp.
#'
#' @export
#'
get_FDR_empirical = function(de_output_tidy, cluster_id_oi, contrast_oi, plot = FALSE){
  de_oi = de_output_tidy %>% filter(cluster_id == cluster_id_oi & contrast == contrast_oi)
  emp_res = p.adjust_empirical(de_oi %>% pull(p_val), de_oi %>% pull(logFC), plot = plot, celltype = cluster_id_oi, contrast = contrast_oi)
  de_oi = de_oi %>% mutate(p_emp = emp_res$pval, p_adj_emp = emp_res$FDR)
}

#' @title Add empirical p-values and adjusted p-values to the DE output table.
#'
#' @description \code{add_empirical_pval_fdr} Add empirical p-values and adjusted p-values to the DE output of Muscat. Credits to Jeroen Gillis (cf satuRn package)
#' @usage add_empirical_pval_fdr(de_output_tidy, plot = FALSE)
#'
#' @inheritParams get_FDR_empirical
#' @return de_output_tidy dataframe with two columns added: p_emp and p_adj_emp.
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
#' de_output_tidy = muscat::resDS(celltype_de$sce, celltype_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()
#' de_output_tidy = add_empirical_pval_fdr(de_output_tidy)
#'}
#'
#' @export
#'
add_empirical_pval_fdr = function(de_output_tidy, plot = FALSE){
  requireNamespace("dplyr")
  
  all_celltypes = de_output_tidy$cluster_id %>% unique()
  all_contrasts = de_output_tidy$contrast %>% unique()
    
  de_output_tidy_new = all_celltypes %>% lapply(function(cluster_id_oi, de_output_tidy){
    all_contrasts %>% lapply(function(contrast_oi, de_output_tidy) {
      de_output_subset = get_FDR_empirical(de_output_tidy, cluster_id_oi,  contrast_oi, plot = plot)
    },de_output_tidy) %>% bind_rows()
  },de_output_tidy) %>% bind_rows()
  return(de_output_tidy_new)
}
#' @title Get diagnostic plots of the empirical null.
#' @description \code{get_FDR_empirical_plots}  Get diagnostic plots of the empirical null.. This is the function that works under  the hood of `get_FDR_empirical_plots_all`. Credits to Jeroen Gillis (cf satuRn package)
#' @usage get_FDR_empirical_plots(de_output_tidy, cluster_id_oi, contrast_oi)
#'
#' @param de_output_tidy Data frame of DE results, containing at least the following columns: cluster_id, contrast, p_val, logFC.
#' @param cluster_id_oi Indicate which celltype DE results should be filtered for.
#' @param contrast_oi Indicate which contrast DE results should be filtered for.
#' @return plot object
#'
#' @export
#'
get_FDR_empirical_plots = function(de_output_tidy, cluster_id_oi, contrast_oi){
  de_oi = de_output_tidy %>% filter(cluster_id == cluster_id_oi & contrast == contrast_oi)
  emp_res = p.adjust_empirical(de_oi %>% pull(p_val), de_oi %>% pull(logFC), plot = TRUE, celltype = cluster_id_oi, contrast = contrast_oi)
  emp_res$plot_output
}
#' @title Get diagnostic plots of the empirical null.
#' @description \code{get_FDR_empirical_plots_all}  Get diagnostic plots of the empirical null. Credits to Jeroen Gillis (cf satuRn package)
#' @usage get_FDR_empirical_plots_all(de_output_tidy)
#'
#' @param de_output_tidy Data frame of DE results, containing at least the following columns: cluster_id, contrast, p_val, logFC.
#' @return list of plots
#'
#' @importFrom magrittr set_names
#' @export
#'
get_FDR_empirical_plots_all = function(de_output_tidy){
  requireNamespace("dplyr")
  
  all_celltypes = de_output_tidy$cluster_id %>% unique()
  all_contrasts = de_output_tidy$contrast %>% unique()
  
  all_plots = all_celltypes %>% lapply(function(cluster_id_oi, de_output_tidy){
    all_contrasts %>% lapply(function(contrast_oi, de_output_tidy) {
      de_output_subset = get_FDR_empirical_plots(de_output_tidy, cluster_id_oi,  contrast_oi)
    },de_output_tidy) %>% magrittr::set_names(all_contrasts)
  },de_output_tidy) %>% magrittr::set_names(all_celltypes)
  return(all_plots %>% unlist(recursive = F))
  # plots = get_FDR_empirical_plots_all(output$receiver_de)
  # for(plot_oi in plots) {
  #   #dev.control("enable")
  #   replayPlot(plot_oi)
  #   dev.control("inhibit")
  # }
  
}

