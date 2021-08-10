#' @title get_muscat_exprs_frac
#'
#' @description \code{get_muscat_exprs_frac} Calculate sample- and group-average of fraction of cells in a cell type expressing a gene.
#' @usage get_muscat_exprs_frac(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA")
#'
#' @inheritParams multi_nichenet_analysis_combined
#'
#' @return List with two dataframes: one with fraction of cells in a cell type expressing a gene, averaged per sample; and one averaged per group.
#'
#' @import Seurat
#' @import dplyr
#' @import muscat
#' @import tibble
#' @import tidyr
#' @importFrom purrr map
#' @importFrom magrittr set_names
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' muscat_exprs_frac = get_muscat_exprs_frac(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_muscat_exprs_frac = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA"){

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
                        drop = FALSE)  #

  samples = sce$sample_id %>% unique() %>% as.character()
  groups = sce$group_id  %>% unique() %>% as.character()

  frq = muscat::calcExprFreqs(sce, assay = "counts", th = 0) # gives NaN sometimes...
  all_genes = rownames(frq)
  
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
#' @usage get_muscat_exprs_avg(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA")
#'
#' @inheritParams multi_nichenet_analysis_combined
#'
#' @return Data frame with average gene expression per sample and per group.
#'
#' @import Seurat
#' @import dplyr
#' @import muscat
#' @import tibble
#' @import tidyr
#' @importFrom scater logNormCounts
#' @importFrom purrr map
#' @importFrom scran quickCluster computeSumFactors
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' muscat_exprs_avg = get_muscat_exprs_avg(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_muscat_exprs_avg = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA"){

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
                        drop = FALSE)  #

  # sce = scater::computeLibraryFactors(sce)
  # sce = scater::logNormCounts(sce)
  # scater::calculateAverage()
  # sce = scater::calculateCPM(sce)
  # scater::normalizeCounts()
  
  set.seed(1919) # seed for reproducibility of the scran::quickCluster function (see https://bioconductor.org/books/release/OSCA/normalization.html)
  clusters = scran::quickCluster(sce)
  sce = scran::computeSumFactors(sce, clusters=clusters)
  sce = scater::logNormCounts(sce)
  
  avg = muscat::aggregateData(sce, assay = "logcounts", fun = "mean", by = c("cluster_id", "sample_id"))

  avg_df = sce$cluster_id %>% unique() %>% lapply(function(celltype_oi, avg){
    avg_celltype = avg@assays@data[[celltype_oi]]

    avg_celltype_samples = avg_celltype[,sce$sample_id %>% unique()]
    avg_celltype_samples = avg_celltype_samples %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, average_sample, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)

  },avg) %>% dplyr::bind_rows()

  return(avg_df)
}

#' @title fix_frq_df
#'
#' @description \code{fix_frq_df}  Fix muscat-bug in fraction calculation in case that expression fraction would be NA / NaN. Change NA to 0.
#' @usage fix_frq_df(seurat_obj, frq_celltype_samples)
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @param frq_celltype_samples Sample-average data frame output of `get_muscat_exprs_frac`
#'
#' @return Fixed data frame with fraction of cells expressing a gene.
#'
#' @import Seurat
#' @import dplyr
#' @importFrom magrittr set_names
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' frq_df = get_muscat_exprs_frac(seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, assay_oi_sce = assay_oi) %>% .$frq_celltype_samples
#' if(nrow(frq_df %>% dplyr::filter(is.na(fraction_sample))) > 0 | nrow(frq_df %>% dplyr::filter(is.nan(fraction_sample))) > 0) {
#'   frq_df = fix_frq_df(seurat_obj, frq_df)
#'   }
#' }
#'
#' @export
#'
fix_frq_df = function(seurat_obj, frq_celltype_samples){
  
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  
  genes = seurat_obj@assays$RNA@data %>% rownames()
  gene_mapping = genes %>% magrittr::set_names(seq(length(genes)))

  frq_celltype_samples_OK = frq_celltype_samples %>% dplyr::filter(gene %in% genes)

  frq_celltype_samples_FIX = frq_celltype_samples %>% dplyr::filter(!gene %in% genes)

  frq_celltype_samples_FIX = frq_celltype_samples_FIX %>% dplyr::mutate(gene = gene_mapping[gene])

  frq_celltype_samples_FIX = frq_celltype_samples_FIX %>% dplyr::mutate(fraction_sample = 0)

  frq_celltype_samples = frq_celltype_samples_OK %>% dplyr::bind_rows(frq_celltype_samples_FIX)

  return(frq_celltype_samples)

}

#' @title get_avg_frac_exprs_abund
#'
#' @description \code{get_avg_frac_exprs_abund}  Calculate the average and fraction of expression of each gene per sample and per group. Calculate relative abundances of cell types as well.
#' @usage get_avg_frac_exprs_abund(seurat_obj, sample_id, celltype_id, group_id, assay_oi = "RNA")
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @param assay_oi Indicates which assay of the Seurat object should be used. Default: "RNA". See: `Seurat::as.SingleCellExperiment`.
#' 
#' @return List containing data frames with average and fraction of expression per sample and per group, and relative cell type abundances as well.
#'
#' @import Seurat
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_frac_exprs_abund(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_avg_frac_exprs_abund = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi = "RNA"){
  
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  
  # input checks
  
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
    is_make_names = levels(seurat_obj@meta.data[,celltype_id]) == make.names(levels(seurat_obj@meta.data[,celltype_id]))
    if(sum(is_make_names) != length(levels(seurat_obj@meta.data[,celltype_id]))){
      stop("The levels of the factor seurat_obj@meta.data[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  }
  if(is.factor(seurat_obj@meta.data[,group_id])){
    is_make_names = levels(seurat_obj@meta.data[,group_id]) == make.names(levels(seurat_obj@meta.data[,group_id]))
    if(sum(is_make_names) != length(levels(seurat_obj@meta.data[,group_id]))){
      stop("The levels of the factor seurat_obj@meta.data[,group_id] should be a syntactically valid R names - see make.names")
    }
  }
  if(is.factor(seurat_obj@meta.data[,sample_id])){
    is_make_names = levels(seurat_obj@meta.data[,sample_id]) == make.names(levels(seurat_obj@meta.data[,sample_id]))
    if(sum(is_make_names) != length(levels(seurat_obj@meta.data[,sample_id]))){
      stop("The levels of the factor seurat_obj@meta.data[,sample_id] should be a syntactically valid R names - see make.names")
    }
  }
  if(!is.character(assay_oi)){
    stop("assay_oi should be a character vector")
  } else {
    if(assay_oi != "RNA"){
      warning("are you sure you don't want to use the RNA assay?")
    }
  }
  
  ## calculate averages, fractions, relative abundance of a cell type in a group

  # calculate average expression
  avg_df = get_muscat_exprs_avg(seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, assay_oi_sce = assay_oi)

  # calculate fraction of expression
  frq_df = get_muscat_exprs_frac(seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, assay_oi_sce = assay_oi) %>% .$frq_celltype_samples

  # check whether something needs to be fixed
  if(nrow(avg_df %>% dplyr::filter(is.na(average_sample))) > 0 | nrow(avg_df %>% dplyr::filter(is.nan(average_sample))) > 0) {
    warning("There are some genes with NA average expression.")
  }
  if(nrow(frq_df %>% dplyr::filter(is.na(fraction_sample))) > 0 | nrow(frq_df %>% dplyr::filter(is.nan(fraction_sample))) > 0) {
    warning("There are some genes with NA fraction of expression. This is the result of the muscat function `calcExprFreqs` which will give NA when there are no cells of a particular cell type in a particular group. As a temporary fix, we give all these genes an expression fraction of 0 in that group for that cell type")
    frq_df = fix_frq_df(seurat_obj, frq_df)
  }

  # prepare grouping to get group averages
  metadata = seurat_obj@meta.data
  if('sample_id' != sample_id){
    metadata$sample_id = metadata[[sample_id]]
  }
  if('group_id' != sample_id){
    metadata$group_id = metadata[[group_id]]
  }
  if('celltype_id' != celltype_id){
    metadata$celltype_id = metadata[[celltype_id]]
  }

  grouping_df = metadata %>% dplyr::select(sample_id, group_id) %>% tibble::as_tibble() %>% dplyr::distinct() %>% dplyr::rename(sample = sample_id, group = group_id)

  avg_df_group = avg_df %>% dplyr::inner_join(grouping_df) %>% dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(average_group = mean(average_sample))
  frq_df_group = frq_df %>% dplyr::inner_join(grouping_df) %>% dplyr::group_by(group, celltype, gene) %>% dplyr::summarise(fraction_group = mean(fraction_sample))

  # calculate relative abundance
  n_celltypes = metadata$celltype_id %>% unique() %>% length()
  if(n_celltypes > 1){
    rel_abundance_celltype_vs_celltype = table(metadata$celltype_id, metadata$group_id) %>% apply(2, function(x){x/sum(x)})
    rel_abundance_celltype_vs_group = rel_abundance_celltype_vs_celltype %>% apply(1, function(x){x/sum(x)})
  } else {
    rel_abundance_celltype_vs_group = table(metadata$celltype_id, metadata$group_id) %>% apply(1, function(x){x/sum(x)})
  }

  # rel_ab_mean = rel_abundance_celltype_vs_group %>% apply(2, mean, na.rm = TRUE)
  # rel_ab_sd = rel_abundance_celltype_vs_group %>% apply(2, sd, na.rm = TRUE)
  # rel_ab_z = (rel_abundance_celltype_vs_group - rel_ab_mean) / rel_ab_sd
  # rel_abundance_df = rel_ab_z %>% data.frame() %>% tibble::rownames_to_column("group") %>% tidyr::gather(celltype, rel_abundance_scaled, -group) %>% tibble::as_tibble()
  rel_abundance_df = rel_abundance_celltype_vs_group %>% data.frame() %>% tibble::rownames_to_column("group") %>% tidyr::gather(celltype, rel_abundance_scaled, -group) %>% tibble::as_tibble() %>% dplyr::mutate(rel_abundance_scaled = scale_quantile_adapted(rel_abundance_scaled))
  return(list(avg_df = avg_df, frq_df = frq_df, avg_df_group = avg_df_group, frq_df_group = frq_df_group, rel_abundance_df = rel_abundance_df))

}

#' @title process_info_to_ic
#'
#' @description \code{process_info_to_ic}  Process cell type expression information into intercellular communication focused information. Only keep information of ligands for the sender cell type setting, and information of receptors for the receiver cell type.
#' @usage process_info_to_ic(info_object, ic_type = "sender", lr_network)
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @param info_object Output of `get_avg_frac_exprs_abund`
#' @param ic_type "sender" or "receiver": indicates whether we should keep ligands or receptors respectively.
#'
#' @return List with expression information of ligands (sender case) or receptors (receiver case) - similar to output of `get_avg_frac_exprs_abund`.
#'
#' @import dplyr
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
#' celltype_info = get_avg_frac_exprs_abund(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
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

    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, avg_ligand_group = average_group)
    frq_df_group = info_object$frq_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, fraction_ligand_group = fraction_group)

    rel_abundance_df = info_object$rel_abundance_df %>% dplyr::rename(sender = celltype, rel_abundance_scaled_sender = rel_abundance_scaled)
  }
  if(ic_type == "receiver"){
    avg_df = info_object$avg_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, avg_receptor = average_sample)
    frq_df = info_object$frq_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor = fraction_sample)

    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, avg_receptor_group = average_group)
    frq_df_group = info_object$frq_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor_group = fraction_group)

    rel_abundance_df = info_object$rel_abundance_df %>% dplyr::rename(receiver = celltype, rel_abundance_scaled_receiver = rel_abundance_scaled)
  }

  return(list(avg_df = avg_df, frq_df = frq_df, avg_df_group = avg_df_group, frq_df_group = frq_df_group, rel_abundance_df = rel_abundance_df))
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
#' library(Seurat)
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' metadata_abundance = seurat_obj@meta.data[,c(sample_id, group_id, celltype_id)]
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
#' @inheritParams multi_nichenet_analysis_combined
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
#' library(Seurat)
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_frac_exprs_abund(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = Idents(seurat_obj) %>% unique()
#' receivers_oi = Idents(seurat_obj) %>% unique()
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

  # combine relative abundances
  rel_abundance_df_sender = sender_info$rel_abundance_df %>% dplyr::filter(sender %in% senders_oi)
  rel_abundance_df_receiver = receiver_info$rel_abundance_df %>% dplyr::filter(receiver %in% receivers_oi)

  rel_abundance_df_sender_receiver = rel_abundance_df_sender %>%  dplyr::inner_join(rel_abundance_df_receiver, by = "group") %>% dplyr::mutate(sender_receiver_rel_abundance_avg = 0.5*(rel_abundance_scaled_sender  + rel_abundance_scaled_receiver))

  # return
  return(list(avg_df = avg_df_sender_receiver, frq_df = frq_df_sender_receiver, avg_df_group = avg_df_group_sender_receiver, frq_df_group = frq_df_group_sender_receiver, rel_abundance_df = rel_abundance_df_sender_receiver))
}

#' @title combine_sender_receiver_de
#'
#' @description \code{combine_sender_receiver_de}  Combine Muscat differential expression output for senders and receivers by linkgin ligands to receptors based on the prior knowledge ligand-receptor network.
#' @usage combine_sender_receiver_de(sender_de, receiver_de, senders_oi, receivers_oi, lr_network)
#'
#' @inheritParams multi_nichenet_analysis_combined
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
  
  de_output_tidy_sender = sender_de %>% filter(cluster_id %in% senders_oi)
  de_output_tidy_receiver = receiver_de %>% filter(cluster_id %in% receivers_oi)
  
  de_output_tidy_sender = de_output_tidy_sender %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast) %>% dplyr::filter(cluster_id %in% senders_oi) %>% dplyr::rename(ligand = gene, lfc_ligand = logFC, p_val_ligand = p_val,  p_adj_ligand = p_adj, sender = cluster_id)
  de_output_tidy_receiver =  de_output_tidy_receiver %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast) %>% dplyr::filter(cluster_id %in% receivers_oi) %>% dplyr::rename(receptor = gene, lfc_receptor = logFC, p_val_receptor = p_val,  p_adj_receptor = p_adj, receiver = cluster_id)

  de_tbl_sender_receiver = de_output_tidy_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(de_output_tidy_receiver, by = c("receptor","contrast"))
  de_tbl_sender_receiver = de_tbl_sender_receiver %>% dplyr::mutate(ligand_receptor_lfc_avg = (lfc_receptor + lfc_ligand)/2) %>% dplyr::arrange(-ligand_receptor_lfc_avg) %>% dplyr::select(contrast, sender, receiver, ligand, receptor, lfc_ligand, lfc_receptor, ligand_receptor_lfc_avg, p_val_ligand, p_adj_ligand, p_val_receptor, p_adj_receptor) %>% dplyr::distinct()

  return(de_tbl_sender_receiver)
}
