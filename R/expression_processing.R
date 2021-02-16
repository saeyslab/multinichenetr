#' @title get_muscat_exprs_frac
#'
#' @description \code{get_muscat_exprs_frac}  XXXX
#' @usage get_muscat_exprs_frac(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA")
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
get_muscat_exprs_frac = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA"){

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

  celltypes = sce$cluster_id %>% unique() %>% as.character()

  frq_lists = celltypes %>% lapply(function(celltype_oi, frq){
    frq_celltype = frq@assays@data[[celltype_oi]]

    frq_celltype_samples = frq_celltype[,samples]
    frq_celltype_samples = frq_celltype_samples %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(sample, fraction_sample, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)

    frq_celltype_groups = frq_celltype[,groups]
    frq_celltype_groups = frq_celltype_groups %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tidyr::gather(group, fraction_group, -gene) %>% tibble::as_tibble() %>% dplyr::mutate(celltype = celltype_oi)

    return(list(frq_celltype_samples = frq_celltype_samples, frq_celltype_groups = frq_celltype_groups))
  },frq) %>% magrittr::set_names(sce$cluster_id %>% unique())

  frq_celltype_samples = frq_lists %>% purrr::map("frq_celltype_samples") %>% dplyr::bind_rows()
  frq_celltype_groups = frq_lists %>% purrr::map("frq_celltype_groups") %>% dplyr::bind_rows()
  rm(frq_lists)

  return(list(frq_celltype_samples = frq_celltype_samples, frq_celltype_groups = frq_celltype_groups))
}

#' @title get_muscat_exprs_avg
#'
#' @description \code{get_muscat_exprs_avg}  XXXX
#' @usage get_muscat_exprs_avg(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA")
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
get_muscat_exprs_avg = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi_sce = "RNA"){

  # convert seurat to SCE object
  sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = assay_oi_sce)

  # prepare SCE for the muscat pseudobulk analysis
  sce$id = sce[[sample_id]]
  sce = muscat::prepSCE(sce,
                        kid = celltype_id, # subpopulation assignments
                        gid = group_id,  # group IDs (ctrl/stim)
                        sid = "id",   # sample IDs (ctrl/stim.1234)
                        drop = FALSE)  #

  sce = scater::computeLibraryFactors(sce)
  sce = scater::logNormCounts(sce)
  # scater::calculateAverage()
  # sce = scater::calculateCPM(sce)
  # scater::normalizeCounts()
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
#' @description \code{fix_frq_df}  XXXX
#' @usage fix_frq_df(seurat_obj, frq_celltype_samples)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @param frq_celltype_samples XXXX
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
fix_frq_df = function(seurat_obj, frq_celltype_samples){
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
#' @description \code{get_avg_frac_exprs_abund}  XXXX
#' @usage get_avg_frac_exprs_abund(seurat_obj, sample_id, celltype_id, group_id, assay_oi = "RNA")
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
get_avg_frac_exprs_abund = function(seurat_obj, sample_id, celltype_id, group_id, assay_oi = "RNA"){
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
#' @description \code{process_info_to_ic}  XXXX
#' @usage process_info_to_ic(info_object, ic_type = "sender", lr_network)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @param info_object XXX
#' @param ic_type XXX
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
process_info_to_ic = function(info_object, ic_type = "sender", lr_network){

  ligands = lr_network %>% pull(ligand) %>% unique()
  receptors = lr_network %>% pull(receptor) %>% unique()

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

#' @title combine_sender_receiver_info_ic
#'
#' @description \code{combine_sender_receiver_info_ic}  XXXX
#' @usage combine_sender_receiver_info_ic(sender_info, receiver_info, senders_oi, receivers_oi, lr_network)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @param sender_info XXX
#' @param receiver_info XXX
#' @param senders_oi XXX
#' @param receivers_oi XXX
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
combine_sender_receiver_info_ic = function(sender_info, receiver_info, senders_oi, receivers_oi, lr_network){

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
#' @description \code{combine_sender_receiver_de}  XXXX
#' @usage combine_sender_receiver_de(sender_de, receiver_de, senders_oi, receivers_oi, lr_network)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @inheritParams combine_sender_receiver_info_ic
#' @param sender_de XXX
#' @param receiver_de XXX
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
combine_sender_receiver_de = function(sender_de, receiver_de, senders_oi, receivers_oi, lr_network){

  de_output_tidy_sender = muscat::resDS(sender_de$sce, sender_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()
  de_output_tidy_receiver = muscat::resDS(receiver_de$sce, receiver_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()

  de_output_tidy_sender = de_output_tidy_sender %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj.loc, contrast) %>% dplyr::filter(cluster_id %in% senders_oi) %>% dplyr::rename(ligand = gene, lfc_ligand = logFC, p_val_ligand = p_val,  p_adj_ligand = p_adj.loc, sender = cluster_id)
  de_output_tidy_receiver =  de_output_tidy_receiver %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj.loc, contrast) %>% dplyr::filter(cluster_id %in% receivers_oi) %>% dplyr::rename(receptor = gene, lfc_receptor = logFC, p_val_receptor = p_val,  p_adj_receptor = p_adj.loc, receiver = cluster_id)

  de_tbl_sender_receiver = de_output_tidy_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(de_output_tidy_receiver, by = c("receptor","contrast"))
  de_tbl_sender_receiver = de_tbl_sender_receiver %>% dplyr::mutate(ligand_receptor_lfc_avg = (lfc_receptor + lfc_ligand)/2) %>% dplyr::arrange(-ligand_receptor_lfc_avg) %>% dplyr::select(contrast, sender, receiver, ligand, receptor, lfc_ligand, lfc_receptor, ligand_receptor_lfc_avg, p_val_ligand, p_adj_ligand, p_val_receptor, p_adj_receptor) %>% dplyr::distinct()

  return(de_tbl_sender_receiver)
}
