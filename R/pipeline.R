#' @title ms_mg_nichenet_analysis
#'
#' @description \code{ms_mg_nichenet_analysis}  XXXX
#' @usage ms_mg_nichenet_analysis(sender_receiver_separate = TRUE, ...)
#'
#' @param sender_receiver_separate XXXX
#' @param ... Arguments to `ms_mg_nichenet_analysis_separate` (default; when `sender_receiver_separate = TRUE`) or `ms_mg_nichenet_analysis_combined` (when `sender_receiver_separate = FALSE`)
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
ms_mg_nichenet_analysis = function(sender_receiver_separate = TRUE, ...){

  if(!is.logical(sender_receiver_separate)) {
    stop("The sender_receiver_separate argument should be TRUE or FALSE")
  }

  if(sender_receiver_separate == TRUE){
    output = ms_mg_nichenet_analysis_separate(...)
  } else {
    output = ms_mg_nichenet_analysis_combined(...)
  }
  return(output)
}

#' @title ms_mg_nichenet_analysis_separate
#'
#' @description \code{ms_mg_nichenet_analysis_separate}  XXXX
#' @usage ms_mg_nichenet_analysis_separate(
#' seurat_obj_receiver,seurat_obj_sender,celltype_id_receiver,celltype_id_sender,sample_id,group_id,lr_network,ligand_target_matrix,contrasts_oi,contrast_tbl,
#' prioritizing_weights = c("scaled_lfc_ligand" = 1, "scaled_p_val_ligand" = 1, "scaled_lfc_receptor" = 1, "scaled_p_val_receptor" = 1, "scaled_activity_scaled" = 1.5,
#' "scaled_activity" = 0.5,"scaled_avg_exprs_ligand" = 1,"scaled_avg_frq_ligand" = 1,"scaled_avg_exprs_receptor" = 1, "scaled_avg_frq_receptor" = 1,
#' "fraction_expressing_ligand_receptor" = 1,"scaled_abundance_sender" = 0, "scaled_abundance_receiver" = 0),
#' assay_oi_sce = "RNA",assay_oi_pb ="counts",fun_oi_pb = "sum",de_method_oi = "edgeR",min_cells = 10,logFC_threshold = 0.25,p_val_threshold = 0.05,frac_cutoff = 0.05,p_val_adj = FALSE,top_n_target = 250, verbose = TRUE)
#'
#' @param seurat_obj_receiver XXXX
#' @param seurat_obj_sender XXXX
#' @param celltype_id_receiver XXXX
#' @param celltype_id_sender XXXX
#' @param sample_id XXXX
#' @param group_id XXXX
#' @param lr_network XXXX
#' @param ligand_target_matrix XXXX
#' @param contrasts_oi XXXX
#' @param contrast_tbl XXXX
#' @param prioritizing_weights XXXX
#' @param assay_oi_sce XXXX
#' @param assay_oi_pb XXXX
#' @param fun_oi_pb XXXX
#' @param de_method_oi XXXX
#' @param min_cells XXXX
#' @param logFC_threshold XXXX
#' @param p_val_threshold XXXX
#' @param frac_cutoff XXXX
#' @param p_val_adj XXXX
#' @param top_n_target XXXX
#' @param verbose XXXX
#'
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
ms_mg_nichenet_analysis_separate = function(seurat_obj_receiver,
                                            seurat_obj_sender,
                                            celltype_id_receiver,
                                            celltype_id_sender,
                                            sample_id,
                                            group_id,
                                            lr_network,
                                            ligand_target_matrix,
                                            contrasts_oi,
                                            contrast_tbl,
                                            prioritizing_weights = c("scaled_lfc_ligand" = 1,
                                                                     "scaled_p_val_ligand" = 1,
                                                                     "scaled_lfc_receptor" = 1,
                                                                     "scaled_p_val_receptor" = 1,
                                                                     "scaled_activity_scaled" = 1.5,
                                                                     "scaled_activity" = 0.5,
                                                                     "scaled_avg_exprs_ligand" = 1,
                                                                     "scaled_avg_frq_ligand" = 1,
                                                                     "scaled_avg_exprs_receptor" = 1,
                                                                     "scaled_avg_frq_receptor" = 1,
                                                                     "fraction_expressing_ligand_receptor" = 1,
                                                                     "scaled_abundance_sender" = 0,
                                                                     "scaled_abundance_receiver" = 0),
                                            assay_oi_sce = "RNA",
                                            assay_oi_pb ="counts",
                                            fun_oi_pb = "sum",
                                            de_method_oi = "edgeR",
                                            min_cells = 10,
                                            logFC_threshold = 0.25,
                                            p_val_threshold = 0.05,
                                            frac_cutoff = 0.05,
                                            p_val_adj = FALSE,
                                            top_n_target = 250, verbose = TRUE){


  # input checks

  if (class(seurat_obj_receiver) != "Seurat") {
    stop("seurat_obj_receiver should be a Seurat object")
  }
  if (class(seurat_obj_sender) != "Seurat") {
    stop("seurat_obj_sender should be a Seurat object")
  }
  if (!celltype_id_receiver %in% colnames(seurat_obj_receiver@meta.data)) {
    stop("celltype_id_receiver should be a column name in the metadata dataframe of seurat_obj_receiver")
  }
  if (!celltype_id_sender %in% colnames(seurat_obj_sender@meta.data)) {
    stop("celltype_id_sender should be a column name in the metadata dataframe of seurat_obj_sender")
  }
  if (!sample_id %in% colnames(seurat_obj_receiver@meta.data)) {
    stop("sample_id should be a column name in the metadata dataframe of seurat_obj_receiver")
  }
  if (!sample_id %in% colnames(seurat_obj_sender@meta.data)) {
    stop("sample_id should be a column name in the metadata dataframe of seurat_obj_sender")
  }
  if (!group_id %in% colnames(seurat_obj_receiver@meta.data)) {
    stop("group_id should be a column name in the metadata dataframe of seurat_obj_receiver")
  }
  if (!group_id %in% colnames(seurat_obj_sender@meta.data)) {
    stop("group_id should be a column name in the metadata dataframe of seurat_obj_sender")
  }
  if(!is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector")
  }
  if(!is.data.frame(contrast_tbl)){
    stop("contrast_tbl should be a data frame / tibble")
  }
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = seurat_obj_receiver@meta.data[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_oi_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    warning("conditions written in contrasts_oi should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  if (sum(contrasts_oi_simplified %in% unique(contrast_tbl$contrast)) != length(contrasts_oi_simplified)) {
    warning("conditions written in contrasts_oi should be in the contrast column of contrast_tbl column! This is not the case, which can lead to errors downstream.")
  }

  # Check concordance ligand-receptor network and ligand-target network - and concordance with Seurat object features
  if (!is.matrix(ligand_target_matrix)){
    stop("ligand_target_matrix should be a matrix")
  }
  if (!is.data.frame(lr_network)){
    stop("lr_network should be a data frame / tibble")
  }

  if(! "ligand" %in% colnames(lr_network)){
    lr_network = lr_network %>% dplyr::rename(ligand = from)
  }
  if(! "receptor" %in% colnames(lr_network)){
    lr_network = lr_network %>% dplyr::rename(receptor = to)
  }
  ligands_lrnetwork = lr_network$ligand %>% unique()
  receptors_lrnetwork = lr_network$receptor %>% unique()
  ligands_ligand_target_matrix = colnames(ligand_target_matrix)

  if(length(intersect(ligands_lrnetwork, ligands_ligand_target_matrix)) != length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }
  if(length(intersect(ligands_lrnetwork, ligands_ligand_target_matrix)) != length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }

  if(length(seurat_obj_sender@assays$RNA@counts %>% rownames() %>% generics::intersect(ligands_lrnetwork)) < 25 ){
    warning("Less than 25 ligands from your ligand-receptor network are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(seurat_obj_sender@assays$RNA@counts %>% rownames() %>% generics::intersect(ligands_ligand_target_matrix)) < 25 ){
    warning("Less than 25 ligands from your ligand-target matrix are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(seurat_obj_receiver@assays$RNA@counts %>% rownames() %>% generics::intersect(receptors_lrnetwork)) < 25 ){
    warning("Less than 25 receptors from your ligand-receptor network are in your expression matrix of the receiver cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }

  if(length(prioritizing_weights) != 13 | !is.double(prioritizing_weights)) {
    stop("prioritizing_weights should be a numeric vector with length 13")
  }
  names_prioritizing_weights = c("scaled_lfc_ligand",
                                 "scaled_p_val_ligand",
                                 "scaled_lfc_receptor",
                                 "scaled_p_val_receptor",
                                 "scaled_activity_scaled",
                                 "scaled_activity",
                                 "scaled_avg_exprs_ligand",
                                 "scaled_avg_frq_ligand",
                                 "scaled_avg_exprs_receptor",
                                 "scaled_avg_frq_receptor",
                                 "fraction_expressing_ligand_receptor",
                                 "scaled_abundance_sender",
                                 "scaled_abundance_receiver")
  if(sum(names_prioritizing_weights %in% names(prioritizing_weights)) != length(names_prioritizing_weights)) {
    stop("prioritizing_weights should be have the correct names. Check the vignettes and code documentation")
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
  if(!is.double(frac_cutoff)){
    stop("frac_cutoff should be numeric")
  } else {
    if(frac_cutoff <= 0 | frac_cutoff > 1) {
      warning("frac_cutoff is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.25, 0 excluded.")
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

  if(verbose == TRUE){
    print("Extract expression information from receiver")
  }

  receiver_info = suppressMessages(get_avg_frac_exprs_abund(
    seurat_obj = seurat_obj_receiver,
    sample_id = sample_id,
    celltype_id =  celltype_id_receiver,
    group_id = group_id))
  receiver_info_ic = suppressMessages(process_info_to_ic(
    info_object = receiver_info,
    ic_type = "receiver",
    lr_network = lr_network))

  if(verbose == TRUE){
    print("Extract expression information from sender")
  }

  sender_info = suppressMessages(get_avg_frac_exprs_abund(
    seurat_obj = seurat_obj_sender,
    sample_id = sample_id,
    celltype_id =  celltype_id_sender,
    group_id = group_id))
  sender_info_ic = suppressMessages(process_info_to_ic(
    info_object = sender_info,
    ic_type = "sender",
    lr_network = lr_network))

  senders_oi = Idents(seurat_obj_sender) %>% unique()
  receivers_oi = Idents(seurat_obj_receiver) %>% unique()

  sender_receiver_info = suppressMessages(combine_sender_receiver_info_ic(
    sender_info = sender_info_ic,
    receiver_info = receiver_info_ic,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network))

  ### Perform the DE analysis ----------------------------------------------------------------

  # best: pool samples that do not belong to contrasts of interest in advance!
  # necessary to have the same condition indications for both sender and receiver objects
  # so: change your condition/group metadata column!

  if(verbose == TRUE){
    print("Calculate differential expression in receiver")
  }
  receiver_de = perform_muscat_de_analysis(
    seurat_obj = seurat_obj_receiver,
    sample_id = sample_id,
    celltype_id = celltype_id_receiver,
    group_id = group_id,
    covariates = covariates,
    contrasts = contrasts_oi,
    assay_oi_sce = assay_oi_sce,
    assay_oi_pb = assay_oi_pb,
    fun_oi_pb = fun_oi_pb,
    de_method_oi = de_method_oi,
    min_cells = min_cells)

  if(verbose == TRUE){
    print("Calculate differential expression in sender")
  }
  sender_de = perform_muscat_de_analysis(
    seurat_obj = seurat_obj_sender,
    sample_id = sample_id,
    celltype_id = celltype_id_sender,
    group_id = group_id,
    covariates = covariates,
    contrasts = contrasts_oi,
    assay_oi_sce = assay_oi_sce,
    assay_oi_pb = assay_oi_pb,
    fun_oi_pb = fun_oi_pb,
    de_method_oi = de_method_oi,
    min_cells = min_cells)


  sender_receiver_de = suppressMessages(combine_sender_receiver_de(
    sender_de = sender_de,
    receiver_de = receiver_de,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network
  ))

  ### Use the DE analysis for defining the DE genes in the receiver cell type and perform NicheNet ligand activity and ligand-target inference ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Calculate NicheNet ligand activities and ligand-target links for receiver")
  }
  ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
    receiver_de = receiver_de,
    receivers_oi = receivers_oi,
    receiver_frq_df_group = receiver_info$frq_df_group,
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    frac_cutoff = frac_cutoff,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target
  )))

  rm(receiver_info_ic)
  rm(sender_info_ic)

  ### Combine the three types of information calculated above to prioritize ligand-receptor interactions ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Combine all the information in prioritization tables")
  }
  ### Remove types of information that we don't need anymore:

  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

  metadata_combined = seurat_obj_receiver@meta.data %>% dplyr::bind_rows(seurat_obj_sender@meta.data) %>% tibble::as_tibble()

  if(!is.na(covariates)){
    grouping_tbl = metadata_combined[,c(sample_id, group_id, covariates)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group",covariates)
  } else {
    grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group")
  }

  # print(grouping_tbl)

  prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    prioritizing_weights = prioritizing_weights,
    fraction_cutoff = frac_cutoff
  ))

  # Prepare Unsupervised analysis of samples! ------------------------------------------------------------------------------------------------------------
  if(verbose == TRUE){
    print("Prepare the ligand-receptor expression product matrix to be used for unsupervised analyses")
  }
  ids_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% pull(id) %>% unique()

  lr_prod_df = sender_receiver_info$avg_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sample, id, ligand_receptor_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(id, ligand_receptor_prod)
  lr_prod_mat = lr_prod_df %>% dplyr::select(-sample) %>% data.frame() %>% as.matrix()
  rownames(lr_prod_mat) = lr_prod_df$sample

  col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
  row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0)) %>% .[. == 0] %>% names()

  lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(col_remove),colnames(.) %>% generics::setdiff(col_remove)]


  return(
    list(
      sender_info = sender_info,
      receiver_info = receiver_info,
      sender_de = sender_de,
      receiver_de = receiver_de,
      sender_receiver_info = sender_receiver_info,
      sender_receiver_de =  sender_receiver_de,
      ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
      prioritization_tables = prioritization_tables,
      lr_prod_mat = lr_prod_mat,
      grouping_tbl = grouping_tbl
    )
  )
}
#' @title ms_mg_nichenet_analysis_combined
#'
#' @description \code{ms_mg_nichenet_analysis_combined}  XXXX
#' @usage ms_mg_nichenet_analysis_combined(
#' seurat_obj, celltype_id, sample_id,group_id,lr_network,ligand_target_matrix,contrasts_oi,contrast_tbl,
#' prioritizing_weights = c("scaled_lfc_ligand" = 1, "scaled_p_val_ligand" = 1, "scaled_lfc_receptor" = 1, "scaled_p_val_receptor" = 1, "scaled_activity_scaled" = 1.5,
#' "scaled_activity" = 0.5,"scaled_avg_exprs_ligand" = 1,"scaled_avg_frq_ligand" = 1,"scaled_avg_exprs_receptor" = 1, "scaled_avg_frq_receptor" = 1,
#' "fraction_expressing_ligand_receptor" = 1,"scaled_abundance_sender" = 0, "scaled_abundance_receiver" = 0),
#' assay_oi_sce = "RNA",assay_oi_pb ="counts",fun_oi_pb = "sum",de_method_oi = "edgeR",min_cells = 10,logFC_threshold = 0.25,p_val_threshold = 0.05,frac_cutoff = 0.05,p_val_adj = FALSE,top_n_target = 250, verbose = TRUE)
#'
#' @param seurat_obj XXXX
#' @param celltype_id XXXX
#' @inheritParams ms_mg_nichenet_analysis_separate
#'
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
ms_mg_nichenet_analysis_combined = function(seurat_obj,
                                            celltype_id,
                                            sample_id,
                                            group_id,
                                            lr_network,
                                            ligand_target_matrix,
                                            contrasts_oi,
                                            contrast_tbl,
                                            prioritizing_weights = c("scaled_lfc_ligand" = 1,
                                                                     "scaled_p_val_ligand" = 1,
                                                                     "scaled_lfc_receptor" = 1,
                                                                     "scaled_p_val_receptor" = 1,
                                                                     "scaled_activity_scaled" = 1.5,
                                                                     "scaled_activity" = 0.5,
                                                                     "scaled_avg_exprs_ligand" = 1,
                                                                     "scaled_avg_frq_ligand" = 1,
                                                                     "scaled_avg_exprs_receptor" = 1,
                                                                     "scaled_avg_frq_receptor" = 1,
                                                                     "fraction_expressing_ligand_receptor" = 1,
                                                                     "scaled_abundance_sender" = 0,
                                                                     "scaled_abundance_receiver" = 0),
                                            assay_oi_sce = "RNA",
                                            assay_oi_pb ="counts",
                                            fun_oi_pb = "sum",
                                            de_method_oi = "edgeR",
                                            min_cells = 10,
                                            logFC_threshold = 0.25,
                                            p_val_threshold = 0.05,
                                            frac_cutoff = 0.05,
                                            p_val_adj = FALSE,
                                            top_n_target = 250, verbose = TRUE){


  # input checks

  if (class(seurat_obj) != "Seurat") {
    stop("seurat_obj should be a Seurat object")
  }
  if (!celltype_id %in% colnames(seurat_obj@meta.data)) {
    stop("celltype_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if (!sample_id %in% colnames(seurat_obj@meta.data)) {
    stop("sample_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if (!group_id %in% colnames(seurat_obj@meta.data)) {
    stop("group_id should be a column name in the metadata dataframe of seurat_obj")
  }
  if(!is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector")
  }
  if(!is.data.frame(contrast_tbl)){
    stop("contrast_tbl should be a data frame / tibble")
  }
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = seurat_obj@meta.data[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_oi_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    warning("conditions written in contrasts_oi should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  if (sum(contrasts_oi_simplified %in% unique(contrast_tbl$contrast)) != length(contrasts_oi_simplified)) {
    warning("conditions written in contrasts_oi should be in the contrast column of contrast_tbl column! This is not the case, which can lead to errors downstream.")
  }

  # Check concordance ligand-receptor network and ligand-target network - and concordance with Seurat object features
  if (!is.matrix(ligand_target_matrix)){
    stop("ligand_target_matrix should be a matrix")
  }
  if (!is.data.frame(lr_network)){
    stop("lr_network should be a data frame / tibble")
  }

  if(! "ligand" %in% colnames(lr_network)){
    lr_network = lr_network %>% dplyr::rename(ligand = from)
  }
  if(! "receptor" %in% colnames(lr_network)){
    lr_network = lr_network %>% dplyr::rename(receptor = to)
  }
  ligands_lrnetwork = lr_network$ligand %>% unique()
  receptors_lrnetwork = lr_network$receptor %>% unique()
  ligands_ligand_target_matrix = colnames(ligand_target_matrix)

  if(length(intersect(ligands_lrnetwork, ligands_ligand_target_matrix)) != length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }
  if(length(intersect(ligands_lrnetwork, ligands_ligand_target_matrix)) != length(ligands_lrnetwork)){
    warning("Not all Ligands from your ligand-receptor network are in the ligand-target matrix")
  }

  if(length(seurat_obj@assays$RNA@counts %>% rownames() %>% generics::intersect(ligands_lrnetwork)) < 25 ){
    warning("Less than 25 ligands from your ligand-receptor network are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(seurat_obj@assays$RNA@counts %>% rownames() %>% generics::intersect(ligands_ligand_target_matrix)) < 25 ){
    warning("Less than 25 ligands from your ligand-target matrix are in your expression matrix of the sender cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }
  if(length(seurat_obj@assays$RNA@counts %>% rownames() %>% generics::intersect(receptors_lrnetwork)) < 25 ){
    warning("Less than 25 receptors from your ligand-receptor network are in your expression matrix of the receiver cell.\nDid you convert the gene symbols of the ligand-receptor network and the ligand-target matrix if your data is not from human?")
  }

  if(length(prioritizing_weights) != 13 | !is.double(prioritizing_weights)) {
    stop("prioritizing_weights should be a numeric vector with length 13")
  }
  names_prioritizing_weights = c("scaled_lfc_ligand",
                                 "scaled_p_val_ligand",
                                 "scaled_lfc_receptor",
                                 "scaled_p_val_receptor",
                                 "scaled_activity_scaled",
                                 "scaled_activity",
                                 "scaled_avg_exprs_ligand",
                                 "scaled_avg_frq_ligand",
                                 "scaled_avg_exprs_receptor",
                                 "scaled_avg_frq_receptor",
                                 "fraction_expressing_ligand_receptor",
                                 "scaled_abundance_sender",
                                 "scaled_abundance_receiver")
  if(sum(names_prioritizing_weights %in% names(prioritizing_weights)) != length(names_prioritizing_weights)) {
    stop("prioritizing_weights should be have the correct names. Check the vignettes and code documentation")
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
  if(!is.double(frac_cutoff)){
    stop("frac_cutoff should be numeric")
  } else {
    if(frac_cutoff <= 0 | frac_cutoff > 1) {
      warning("frac_cutoff is now 0 or smaller; or higher than 1. We recommend setting this parameter between 0 and 1 - preferably between 0 and 0.25, 0 excluded.")
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


  if(verbose == TRUE){
    print("Extract expression information from all cell types")
  }

  celltype_info = suppressMessages(get_avg_frac_exprs_abund(
    seurat_obj = seurat_obj,
    sample_id = sample_id,
    celltype_id =  celltype_id,
    group_id = group_id))

  receiver_info_ic = suppressMessages(process_info_to_ic(
    info_object = celltype_info,
    ic_type = "receiver",
    lr_network = lr_network))

  sender_info_ic = suppressMessages(process_info_to_ic(
    info_object = celltype_info,
    ic_type = "sender",
    lr_network = lr_network))

  senders_oi = Idents(seurat_obj) %>% unique()
  receivers_oi = Idents(seurat_obj) %>% unique()

  sender_receiver_info = suppressMessages(combine_sender_receiver_info_ic(
    sender_info = sender_info_ic,
    receiver_info = receiver_info_ic,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network))

  ### Perform the DE analysis ----------------------------------------------------------------

  # best: pool samples that do not belong to contrasts of interest in advance!
  # necessary to have the same condition indications for both sender and receiver objects
  # so: change your condition/group metadata column!

  if(verbose == TRUE){
    print("Calculate differential expression for all cell types")
  }
  celltype_de = perform_muscat_de_analysis(
    seurat_obj = seurat_obj,
    sample_id = sample_id,
    celltype_id = celltype_id,
    group_id = group_id,
    covariates = covariates,
    contrasts = contrasts_oi,
    assay_oi_sce = assay_oi_sce,
    assay_oi_pb = assay_oi_pb,
    fun_oi_pb = fun_oi_pb,
    de_method_oi = de_method_oi,
    min_cells = min_cells)

  sender_receiver_de = suppressMessages(combine_sender_receiver_de(
    sender_de = celltype_de,
    receiver_de = celltype_de,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network
  ))

  ### Use the DE analysis for defining the DE genes in the receiver cell type and perform NicheNet ligand activity and ligand-target inference ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Calculate NicheNet ligand activities and ligand-target links")
  }
  ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = receivers_oi,
    receiver_frq_df_group = celltype_info$frq_df_group,
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    frac_cutoff = frac_cutoff,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target
  )))

  rm(receiver_info_ic)
  rm(sender_info_ic)

  ### Combine the three types of information calculated above to prioritize ligand-receptor interactions ----------------------------------------------------------------
  if(verbose == TRUE){
    print("Combine all the information in prioritization tables")
  }
  ### Remove types of information that we don't need anymore:

  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

  metadata_combined = seurat_obj@meta.data %>% tibble::as_tibble()

  if(!is.na(covariates)){
    grouping_tbl = metadata_combined[,c(sample_id, group_id, covariates)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group",covariates)
  } else {
    grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group")
  }

  # print(grouping_tbl)

  prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    prioritizing_weights = prioritizing_weights,
    fraction_cutoff = frac_cutoff
  ))

  # Prepare Unsupervised analysis of samples! ------------------------------------------------------------------------------------------------------------
  if(verbose == TRUE){
    print("Prepare the ligand-receptor expression product matrix to be used for unsupervised analyses")
  }
  ids_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% pull(id) %>% unique()

  lr_prod_df = sender_receiver_info$avg_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sample, id, ligand_receptor_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(id, ligand_receptor_prod)
  lr_prod_mat = lr_prod_df %>% dplyr::select(-sample) %>% data.frame() %>% as.matrix()
  rownames(lr_prod_mat) = lr_prod_df$sample

  col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
  row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0)) %>% .[. == 0] %>% names()

  lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(col_remove),colnames(.) %>% generics::setdiff(col_remove)]


  return(
    list(
      celltype_info = celltype_info,
      celltype_de = celltype_de,
      sender_receiver_info = sender_receiver_info,
      sender_receiver_de =  sender_receiver_de,
      ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
      prioritization_tables = prioritization_tables,
      lr_prod_mat = lr_prod_mat,
      grouping_tbl = grouping_tbl
    )
  )
}
