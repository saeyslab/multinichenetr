scale_quantile_adapted = function(x, outlier_cutoff = 0){
  y = nichenetr::scale_quantile(x,outlier_cutoff  = outlier_cutoff)
  y = y + 0.001
  return(y)
}
#' @title generate_prioritization_tables
#'
#' @description \code{generate_prioritization_tables}  Perform the MultiNicheNet prioritization of cell-cell interactions. 
#' User can choose the importance attached to each of the following prioritization criteria: differential expression of ligand and receptor, cell-type-condition-specificity of expression of ligand and receptor, NicheNet ligand activity, fraction of samples in a group that express a senderLigand-receiverReceptor pair, relative cell type abundance of sender/receiver.
#' @usage generate_prioritization_tables(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights, fraction_cutoff, abundance_data_receiver, abundance_data_sender)
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @inheritParams combine_sender_receiver_info_ic
#' @param sender_receiver_info Output of `combine_sender_receiver_info_ic`
#' @param sender_receiver_de Output of `combine_sender_receiver_de`
#' @param ligand_activities_targets_DEgenes Output of `get_ligand_activities_targets_DEgenes`
#' @param sender_receiver_tbl Data frame with all sender-receiver cell type combinations (columns: sender and receiver)
#' @param grouping_tbl Data frame showing the groups of each sample (and batches per sample if applicable) (columns: sample and group; and if applicable all batches of interest)
#' @param abundance_data_receiver Data frame with number of cells per cell type - sample combination;  output of `process_info_to_ic`
#' @param abundance_data_sender Data frame with number of cells per cell type - sample combination; output of `process_info_to_ic`
#' 
#' @return List containing multiple data frames prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
#'
#' @import dplyr
#' @import nichenetr
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
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' 
#' metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] 
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' abundance_data_receiver = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' abundance_data_sender = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#' 
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' 
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#' 
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    contrasts = contrasts_oi)
#'    
#' sender_receiver_de = combine_sender_receiver_de(
#'  sender_de = celltype_de,
#'  receiver_de = celltype_de,
#'  senders_oi = senders_oi,
#'  receivers_oi = receivers_oi,
#'  lr_network = lr_network)
#'  
#' ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
#'    receiver_de = celltype_de,
#'    receivers_oi = receivers_oi,
#'    receiver_frq_df_group = celltype_info$frq_df_group,
#'    ligand_target_matrix = ligand_target_matrix)
#' 
#' 
#' sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
#' metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() 
#' grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
#' colnames(grouping_tbl) = c("sample","group") 
#' 
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 1,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0)
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     prioritizing_weights = prioritizing_weights,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#' }
#'
#' @export
#'
#'
#'
#'
generate_prioritization_tables = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights, fraction_cutoff, abundance_data_receiver, abundance_data_sender){

  requireNamespace("dplyr")
  
  # Group prioritization table -------------------------------------------------------------------------------------------------------------------------------------------

  # receiver-focused prioritization for receptor: contrast - receiver - receptor - lfc_receptor - p_adj_receptor: group by contrast and receiver: score each receptor based on those rankings
  # receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  # receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = scale_quantile_adapted(lfc_receptor, outlier_cutoff = 0.001), scaled_p_val_receptor = scale_quantile_adapted(-p_val_receptor, outlier_cutoff = 0.001), scaled_lfc_pval_receptor = scale_quantile_adapted(lfc_pval_receptor, outlier_cutoff = 0.001), scaled_p_val_receptor_adapted = scale_quantile_adapted(p_val_receptor_adapted, outlier_cutoff = 0.001)) %>% dplyr::arrange(-lfc_pval_receptor)
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
  receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up= activity, activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up = scale_quantile_adapted(activity_scaled_up, outlier_cutoff = 0.01), scaled_activity_up = scale_quantile_adapted(activity_up, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up)
  receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "down") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_down = activity, activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_down = scale_quantile_adapted(activity_scaled_down, outlier_cutoff = 0.01), scaled_activity_down = scale_quantile_adapted(activity_down, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_down)
  
  # sender-focused prioritization: contrast - sender - ligand - lfc_ligand - p_adj_ligand: group by contrast and sender: score each ligand based on those rankings
  # sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  # sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = scale_quantile_adapted(lfc_ligand, outlier_cutoff = 0.001), scaled_p_val_ligand = scale_quantile_adapted(-p_val_ligand, outlier_cutoff = 0.001), scaled_lfc_pval_ligand = scale_quantile_adapted(lfc_pval_ligand, outlier_cutoff = 0.001), scaled_p_val_ligand_adapted = scale_quantile_adapted(p_val_ligand_adapted, outlier_cutoff = 0.001)) %>% dplyr::arrange(-lfc_pval_ligand)
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  
  # cell-type and condition specificity of expression of ligand:  per ligand: score each sender-condition combination based on expression and fraction
  ligand_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, avg_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand_group)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  ligand_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, fraction_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_frq_ligand = scale_quantile_adapted(fraction_ligand_group)) %>% dplyr::arrange(-scaled_avg_frq_ligand)
  ligand_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, pb_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_pb_ligand = scale_quantile_adapted(pb_ligand_group)) %>% dplyr::arrange(-scaled_pb_ligand)
  
  # cell-type and condition specificity of expression of receptor:  per receptor: score each receiver-condition combination based on expression and fraction
  receptor_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, avg_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor_group)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  receptor_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, fraction_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_frq_receptor = scale_quantile_adapted(fraction_receptor_group)) %>% dplyr::arrange(-scaled_avg_frq_receptor)
  receptor_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, pb_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_pb_receptor = scale_quantile_adapted(pb_receptor_group)) %>% dplyr::arrange(-scaled_pb_receptor)
  
  # both receptor and ligand should be expressed!
  ligand_receptor_expressed_prioritization = sender_receiver_info$frq_df %>% dplyr::inner_join(grouping_tbl)  %>% dplyr::ungroup() %>% dplyr::select(sample, group, sender, receiver, ligand, receptor, fraction_ligand, fraction_receptor ) %>% dplyr::distinct()  %>%
    dplyr::group_by(ligand, receptor, sender, receiver, group) %>%
    dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_ligand > fraction_cutoff & fraction_receptor > fraction_cutoff)) %>%
    dplyr::mutate(fraction_expressing_ligand_receptor = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_ligand_receptor) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup()

  # sender-focused prioritization of cell abundance: contrast - sender - rel abundance
  sender_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, rel_abundance_scaled_sender) %>% dplyr::distinct()  %>% dplyr::arrange(-rel_abundance_scaled_sender )

  # receiver-focused prioritization of cell abundance: contrast - receiver - rel abundance
  receiver_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, rel_abundance_scaled_receiver) %>% dplyr::distinct() %>% dplyr::arrange(-rel_abundance_scaled_receiver )


  # final group-based prioritization
  group_prioritization_tbl = contrast_tbl %>%
    dplyr::inner_join(sender_receiver_de) %>%
    dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
    dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
    dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
    dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
    dplyr::inner_join(sender_receiver_info$rel_abundance_df) %>%
    dplyr::inner_join(sender_ligand_prioritization) %>%
    dplyr::inner_join(receiver_receptor_prioritization) %>%
    dplyr::inner_join(receiver_ligand_activity_prioritization_up) %>%
    dplyr::inner_join(receiver_ligand_activity_prioritization_down) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization_pb) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization_pb) %>%
    dplyr::inner_join(ligand_receptor_expressed_prioritization) %>%
    dplyr::inner_join(sender_abundance_prioritization) %>%
    dplyr::inner_join(receiver_abundance_prioritization) %>% 
    mutate(max_scaled_activity = pmax(scaled_activity_scaled_up, scaled_activity_scaled_down), na.rm = TRUE)

  # have a weighted average the final score (no product!!)
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
                      (prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
                      (prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
                      (prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
                      (prioritizing_weights["activity_scaled"] * max_scaled_activity) +
                      # (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_up) +
                      # (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_down) +
                       # (prioritizing_weights["exprs_ligand"] * (scaled_avg_exprs_ligand + scaled_avg_frq_ligand)/2 ) +
                       # (prioritizing_weights["exprs_receptor"] * (scaled_avg_exprs_receptor + scaled_avg_frq_receptor)/2 ) +
                      (prioritizing_weights["exprs_ligand"] * scaled_pb_ligand) + # batch-effect corrected if needed
                      (prioritizing_weights["exprs_receptor"] * scaled_pb_receptor) + # batch-effect corrected if needed
                      (prioritizing_weights["frac_exprs_ligand_receptor"] * fraction_expressing_ligand_receptor) +
                      (prioritizing_weights["abund_sender"] * rel_abundance_scaled_sender) +
                      (prioritizing_weights["abund_receiver"] * rel_abundance_scaled_receiver)
                    )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)

  
  # Sample-based Prioritization ----------------------------------------------- ----------------------------------------------------------------
  # sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score))
  sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(sender_receiver_info$pb_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::left_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score)) #maybe NA there if it is not considered to be expressed

  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_"))
  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::group_by(id) %>% dplyr::mutate(scaled_LR_prod = nichenetr::scaling_zscore(ligand_receptor_prod), scaled_LR_frac = nichenetr::scaling_zscore(ligand_receptor_fraction_prod), scaled_LR_pb_prod = nichenetr::scaling_zscore(ligand_receptor_pb_prod)) %>% dplyr::ungroup()

  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::left_join(abundance_data_receiver) %>% dplyr::left_join(abundance_data_sender) 
  
  
  sample_prioritization_tbl$n_cells_sender[is.na(sample_prioritization_tbl$n_cells_sender)] = 0
  sample_prioritization_tbl$n_cells_receiver[is.na(sample_prioritization_tbl$n_cells_receiver)] = 0

  sample_prioritization_tbl$keep_sender[is.na(sample_prioritization_tbl$keep_sender)] = 0
  sample_prioritization_tbl$keep_receiver[is.na(sample_prioritization_tbl$keep_receiver)] = 0
  
  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(keep_sender_receiver = keep_receiver + keep_sender)
  
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 0] = "Sender & Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 1 & sample_prioritization_tbl$keep_receiver == 0] = "Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 1 & sample_prioritization_tbl$keep_sender == 0] = "Sender absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 2] = "Sender & Receiver present"
  sample_prioritization_tbl = sample_prioritization_tbl %>% mutate(keep_sender_receiver = factor(keep_sender_receiver, levels = c("Sender & Receiver absent",  "Receiver absent", "Sender absent", "Sender & Receiver present")))
  
  # ligand-target information  -----------------------------------------------
  ligand_activities_target_de_tbl = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::inner_join(ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::rename(target = gene, p_val_adj = p_adj)) %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled, target, ligand_target_weight, logFC, p_val, p_val_adj, direction_regulation) %>% dplyr::distinct()

  ligand_activities_target_de_tbl = ligand_activities_target_de_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
  
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(id, group, prioritization_score) %>% dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% dplyr::mutate(top_group = group) %>% dplyr::distinct(id, top_group) %>% dplyr::ungroup())
  
  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl))

}
#' @title get_top_n_lr_pairs
#'
#' @description \code{get_top_n_lr_pairs}  Get top n ligand-receptor pairs based on MultiNicheNet prioritization. Top pairs can be filtered by group, senders, receivers en based on top_n or cutoff based on the prioritization score.
#' 
#' @usage get_top_n_lr_pairs(prioritization_tables, top_n, groups_oi = NULL, senders_oi = NULL, receivers_oi = NULL, rank_per_group = TRUE)
#'
#' @param prioritization_tables output of `generate_prioritization_tables`
#' @param top_n Indicates how many top ligand-receptor pairs need to be returned
#' @param groups_oi character vector indicating the groups for which top pairs need to be returned. Default: NULL: all groups are considered.
#' @param senders_oi character vector indicating the senders for which top pairs need to be returned. Default: NULL: all senders are considered.
#' @param receivers_oi character vector indicating the receivers for which top pairs need to be returned. Default: NULL: all receivers are considered.
#' @param rank_per_group Should top_n be given per group (TRUE, default) or over all groups (FALSE)
#' 
#' @return Tibble that shows the top-ranked ligand-receptor pairs for the groups and cell types of interest
#'
#' @import dplyr
#' @import nichenetr
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
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' 
#' metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] 
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' abundance_data_receiver = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' abundance_data_sender = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#' 
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' 
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#' 
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    contrasts = contrasts_oi)
#'    
#' sender_receiver_de = combine_sender_receiver_de(
#'  sender_de = celltype_de,
#'  receiver_de = celltype_de,
#'  senders_oi = senders_oi,
#'  receivers_oi = receivers_oi,
#'  lr_network = lr_network)
#'  
#' ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
#'    receiver_de = celltype_de,
#'    receivers_oi = receivers_oi,
#'    receiver_frq_df_group = celltype_info$frq_df_group,
#'    ligand_target_matrix = ligand_target_matrix)
#' 
#' 
#' sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
#' metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() 
#' grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
#' colnames(grouping_tbl) = c("sample","group") 
#' 
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 1,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0)
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     prioritizing_weights = prioritizing_weights,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#'     
#' top50_tbl = get_top_n_lr_pairs(prioritization_tables, 50)
#' }
#'
#' @export
#'
#'
#'
#'

get_top_n_lr_pairs = function(prioritization_tables, top_n, groups_oi = NULL, senders_oi = NULL, receivers_oi = NULL, rank_per_group = TRUE){
  prioritization_tbl_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(group == top_group & fraction_expressing_ligand_receptor > 0) %>% dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score)
  if(!is.null(groups_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(group %in% groups_oi)
  }
  if(!is.null(senders_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(sender %in% senders_oi)
  }
  if(!is.null(receivers_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(receiver %in% receivers_oi)
  }
  if(rank_per_group == TRUE){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::group_by(group) %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% dplyr::filter(prioritization_rank <= top_n)
  } else {
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% dplyr::filter(prioritization_rank <= top_n)
  }
  return(prioritization_tbl_oi)
}
# generate_prioritization_tables_rankings = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights, fraction_cutoff, abundance_data_receiver, abundance_data_sender){
#   
#   requireNamespace("dplyr")
#   
#   # Group prioritization table -------------------------------------------------------------------------------------------------------------------------------------------
#   
#   # receiver-focused prioritization for receptor: contrast - receiver - receptor - lfc_receptor - p_adj_receptor: group by contrast and receiver: score each receptor based on those rankings
#   receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
#   receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
#   
#   # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
#   receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up= activity, activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up = rank(activity_scaled_up, ties.method = "average", na.last = FALSE)/max(rank(activity_scaled_up, ties.method = "average", na.last = FALSE)), scaled_activity_up = rank(activity_up, ties.method = "average", na.last = FALSE)/max(rank(activity_up, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-activity_scaled_up)
#   receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "down") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_down = activity, activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_down = rank(activity_scaled_down, ties.method = "average", na.last = FALSE)/max(rank(activity_scaled_down, ties.method = "average", na.last = FALSE)), scaled_activity_down = rank(activity_down, ties.method = "average", na.last = FALSE)/max(rank(activity_down, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-activity_scaled_down)
#   
#   # sender-focused prioritization: contrast - sender - ligand - lfc_ligand - p_adj_ligand: group by contrast and sender: score each ligand based on those rankings
#   sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
#   sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
#   
#   # cell-type and condition specificity of expression of ligand:  per ligand: score each sender-condition combination based on expression and fraction
#   ligand_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, avg_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_exprs_ligand = rank(avg_ligand_group, ties.method = "average", na.last = FALSE)/max(rank(avg_ligand_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
#   ligand_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, fraction_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_frq_ligand = rank(fraction_ligand_group, ties.method = "average", na.last = FALSE)/max(rank(fraction_ligand_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_avg_frq_ligand)
#   ligand_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, pb_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_pb_ligand = rank(pb_ligand_group, ties.method = "average", na.last = FALSE)/max(rank(pb_ligand_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_pb_ligand)
#   
#   # cell-type and condition specificity of expression of receptor:  per receptor: score each receiver-condition combination based on expression and fraction
#   receptor_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, avg_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_exprs_receptor = rank(avg_receptor_group, ties.method = "average", na.last = FALSE)/max(rank(avg_receptor_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
#   receptor_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, fraction_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_frq_receptor = rank(fraction_receptor_group, ties.method = "average", na.last = FALSE)/max(rank(fraction_receptor_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_avg_frq_receptor)
#   receptor_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, pb_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_pb_receptor = rank(pb_receptor_group, ties.method = "average", na.last = FALSE)/max(rank(pb_receptor_group, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-scaled_pb_receptor)
#   
#   # both receptor and ligand should be expressed!
#   ligand_receptor_expressed_prioritization = sender_receiver_info$frq_df %>% dplyr::inner_join(grouping_tbl)  %>% dplyr::ungroup() %>% dplyr::select(sample, group, sender, receiver, ligand, receptor, fraction_ligand, fraction_receptor ) %>% dplyr::distinct()  %>%
#     dplyr::group_by(ligand, receptor, sender, receiver, group) %>%
#     dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_ligand > fraction_cutoff & fraction_receptor > fraction_cutoff)) %>%
#     dplyr::mutate(fraction_expressing_ligand_receptor = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_ligand_receptor) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup()
#   
#   # sender-focused prioritization of cell abundance: contrast - sender - rel abundance
#   sender_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, rel_abundance_scaled_sender) %>% dplyr::distinct()  %>% dplyr::arrange(-rel_abundance_scaled_sender )
#   
#   # receiver-focused prioritization of cell abundance: contrast - receiver - rel abundance
#   receiver_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, rel_abundance_scaled_receiver) %>% dplyr::distinct() %>% dplyr::arrange(-rel_abundance_scaled_receiver )
#   
#   
#   # final group-based prioritization
#   group_prioritization_tbl = contrast_tbl %>%
#     dplyr::inner_join(sender_receiver_de) %>%
#     dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
#     dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
#     dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
#     dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
#     dplyr::inner_join(sender_receiver_info$rel_abundance_df) %>%
#     dplyr::inner_join(sender_ligand_prioritization) %>%
#     dplyr::inner_join(receiver_receptor_prioritization) %>%
#     dplyr::inner_join(receiver_ligand_activity_prioritization_up) %>%
#     dplyr::inner_join(receiver_ligand_activity_prioritization_down) %>%
#     dplyr::inner_join(ligand_celltype_specificity_prioritization) %>%
#     dplyr::inner_join(ligand_celltype_specificity_prioritization_frq) %>%
#     dplyr::inner_join(ligand_celltype_specificity_prioritization_pb) %>%
#     dplyr::inner_join(receptor_celltype_specificity_prioritization) %>%
#     dplyr::inner_join(receptor_celltype_specificity_prioritization_frq) %>%
#     dplyr::inner_join(receptor_celltype_specificity_prioritization_pb) %>%
#     dplyr::inner_join(ligand_receptor_expressed_prioritization) %>%
#     dplyr::inner_join(sender_abundance_prioritization) %>%
#     dplyr::inner_join(receiver_abundance_prioritization)
#   
#   # have a weighted average the final score (no product!!)
#   group_prioritization_tbl = group_prioritization_tbl %>%
#     dplyr::mutate(prioritization_score =
#                     (
#                       (prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
#                         (prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
#                         (prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
#                         (prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
#                         (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_up) +
#                         (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_down) +
#                         # (prioritizing_weights["exprs_ligand"] * (scaled_avg_exprs_ligand + scaled_avg_frq_ligand)/2 ) +
#                         # (prioritizing_weights["exprs_receptor"] * (scaled_avg_exprs_receptor + scaled_avg_frq_receptor)/2 ) +
#                         (prioritizing_weights["exprs_ligand"] * scaled_pb_ligand) + # batch-effect corrected if needed
#                         (prioritizing_weights["exprs_receptor"] * scaled_pb_receptor) + # batch-effect corrected if needed
#                         (prioritizing_weights["frac_exprs_ligand_receptor"] * fraction_expressing_ligand_receptor) +
#                         (prioritizing_weights["abund_sender"] * rel_abundance_scaled_sender) +
#                         (prioritizing_weights["abund_receiver"] * rel_abundance_scaled_receiver)
#                     )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)
#   
#   
#   # Sample-based Prioritization ----------------------------------------------- ----------------------------------------------------------------
#   # sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score))
#   sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(sender_receiver_info$pb_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::left_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score)) #maybe NA there if it is not considered to be expressed
#   
#   sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_"))
#   sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::group_by(id) %>% dplyr::mutate(scaled_LR_prod = nichenetr::scaling_zscore(ligand_receptor_prod), scaled_LR_frac = nichenetr::scaling_zscore(ligand_receptor_fraction_prod), scaled_LR_pb_prod = nichenetr::scaling_zscore(ligand_receptor_pb_prod)) %>% dplyr::ungroup()
#   
#   sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::left_join(abundance_data_receiver) %>% dplyr::left_join(abundance_data_sender) 
#   
#   
#   sample_prioritization_tbl$n_cells_sender[is.na(sample_prioritization_tbl$n_cells_sender)] = 0
#   sample_prioritization_tbl$n_cells_receiver[is.na(sample_prioritization_tbl$n_cells_receiver)] = 0
#   
#   sample_prioritization_tbl$keep_sender[is.na(sample_prioritization_tbl$keep_sender)] = 0
#   sample_prioritization_tbl$keep_receiver[is.na(sample_prioritization_tbl$keep_receiver)] = 0
#   
#   sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(keep_sender_receiver = keep_receiver + keep_sender)
#   
#   sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 0] = "Sender & Receiver absent"
#   sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 1 & sample_prioritization_tbl$keep_receiver == 0] = "Receiver absent"
#   sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 1 & sample_prioritization_tbl$keep_sender == 0] = "Sender absent"
#   sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 2] = "Sender & Receiver present"
#   sample_prioritization_tbl = sample_prioritization_tbl %>% mutate(keep_sender_receiver = factor(keep_sender_receiver, levels = c("Sender & Receiver absent",  "Receiver absent", "Sender absent", "Sender & Receiver present")))
#   
#   # ligand-target information  -----------------------------------------------
#   ligand_activities_target_de_tbl = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::inner_join(ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::rename(target = gene, p_val_adj = p_adj)) %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled, target, ligand_target_weight, logFC, p_val, p_val_adj, direction_regulation) %>% dplyr::distinct()
#   
#   ligand_activities_target_de_tbl = ligand_activities_target_de_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
#   group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
#   
#   group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(id, group, prioritization_score) %>% dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% dplyr::mutate(top_group = group) %>% dplyr::distinct(id, top_group) %>% dplyr::ungroup())
#   
#   return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl))
#   
# }
