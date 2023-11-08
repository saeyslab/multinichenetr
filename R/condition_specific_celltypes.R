#' @title generate_prioritization_tables_condition_specific_celltypes_sender
#'
#' @description \code{generate_prioritization_tables_condition_specific_celltypes_sender}  Perform the MultiNicheNet prioritization of cell-cell interactions. Focus on including condition-specific cell types as sender cells. This implies no DE information will be used for prioritization of ligands. 
#' Combine the following prioritization criteria in a single aggregated prioritization score: differential expression of ligand and receptor, cell-type-condition-specificity of expression of ligand and receptor, NicheNet ligand activity, fraction of samples in a group that express a senderLigand-receiverReceptor pair.
#' @usage generate_prioritization_tables_condition_specific_celltypes_sender(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, scenario = "regular", fraction_cutoff, abundance_data_receiver, abundance_data_sender, ligand_activity_down = FALSE)
#'
#' @inheritParams generate_prioritization_tables
#' 
#' @return List containing multiple data frames of prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
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
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables_condition_specific_celltypes_sender(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#' }
#'
#' @export
#'
#'
#'
#'
generate_prioritization_tables_condition_specific_celltypes_sender = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, scenario = "regular", fraction_cutoff, abundance_data_receiver, abundance_data_sender, ligand_activity_down = FALSE){
  
  requireNamespace("dplyr")
  
  if(scenario == "regular"){
    prioritizing_weights = c("de_ligand" = 0,"de_receptor" = 1,"activity_scaled" = 1,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1)
  }
  if(scenario == "lower_DE"){
    prioritizing_weights = c("de_ligand" = 0,"de_receptor" = 0.5,"activity_scaled" = 2,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1)
  }
  
  
  # Group prioritization table -------------------------------------------------------------------------------------------------------------------------------------------
  
  # receiver-focused prioritization for receptor: contrast - receiver - receptor - lfc_receptor - p_adj_receptor: group by contrast and receiver: score each receptor based on those rankings
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = base::rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = base::rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(base::rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = base::rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = base::rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(base::rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
  receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up= activity, activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up = scale_quantile_adapted(activity_scaled_up, outlier_cutoff = 0.01), scaled_activity_up = scale_quantile_adapted(activity_up, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up)
  receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "down") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_down = activity, activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_down = scale_quantile_adapted(activity_scaled_down, outlier_cutoff = 0.01), scaled_activity_down = scale_quantile_adapted(activity_down, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_down)
  
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = base::rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = base::rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(base::rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = base::rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = base::rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(base::rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  
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
  
  
  # final group-based prioritization
  if(ligand_activity_down == TRUE){
    group_prioritization_tbl = contrast_tbl %>%
      dplyr::inner_join(sender_receiver_de) %>%
      dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
      dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
      dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
      dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
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
      mutate(max_scaled_activity = pmax(scaled_activity_scaled_up, scaled_activity_scaled_down, na.rm = TRUE))
  } else {
    group_prioritization_tbl = contrast_tbl %>%
      dplyr::inner_join(sender_receiver_de) %>%
      dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
      dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
      dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
      dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
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
      mutate(max_scaled_activity = scaled_activity_scaled_up)
  }
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = prioritizing_weights["de_ligand"] + prioritizing_weights["de_receptor"] + prioritizing_weights["activity_scaled"] + prioritizing_weights["exprs_ligand"] + prioritizing_weights["exprs_receptor"] + prioritizing_weights["frac_exprs_ligand_receptor"]
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    (
                      (0.5*prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
                        (0.5*prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
                        (0.5*prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
                        (0.5*prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
                        (prioritizing_weights["activity_scaled"] * max_scaled_activity) +
                        (prioritizing_weights["exprs_ligand"] * scaled_pb_ligand) + # batch-effect corrected if needed
                        (prioritizing_weights["exprs_receptor"] * scaled_pb_receptor) + # batch-effect corrected if needed
                        (prioritizing_weights["frac_exprs_ligand_receptor"] * fraction_expressing_ligand_receptor)
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score)
  
  
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
  
  # post-process group_prioritization   -----------------------------------------------
  
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(id, group, prioritization_score) %>% dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% dplyr::mutate(top_group = group) %>% dplyr::distinct(id, top_group) %>% dplyr::ungroup())
  
  columns_to_remove = c("scaled_lfc_ligand","scaled_p_val_ligand_adapted","scaled_lfc_receptor","scaled_p_val_receptor_adapted","max_scaled_activity","scaled_pb_ligand","scaled_pb_receptor","fraction_expressing_ligand_receptor",
                        "avg_ligand_group","avg_receptor_group","ligand_receptor_prod_group", "ligand_receptor_fraction_prod_group", "sender_receiver_rel_abundance_avg", "lfc_pval_ligand", "p_val_ligand_adapted", 
                        "scaled_p_val_ligand", "scaled_lfc_pval_ligand", "lfc_pval_receptor", "p_val_receptor_adapted", "scaled_p_val_receptor", "scaled_lfc_pval_receptor","scaled_activity_scaled_up", "scaled_activity_up", "scaled_activity_scaled_down", 
                        "scaled_activity_down", "scaled_avg_exprs_ligand", "scaled_avg_frq_ligand", "scaled_avg_exprs_receptor", "scaled_avg_frq_receptor")
  columns_to_keep = group_prioritization_tbl %>% colnames() %>% setdiff(columns_to_remove)
  group_prioritization_table_source = group_prioritization_tbl %>% select(columns_to_keep) %>% distinct()
  group_prioritization_tbl = group_prioritization_tbl %>% select(contrast, group, sender, receiver, ligand, receptor, lr_interaction, id, scaled_lfc_ligand, scaled_p_val_ligand_adapted, scaled_lfc_receptor, scaled_p_val_receptor_adapted, max_scaled_activity, scaled_pb_ligand, scaled_pb_receptor, fraction_expressing_ligand_receptor, prioritization_score, top_group) %>% distinct() 
  
  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl, group_prioritization_table_source = group_prioritization_table_source))
  
}
#' @title generate_prioritization_tables_condition_specific_celltypes_receiver
#'
#' @description \code{generate_prioritization_tables_condition_specific_celltypes_receiver}  Perform the MultiNicheNet prioritization of cell-cell interactions. Focus on including condition-specific cell types as receiver cells. This implies no DE information will be used for prioritization of receptors, nor ligand activities for ligands 
#' Combine the following prioritization criteria in a single aggregated prioritization score: differential expression of ligand and receptor, cell-type-condition-specificity of expression of ligand and receptor, NicheNet ligand activity, fraction of samples in a group that express a senderLigand-receiverReceptor pair.
#' @usage generate_prioritization_tables_condition_specific_celltypes_receiver(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, scenario = "regular", fraction_cutoff, abundance_data_receiver, abundance_data_sender, ligand_activity_down = FALSE)
#'
#' @inheritParams generate_prioritization_tables
#' 
#' @return List containing multiple data frames of prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
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
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables_condition_specific_celltypes_receiver(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#' }
#'
#' @export
#'
#'
#'
#'
generate_prioritization_tables_condition_specific_celltypes_receiver = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, scenario = "regular", fraction_cutoff, abundance_data_receiver, abundance_data_sender, ligand_activity_down = FALSE){
  
  requireNamespace("dplyr")
  
  if(scenario == "regular"){
    prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 0,"activity_scaled" = 0,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1)
  }
  if(scenario == "lower_DE"){
    prioritizing_weights = c("de_ligand" = 0.5,"de_receptor" = 0,"activity_scaled" = 0,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1)
  }
  
  
  # Group prioritization table -------------------------------------------------------------------------------------------------------------------------------------------
  
  # receiver-focused prioritization for receptor: contrast - receiver - receptor - lfc_receptor - p_adj_receptor: group by contrast and receiver: score each receptor based on those rankings
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = base::rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = base::rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(base::rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = base::rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = base::rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(base::rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
  receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up= activity, activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up = scale_quantile_adapted(activity_scaled_up, outlier_cutoff = 0.01), scaled_activity_up = scale_quantile_adapted(activity_up, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up)
  receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "down") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_down = activity, activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_down = scale_quantile_adapted(activity_scaled_down, outlier_cutoff = 0.01), scaled_activity_down = scale_quantile_adapted(activity_down, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_down)
  
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = base::rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = base::rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(base::rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = base::rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(base::rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = base::rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(base::rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  
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
  
  # final group-based prioritization
  if(ligand_activity_down == TRUE){
    group_prioritization_tbl = contrast_tbl %>%
      dplyr::inner_join(sender_receiver_de) %>%
      dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
      dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
      dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
      dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
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
      mutate(max_scaled_activity = 0)
  } else {
    group_prioritization_tbl = contrast_tbl %>%
      dplyr::inner_join(sender_receiver_de) %>%
      dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
      dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
      dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
      dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
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
      mutate(max_scaled_activity = 0)
  }

  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = prioritizing_weights["de_ligand"] + prioritizing_weights["de_receptor"] + prioritizing_weights["activity_scaled"] + prioritizing_weights["exprs_ligand"] + prioritizing_weights["exprs_receptor"] + prioritizing_weights["frac_exprs_ligand_receptor"]
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    (
                      (0.5*prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
                        (0.5*prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
                        (0.5*prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
                        (0.5*prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
                        (prioritizing_weights["activity_scaled"] * max_scaled_activity) +
                        (prioritizing_weights["exprs_ligand"] * scaled_pb_ligand) + # batch-effect corrected if needed
                        (prioritizing_weights["exprs_receptor"] * scaled_pb_receptor) + # batch-effect corrected if needed
                        (prioritizing_weights["frac_exprs_ligand_receptor"] * fraction_expressing_ligand_receptor)
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score)
  
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
  
  # post-process group_prioritization   -----------------------------------------------
  
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 

  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(
    group_prioritization_tbl %>% dplyr::distinct(id, group, prioritization_score) %>% dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% dplyr::mutate(top_group = group) %>% dplyr::distinct(id, top_group) %>% dplyr::ungroup()
    )
  
  columns_to_remove = c("scaled_lfc_ligand","scaled_p_val_ligand_adapted","scaled_lfc_receptor","scaled_p_val_receptor_adapted","max_scaled_activity","scaled_pb_ligand","scaled_pb_receptor","fraction_expressing_ligand_receptor",
                        "avg_ligand_group","avg_receptor_group","ligand_receptor_prod_group", "ligand_receptor_fraction_prod_group", "sender_receiver_rel_abundance_avg", "lfc_pval_ligand", "p_val_ligand_adapted", 
                        "scaled_p_val_ligand", "scaled_lfc_pval_ligand", "lfc_pval_receptor", "p_val_receptor_adapted", "scaled_p_val_receptor", "scaled_lfc_pval_receptor","scaled_activity_scaled_up", "scaled_activity_up", "scaled_activity_scaled_down", 
                        "scaled_activity_down", "scaled_avg_exprs_ligand", "scaled_avg_frq_ligand", "scaled_avg_exprs_receptor", "scaled_avg_frq_receptor")
  columns_to_keep = group_prioritization_tbl %>% colnames() %>% setdiff(columns_to_remove)
  group_prioritization_table_source = group_prioritization_tbl %>% select(columns_to_keep) %>% distinct()
  group_prioritization_tbl = group_prioritization_tbl %>% select(contrast, group, sender, receiver, ligand, receptor, lr_interaction, id, scaled_lfc_ligand, scaled_p_val_ligand_adapted, scaled_lfc_receptor, scaled_p_val_receptor_adapted, max_scaled_activity, scaled_pb_ligand, scaled_pb_receptor, fraction_expressing_ligand_receptor, prioritization_score, top_group) %>% distinct() 
  
  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl, group_prioritization_table_source = group_prioritization_table_source))
  
}
#' @title prioritize_condition_specific_sender
#'
#' @description \code{prioritize_condition_specific_sender}  Perform the MultiNicheNet prioritization of cell-cell interactions. Focus on including condition-specific cell types as sender cells. This implies no DE information will be used for prioritization of ligands. 
#' @usage prioritize_condition_specific_sender(abundance_info, abundance_expression_info, condition_specific_celltypes, grouping_tbl, fraction_cutoff, contrast_tbl, sender_receiver_de, lr_network, ligand_activities_targets_DEgenes, scenario = "regular", ligand_activity_down = FALSE)
#'
#' @inheritParams generate_prioritization_tables
#' @param condition_specific_celltypes Character vector of condition-specific cell types
#' @param abundance_expression_info Output from `get_abundance_expression_info`
#' @param abundance_info Output from `make_abundance_plots` 
#' @inheritParams process_info_to_ic
#' 
#' @return List containing multiple data frames of prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
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
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#'     
#' condition_specific_celltypes = "CAFs"
#' abundance_info = make_abundance_plots(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id =  celltype_id, min_cells = 10, senders_oi = senders_oi, receivers_oi = receivers_oi)
#' prioritization_tables_with_condition_specific_celltype_sender = prioritize_condition_specific_sender(
#' abundance_info = abundance_info
#'abundance_expression_info = abundance_expression_info, 
#'condition_specific_celltypes = condition_specific_celltypes, 
#'grouping_tbl = grouping_tbl, 
#'fraction_cutoff = fraction_cutoff, 
#'contrast_tbl = contrast_tbl, 
#'sender_receiver_de = sender_receiver_de, 
#'lr_network = lr_network, 
#'ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'scenario = "regular",
#'ligand_activity_down = FALSE
#') 
#' 
#'    
#' }
#'
#' @export
#'
prioritize_condition_specific_sender <- function(
    abundance_info,
    abundance_expression_info, 
    condition_specific_celltypes, 
    grouping_tbl, 
    fraction_cutoff, 
    contrast_tbl, 
    sender_receiver_de, 
    lr_network, 
    ligand_activities_targets_DEgenes, 
    scenario = "regular",
    ligand_activity_down = FALSE) {
  
  requireNamespace("dplyr")
  requireNamespace("nichenetr")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  
  expessed_lr_df_condition_specific_celltypes = abundance_expression_info$sender_receiver_info$frq_df %>% filter(sender %in% condition_specific_celltypes) %>% dplyr::inner_join(grouping_tbl)  %>% dplyr::ungroup() %>% dplyr::select(sample, group, sender, ligand, fraction_ligand) %>% dplyr::distinct()  %>%
    dplyr::group_by(ligand, sender, group) %>%
    dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_ligand > fraction_cutoff)) %>%
    dplyr::mutate(fraction_expressing_ligand = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_ligand) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup() %>% filter(fraction_expressing_ligand > 0)
  
  sender_ligand_group = cross_join(abundance_info$abundance_data %>% ungroup() %>% distinct(celltype_id) %>% rename(sender = celltype_id), abundance_info$abundance_data %>% ungroup() %>% distinct(group_id) %>% rename(group = group_id)) %>% inner_join(expessed_lr_df_condition_specific_celltypes %>% distinct(sender, ligand))
  
  expessed_lr_df_condition_specific_celltypes = sender_ligand_group %>% left_join(expessed_lr_df_condition_specific_celltypes) 
  expessed_lr_df_condition_specific_celltypes %>% filter(ligand == "IFNG")
  
  sender_receiver_de_adapted = sender_receiver_de %>% distinct(contrast, receiver, receptor, lfc_receptor, p_val_receptor, p_adj_receptor) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% inner_join(contrast_tbl))  %>% ungroup() %>% select(-group) %>% select(-fraction_expressing_ligand)  %>% distinct() %>% filter(!is.na(sender)) 
  
  sender_receiver_de_adapted = sender_receiver_de_adapted %>% mutate(lfc_ligand = NA, ligand_receptor_lfc_avg = NA, p_val_ligand = NA, p_adj_ligand = NA) %>% select(colnames(sender_receiver_de))
  sender_receiver_de_adapted = bind_rows(sender_receiver_de, sender_receiver_de_adapted)
  
  sender_receiver_tbl_adapted = sender_receiver_de_adapted %>% dplyr::distinct(sender, receiver)
  
  abundance_expression_info_adapted = abundance_expression_info
  
  # pb_df_group
  pb_df_group_adapted = abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup() %>% distinct(group, receiver, receptor, pb_receptor_group) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% select(-fraction_expressing_ligand))  %>% ungroup() %>% distinct() %>% filter(!is.na(sender)) 
  abundance_expression_info_adapted$sender_receiver_info$pb_df_group = bind_rows(
    abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup(),
    abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup() %>% right_join(pb_df_group_adapted) 
  ) %>% distinct()
  
  # pb_df
  pb_df_adapted = abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup() %>% distinct(sample, receiver, receptor, pb_receptor) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% inner_join(grouping_tbl) %>% select(-fraction_expressing_ligand) %>% select(-group))  %>% ungroup() %>% distinct() %>% filter(!is.na(sender)) 
  abundance_expression_info_adapted$sender_receiver_info$pb_df = bind_rows(
    abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup(),
    abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup() %>% right_join(pb_df_adapted)
  ) %>% distinct()
  
  prioritization_tables_with_condition_specific_celltype_sender = suppressMessages(generate_prioritization_tables_condition_specific_celltypes_sender(
    sender_receiver_info = abundance_expression_info_adapted$sender_receiver_info, # had to be changed 
    sender_receiver_de = sender_receiver_de_adapted, # had to be changed
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl_adapted, # had to be changed
    grouping_tbl = grouping_tbl,
    scenario = scenario,
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down 
  ))
  
  return(prioritization_tables_with_condition_specific_celltype_sender)
  
}
#' @title prioritize_condition_specific_receiver
#'
#' @description \code{prioritize_condition_specific_receiver}  Perform the MultiNicheNet prioritization of cell-cell interactions. Focus on including condition-specific cell types as receiver cells. This implies no DE information will be used for prioritization of receptors, nor ligand activities for ligands 
#' @usage prioritize_condition_specific_receiver(abundance_info, abundance_expression_info, condition_specific_celltypes, grouping_tbl, fraction_cutoff, contrast_tbl, sender_receiver_de, lr_network, ligand_activities_targets_DEgenes, scenario = "regular", ligand_activity_down = FALSE)
#'
#' @inheritParams prioritize_condition_specific_sender
#' 
#' @return List containing multiple data frames of prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
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
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#'     
#' condition_specific_celltypes = "CAFs"
#' abundance_info = make_abundance_plots(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id =  celltype_id, min_cells = 10, senders_oi = senders_oi, receivers_oi = receivers_oi)
#' prioritization_tables_with_condition_specific_celltype_receiver = prioritize_condition_specific_receiver(
#' abundance_info = abundance_info,
#'abundance_expression_info = abundance_expression_info, 
#'condition_specific_celltypes = condition_specific_celltypes, 
#'grouping_tbl = grouping_tbl, 
#'fraction_cutoff = fraction_cutoff, 
#'contrast_tbl = contrast_tbl, 
#'sender_receiver_de = sender_receiver_de, 
#'lr_network = lr_network, 
#'ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'scenario = "regular",
#'ligand_activity_down = FALSE
#') 
#' 
#'    
#' }
#'
#' @export
#'
#'
#'
#'
prioritize_condition_specific_receiver <- function(
    abundance_info,
    abundance_expression_info, 
    condition_specific_celltypes, 
    grouping_tbl, 
    fraction_cutoff, 
    contrast_tbl, 
    sender_receiver_de, 
    lr_network, 
    ligand_activities_targets_DEgenes, 
    scenario = "regular",
    ligand_activity_down = FALSE) {
  
  requireNamespace("dplyr")
  requireNamespace("nichenetr")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  
  expessed_lr_df_condition_specific_celltypes = abundance_expression_info$sender_receiver_info$frq_df %>% filter(receiver %in% condition_specific_celltypes) %>% dplyr::inner_join(grouping_tbl)  %>% dplyr::ungroup() %>% dplyr::select(sample, group, receiver, receptor, fraction_receptor) %>% dplyr::distinct()  %>%
    dplyr::group_by(receptor, receiver, group) %>%
    dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_receptor > fraction_cutoff)) %>%
    dplyr::mutate(fraction_expressing_receptor = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_receptor) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup() %>% filter(fraction_expressing_receptor > 0)
  
  receiver_receptor_group = cross_join(abundance_info$abundance_data %>% ungroup() %>% distinct(celltype_id) %>% rename(receiver = celltype_id), abundance_info$abundance_data %>% ungroup() %>% distinct(group_id) %>% rename(group = group_id)) %>% inner_join(expessed_lr_df_condition_specific_celltypes %>% distinct(receiver, receptor))
  
  expessed_lr_df_condition_specific_celltypes = receiver_receptor_group %>% left_join(expessed_lr_df_condition_specific_celltypes) 
  
  sender_receiver_de_adapted = sender_receiver_de %>% distinct(contrast, sender, ligand, lfc_ligand, p_val_ligand, p_adj_ligand) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% inner_join(contrast_tbl))  %>% ungroup() %>% select(-group) %>% select(-fraction_expressing_receptor)  %>% distinct() %>% filter(!is.na(receiver)) 
  
  sender_receiver_de_adapted = sender_receiver_de_adapted %>% mutate(lfc_receptor = NA, ligand_receptor_lfc_avg = NA, p_val_receptor = NA, p_adj_receptor = NA) %>% select(colnames(sender_receiver_de))
  sender_receiver_de_adapted = bind_rows(sender_receiver_de, sender_receiver_de_adapted)
  
  sender_receiver_tbl_adapted = sender_receiver_de_adapted %>% dplyr::distinct(sender, receiver)
  
  abundance_expression_info_adapted = abundance_expression_info
  
  # pb_df_group
  pb_df_group_adapted = abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup() %>% distinct(group, sender, ligand, pb_ligand_group) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% select(-fraction_expressing_receptor))  %>% ungroup() %>% distinct() %>% filter(!is.na(receiver)) 
  abundance_expression_info_adapted$sender_receiver_info$pb_df_group = bind_rows(
    abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup(),
    abundance_expression_info_adapted$sender_receiver_info$pb_df_group %>% ungroup() %>% right_join(pb_df_group_adapted) 
  ) %>% distinct()
  
  # pb_df
  pb_df_adapted = abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup() %>% distinct(sample, sender, ligand, pb_ligand) %>% left_join(lr_network) %>% left_join(expessed_lr_df_condition_specific_celltypes %>% inner_join(grouping_tbl) %>% select(-fraction_expressing_receptor) %>% select(-group))  %>% ungroup() %>% distinct() %>% filter(!is.na(receiver)) 
  abundance_expression_info_adapted$sender_receiver_info$pb_df = bind_rows(
    abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup(),
    abundance_expression_info_adapted$sender_receiver_info$pb_df %>% ungroup() %>% right_join(pb_df_adapted)
  ) %>% distinct()
  
  ligand_activities_targets_DEgenes_adapted = ligand_activities_targets_DEgenes
  
  ligand_activities_targets_DEgenes_adapted$ligand_activities = bind_rows(
    ligand_activities_targets_DEgenes_adapted$ligand_activities, 
    ligand_activities_targets_DEgenes_adapted$ligand_activities %>% ungroup() %>% distinct(ligand, contrast, direction_regulation) %>% cross_join(tibble(receiver = condition_specific_celltypes)) %>% mutate(activity = NA, target = NA, ligand_target_weight = NA, activity_scaled = NA) %>% select(colnames(ligand_activities_targets_DEgenes_adapted$ligand_activities)) %>% distinct()
  )
  
  
  ligand_activities_targets_DEgenes_adapted$de_genes_df = bind_rows(
    ligand_activities_targets_DEgenes_adapted$de_genes_df, 
    cross_join(
      tibble(receiver = condition_specific_celltypes),
      tibble(contrast = contrast_tbl$contrast),
    ) %>% mutate(gene = NA, logFC = NA, p_val = NA, p_adj = NA) %>% select(colnames(ligand_activities_targets_DEgenes_adapted$de_genes_df)) %>% distinct()
  )
  
  prioritization_tables_with_condition_specific_celltype_receiver = suppressMessages(generate_prioritization_tables_condition_specific_celltypes_receiver(
    sender_receiver_info = abundance_expression_info_adapted$sender_receiver_info, # had to be changed
    sender_receiver_de = sender_receiver_de_adapted, # had to be changed
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes_adapted,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl_adapted, # had to be changed
    grouping_tbl = grouping_tbl,
    scenario = scenario,
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down 
  ))
  
  return(prioritization_tables_with_condition_specific_celltype_receiver)
  
}
