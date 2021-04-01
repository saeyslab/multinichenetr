scale_quantile_adapted = function(x){
  y = nichenetr::scale_quantile(x,outlier_cutoff  = 0)
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
#' @param grouping_tbl Data frame showing the groups of each sample (and covariates per sample if applicable) (columns: sample and group; and if applicable all covariates of interest)
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
#' library(Seurat)
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' 
#' metadata_abundance = seurat_obj@meta.data[,c(sample_id, group_id, celltype_id)]
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' abundance_data_receiver = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' abundance_data_sender = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#' 
#' celltype_info = get_avg_frac_exprs_abund(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' 
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = Idents(seurat_obj) %>% unique()
#' receivers_oi = Idents(seurat_obj) %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#' 
#' celltype_de = perform_muscat_de_analysis(
#'    seurat_obj = seurat_obj,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    covariates = covariates,
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
#' metadata_combined = seurat_obj@meta.data %>% tibble::as_tibble()
#' grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
#' colnames(grouping_tbl) = c("sample","group") 
#' 
#' prioritizing_weights = c("scaled_lfc_ligand" = 1,"scaled_p_val_ligand" = 1,"scaled_lfc_receptor" = 1,"scaled_p_val_receptor" = 1,"scaled_activity_scaled" = 1.5,"scaled_activity" = 0.5,"scaled_avg_exprs_ligand" = 1,"scaled_avg_frq_ligand" = 1,"scaled_avg_exprs_receptor" = 1,"scaled_avg_frq_receptor" = 1,"fraction_expressing_ligand_receptor" = 1,"scaled_abundance_sender" = 0,"scaled_abundance_receiver" = 0)
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
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(scaled_lfc_receptor = scale_quantile_adapted(lfc_receptor), scaled_p_val_receptor = scale_quantile_adapted(-p_val_receptor)) %>% dplyr::arrange(p_val_receptor )

  # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
  # receiver_ligand_activity_prioritization = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, ligand, activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled = scale_quantile_adapted(activity_scaled)) %>% dplyr::arrange(-activity_scaled )
  receiver_ligand_activity_prioritization = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled = scale_quantile_adapted(activity_scaled), scaled_activity = scale_quantile_adapted(activity)) %>% dplyr::arrange(-activity_scaled )

  # sender-focused prioritization: contrast - sender - ligand - lfc_ligand - p_adj_ligand: group by contrast and sender: score each ligand based on those rankings
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(scaled_lfc_ligand = scale_quantile_adapted(lfc_ligand), scaled_p_val_ligand = scale_quantile_adapted(-p_val_ligand)) %>% dplyr::arrange(p_val_ligand)

  # cell-type and condition specificity of expression of ligand:  per ligand: score each sender-condition combination based on expression and fraction
  ligand_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, avg_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand_group)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  ligand_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, fraction_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_frq_ligand = scale_quantile_adapted(fraction_ligand_group)) %>% dplyr::arrange(-scaled_avg_frq_ligand)

  # cell-type and condition specificity of expression of receptor:  per receptor: score each receiver-condition combination based on expression and fraction
  receptor_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, avg_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor_group)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  receptor_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, fraction_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_frq_receptor = scale_quantile_adapted(fraction_receptor_group)) %>% dplyr::arrange(-scaled_avg_frq_receptor)

  # both receptor and ligand should be expresse!
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
    dplyr::inner_join(receiver_ligand_activity_prioritization) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(ligand_receptor_expressed_prioritization) %>%
    dplyr::inner_join(sender_abundance_prioritization) %>%
    dplyr::inner_join(receiver_abundance_prioritization)

  # have a weighted average the final score (no product!!)
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    ((prioritizing_weights["scaled_lfc_ligand"] * scaled_lfc_ligand) +
                       (prioritizing_weights["scaled_p_val_ligand"] * scaled_p_val_ligand) +
                       (prioritizing_weights["scaled_lfc_receptor"] * scaled_lfc_receptor) +
                       (prioritizing_weights["scaled_p_val_receptor"] * scaled_p_val_receptor) +
                       (prioritizing_weights["scaled_activity_scaled"] * scaled_activity_scaled) +
                       (prioritizing_weights["scaled_activity"] * scaled_activity) +
                       (prioritizing_weights["scaled_avg_exprs_ligand"] * scaled_avg_exprs_ligand) +
                       (prioritizing_weights["scaled_avg_frq_ligand"] * scaled_avg_frq_ligand) +
                       (prioritizing_weights["scaled_avg_exprs_receptor"] * scaled_avg_exprs_receptor) +
                       (prioritizing_weights["scaled_avg_frq_receptor"] * scaled_avg_frq_receptor) +
                       (prioritizing_weights["fraction_expressing_ligand_receptor"] * fraction_expressing_ligand_receptor) +
                       (prioritizing_weights["scaled_abundance_sender"] * rel_abundance_scaled_sender) +
                       (prioritizing_weights["scaled_abundance_receiver"] * rel_abundance_scaled_receiver)
                    )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)


  # Sample-based Prioritization ----------------------------------------------- ----------------------------------------------------------------
  # sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score))
  sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::left_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score)) #maybe NA there if it is not considered to be expressed

  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_"))
  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::group_by(id) %>% dplyr::mutate(scaled_LR_prod = nichenetr::scaling_zscore(ligand_receptor_prod), scaled_LR_frac = nichenetr::scaling_zscore(ligand_receptor_fraction_prod)) %>% dplyr::ungroup()

  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::left_join(abundance_data_receiver) %>% dplyr::left_join(abundance_data_sender) %>% dplyr::mutate(keep_sender_receiver = keep_receiver + keep_sender)
  
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 0] = "Sender & Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 1] = "Sender or Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 2] = "Sender & Receiver present"
  sample_prioritization_tbl = sample_prioritization_tbl %>% mutate(keep_sender_receiver = factor(keep_sender_receiver, levels = c("Sender & Receiver absent",  "Sender or Receiver absent", "Sender & Receiver present")))
  
  # ligand-target information  -----------------------------------------------
  ligand_activities_target_de_tbl = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::inner_join(ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::rename(target = gene, p_val_adj = p_adj.glb)) %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled, target, ligand_target_weight, logFC, p_val, p_val_adj) %>% dplyr::distinct()

  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl))

}
