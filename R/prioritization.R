scale_quantile_adapted = function(x){
  y = nichenetr::scale_quantile(x,outlier_cutoff  = 0)
  y = y + 0.001
  return(y)
}
#' @title generate_prioritization_tables
#'
#' @description \code{generate_prioritization_tables}  XXXX
#' @usage generate_prioritization_tables(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights, fraction_cutoff)
#'
#' @inheritParams ms_mg_nichenet_analysis_combined
#' @inheritParams combine_sender_receiver_info_ic
#' @param sender_receiver_info XXX
#' @param sender_receiver_de XXX
#' @param ligand_activities_targets_DEgenes XXX
#' @param sender_receiver_tbl XXX
#' @param grouping_tbl XXX
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
generate_prioritization_tables = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights, fraction_cutoff){

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
  sample_prioritization_tbl = sender_receiver_info$avg_df %>% dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(grouping_tbl) %>% dplyr::left_join(group_prioritization_tbl %>% dplyr::distinct(group, sender, receiver, ligand, receptor, prioritization_score))

  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_"))
  sample_prioritization_tbl = sample_prioritization_tbl %>% dplyr::group_by(id) %>% dplyr::mutate(scaled_LR_prod = nichenetr::scaling_zscore(ligand_receptor_prod), scaled_LR_frac = nichenetr::scaling_zscore(ligand_receptor_fraction_prod)) %>% dplyr::ungroup()

  # ligand-target information  -----------------------------------------------
  ligand_activities_target_de_tbl = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::inner_join(ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::rename(target = gene, p_val_adj = p_adj.loc)) %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled, target, ligand_target_weight, logFC, p_val, p_val_adj) %>% dplyr::distinct()

  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl))

}
