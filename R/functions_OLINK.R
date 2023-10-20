generate_prioritization_tables_OLINK = function(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, ligand_activities_targets_OLINK, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"de_ligand_OLINK" = 2, "de_receptor_OLINK" = 2, "activity_scaled" = 2, "activity_OLINK" = 1, "exprs_ligand" = 3,"exprs_receptor" = 3, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0), fraction_cutoff, abundance_data_receiver, abundance_data_sender){
  
  requireNamespace("dplyr")
  requireNamespace("nichenetr")
  
  
  olink_df_receptor = sender_receiver_de %>% distinct(receptor) %>% left_join(olink_df %>% rename(receptor = gene)) %>% mutate(contrast = "M-S")
  olink_df_reverse_receptor = olink_df_receptor %>% mutate(contrast = "S-M", logFC = -1*logFC)
  olink_df_receptor = olink_df_receptor %>% bind_rows(olink_df_reverse_receptor)
  olink_df_receptor$logFC[is.na(olink_df_receptor$logFC)] = 0
  olink_df_receptor$pval[is.na(olink_df_receptor$pval)] = 1
  
  
  olink_df_ligand = sender_receiver_de %>% distinct(ligand) %>% left_join(olink_df %>% rename(ligand = gene)) %>% mutate(contrast = "M-S")
  olink_df_reverse_ligand = olink_df_ligand %>% mutate(contrast = "S-M", logFC = -1*logFC)
  olink_df_ligand = olink_df_ligand %>% bind_rows(olink_df_reverse_ligand)
  olink_df_ligand$logFC[is.na(olink_df_ligand$logFC)] = 0
  olink_df_ligand$pval[is.na(olink_df_ligand$pval)] = 1
  # Group prioritization table -------------------------------------------------------------------------------------------------------------------------------------------
  
  # receiver-focused prioritization for receptor: contrast - receiver - receptor - lfc_receptor - p_adj_receptor: group by contrast and receiver: score each receptor based on those rankings
  # receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  # receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = nichenetr::scale_quantile_adapted(lfc_receptor, outlier_cutoff = 0.001), scaled_p_val_receptor = nichenetr::scale_quantile_adapted(-p_val_receptor, outlier_cutoff = 0.001), scaled_lfc_pval_receptor = nichenetr::scale_quantile_adapted(lfc_pval_receptor, outlier_cutoff = 0.001), scaled_p_val_receptor_adapted = nichenetr::scale_quantile_adapted(p_val_receptor_adapted, outlier_cutoff = 0.001)) %>% dplyr::arrange(-lfc_pval_receptor)
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) ) 
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  receiver_receptor_prioritization_OLINK = olink_df_receptor %>% rename(lfc_receptor = logFC, p_val_receptor = pval) %>% dplyr::ungroup() %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor)) 
  receiver_receptor_prioritization_OLINK = receiver_receptor_prioritization_OLINK %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  receiver_receptor_prioritization_OLINK = receiver_receptor_prioritization_OLINK %>% rename(scaled_lfc_receptor_OLINK = scaled_lfc_receptor, scaled_p_val_receptor_adapted_OLINK = scaled_p_val_receptor_adapted) %>% distinct(contrast, receptor, scaled_lfc_receptor_OLINK, scaled_p_val_receptor_adapted_OLINK)
  
  print(receiver_receptor_prioritization_OLINK)
  
  
  # receiver-focused prioritization for ligand: contrast - receiver - ligand - activity_scaled: group by contrast and receiver: score each ligand based on the activity
  receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up= activity, activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up = nichenetr::scale_quantile_adapted(activity_scaled_up, outlier_cutoff = 0.01), scaled_activity_up = nichenetr::scale_quantile_adapted(activity_up, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up)
  
  receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "down") %>% dplyr::select(contrast, receiver, ligand, activity, activity_scaled) %>% dplyr::rename(activity_down = activity, activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_down = nichenetr::scale_quantile_adapted(activity_scaled_down, outlier_cutoff = 0.01), scaled_activity_down = nichenetr::scale_quantile_adapted(activity_down, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_down)
  
  # sender-focused prioritization: contrast - sender - ligand - lfc_ligand - p_adj_ligand: group by contrast and sender: score each ligand based on those rankings
  # sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  # sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = nichenetr::scale_quantile_adapted(lfc_ligand, outlier_cutoff = 0.001), scaled_p_val_ligand = nichenetr::scale_quantile_adapted(-p_val_ligand, outlier_cutoff = 0.001), scaled_lfc_pval_ligand = nichenetr::scale_quantile_adapted(lfc_pval_ligand, outlier_cutoff = 0.001), scaled_p_val_ligand_adapted = nichenetr::scale_quantile_adapted(p_val_ligand_adapted, outlier_cutoff = 0.001)) %>% dplyr::arrange(-lfc_pval_ligand)
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  
  sender_ligand_prioritization_OLINK = olink_df_ligand  %>% rename(lfc_ligand = logFC, p_val_ligand = pval) %>% dplyr::ungroup() %>% dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand)) 
  sender_ligand_prioritization_OLINK = sender_ligand_prioritization_OLINK %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  sender_ligand_prioritization_OLINK = sender_ligand_prioritization_OLINK %>% rename(scaled_lfc_ligand_OLINK = scaled_lfc_ligand, scaled_p_val_ligand_adapted_OLINK = scaled_p_val_ligand_adapted) %>% distinct(contrast, ligand, scaled_lfc_ligand_OLINK, scaled_p_val_ligand_adapted_OLINK)
  
  print(sender_ligand_prioritization_OLINK)
  
  # cell-type and condition specificity of expression of ligand:  per ligand: score each sender-condition combination based on expression and fraction
  # ligand_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, avg_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_exprs_ligand = nichenetr::scale_quantile_adapted(avg_ligand_group)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  # ligand_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, fraction_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_avg_frq_ligand = nichenetr::scale_quantile_adapted(fraction_ligand_group)) %>% dplyr::arrange(-scaled_avg_frq_ligand)
  ligand_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, ligand, pb_ligand_group ) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_pb_ligand = nichenetr::scale_quantile_adapted(pb_ligand_group)) %>% dplyr::arrange(-scaled_pb_ligand)
  
  # cell-type and condition specificity of expression of receptor:  per receptor: score each receiver-condition combination based on expression and fraction
  # receptor_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, avg_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_exprs_receptor = nichenetr::scale_quantile_adapted(avg_receptor_group)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  # receptor_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, fraction_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_avg_frq_receptor = nichenetr::scale_quantile_adapted(fraction_receptor_group)) %>% dplyr::arrange(-scaled_avg_frq_receptor)
  receptor_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, pb_receptor_group ) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% dplyr::mutate(scaled_pb_receptor = nichenetr::scale_quantile_adapted(pb_receptor_group)) %>% dplyr::arrange(-scaled_pb_receptor)
  
  # both receptor and ligand should be expressed!
  ligand_receptor_expressed_prioritization = sender_receiver_info$frq_df %>% dplyr::inner_join(grouping_tbl)  %>% dplyr::ungroup() %>% dplyr::select(sample, group, sender, receiver, ligand, receptor, fraction_ligand, fraction_receptor ) %>% dplyr::distinct()  %>%
    dplyr::group_by(ligand, receptor, sender, receiver, group) %>%
    dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_ligand > fraction_cutoff & fraction_receptor > fraction_cutoff)) %>%
    dplyr::mutate(fraction_expressing_ligand_receptor = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_ligand_receptor) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup()
  
  # sender-focused prioritization of cell abundance: contrast - sender - rel abundance
  sender_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, sender, rel_abundance_scaled_sender) %>% dplyr::distinct()  %>% dplyr::arrange(-rel_abundance_scaled_sender )
  
  # receiver-focused prioritization of cell abundance: contrast - receiver - rel abundance
  receiver_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::ungroup() %>% dplyr::select(group, receiver, rel_abundance_scaled_receiver) %>% dplyr::distinct() %>% dplyr::arrange(-rel_abundance_scaled_receiver )
  
  # OLINK lignand activity prioritization
  OLINK_ligand_activity_prioritization_up = ligand_activities_targets_OLINK$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(direction_regulation == "up") %>% dplyr::select(contrast, ligand, activity, activity_scaled) %>% dplyr::rename(activity_up_OLINK = activity, activity_scaled_up_OLINK = activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(scaled_activity_scaled_up_OLINK = nichenetr::scale_quantile_adapted(activity_scaled_up_OLINK, outlier_cutoff = 0.01), scaled_activity_up_OLINK = nichenetr::scale_quantile_adapted(activity_up_OLINK, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up_OLINK)
  
  # final group-based prioritization
  group_prioritization_tbl = contrast_tbl %>%
    dplyr::inner_join(sender_receiver_de) %>%
    dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::select(-target, -ligand_target_weight) %>% dplyr::distinct()) %>%
    dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
    dplyr::inner_join(sender_receiver_info$avg_df_group) %>%
    dplyr::inner_join(sender_receiver_info$frq_df_group) %>%
    dplyr::inner_join(sender_receiver_info$rel_abundance_df) %>%
    dplyr::inner_join(sender_ligand_prioritization) %>% 
    dplyr::inner_join(sender_ligand_prioritization_OLINK) %>% 
    dplyr::inner_join(receiver_receptor_prioritization) %>%
    dplyr::inner_join(receiver_receptor_prioritization_OLINK) %>%
    dplyr::inner_join(receiver_ligand_activity_prioritization_up) %>%
    dplyr::inner_join(receiver_ligand_activity_prioritization_down) %>%
    # dplyr::inner_join(ligand_celltype_specificity_prioritization) %>%
    # dplyr::inner_join(ligand_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization_pb) %>%
    # dplyr::inner_join(receptor_celltype_specificity_prioritization) %>%
    # dplyr::inner_join(receptor_celltype_specificity_prioritization_frq) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization_pb) %>%
    dplyr::inner_join(ligand_receptor_expressed_prioritization) %>%
    dplyr::inner_join(sender_abundance_prioritization) %>%
    dplyr::inner_join(receiver_abundance_prioritization) %>% 
    dplyr::inner_join(OLINK_ligand_activity_prioritization_up) %>%
    mutate(max_scaled_activity = pmax(scaled_activity_scaled_up, scaled_activity_scaled_down), na.rm = TRUE)
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*prioritizing_weights["de_ligand"] + 2*prioritizing_weights["de_receptor"] + 2*prioritizing_weights["de_ligand_OLINK"] + 2*prioritizing_weights["de_receptor_OLINK"] + prioritizing_weights["activity_scaled"] + prioritizing_weights["activity_OLINK"] + prioritizing_weights["exprs_ligand"] + prioritizing_weights["exprs_receptor"] + prioritizing_weights["frac_exprs_ligand_receptor"] + prioritizing_weights["abund_sender"] + prioritizing_weights["abund_receiver"]
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
                        (prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
                        (prioritizing_weights["de_ligand_OLINK"] * scaled_lfc_ligand_OLINK) +
                        (prioritizing_weights["de_receptor_OLINK"] * scaled_lfc_receptor_OLINK) +
                        (prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
                        (prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
                        (prioritizing_weights["de_ligand_OLINK"] * scaled_p_val_ligand_adapted_OLINK) +
                        # (prioritizing_weights["de_receptor_OLINK"] * scaled_p_val_receptor_adapted_OLINK) +
                        (prioritizing_weights["activity_scaled"] * max_scaled_activity) +
                        (prioritizing_weights["activity_OLINK"] * scaled_activity_scaled_up_OLINK) +
                        # (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_up) +
                        # (prioritizing_weights["activity_scaled"] * scaled_activity_scaled_down) +
                        # (prioritizing_weights["exprs_ligand"] * (scaled_avg_exprs_ligand + scaled_avg_frq_ligand)/2 ) +
                        # (prioritizing_weights["exprs_receptor"] * (scaled_avg_exprs_receptor + scaled_avg_frq_receptor)/2 ) +
                        (prioritizing_weights["exprs_ligand"] * scaled_pb_ligand) + # batch-effect corrected if needed
                        (prioritizing_weights["exprs_receptor"] * scaled_pb_receptor) + # batch-effect corrected if needed
                        (prioritizing_weights["frac_exprs_ligand_receptor"] * fraction_expressing_ligand_receptor) +
                        (prioritizing_weights["abund_sender"] * rel_abundance_scaled_sender) +
                        (prioritizing_weights["abund_receiver"] * rel_abundance_scaled_receiver)
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
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down"))) 
  
  ligand_activities_target_de_tbl_OLINK = ligand_activities_targets_OLINK$ligand_activities %>% dplyr::inner_join(ligand_activities_targets_OLINK$de_genes_df %>% dplyr::rename(target = gene)) %>% dplyr::select(contrast, ligand, activity, activity_scaled, target, ligand_target_weight, logFC, pval, direction_regulation) %>% dplyr::distinct()
  
  
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(group_prioritization_tbl %>% dplyr::distinct(id, group, prioritization_score) %>% dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% dplyr::mutate(top_group = group) %>% dplyr::distinct(id, top_group) %>% dplyr::ungroup())
  
  return(list(group_prioritization_tbl = group_prioritization_tbl, sample_prioritization_tbl = sample_prioritization_tbl, ligand_activities_target_de_tbl = ligand_activities_target_de_tbl, ligand_activities_target_de_tbl_OLINK = ligand_activities_target_de_tbl_OLINK))
  
}

get_ligand_activities_targets_OLINK = function(olink_df,  ligand_target_matrix, logFC_threshold = 1, pval_threshold = 0.05, top_n_target = 250, verbose = FALSE){
  
  requireNamespace("dplyr")

  olink_df = olink_df %>% mutate(contrast = contrast_tbl$contrast %>% .[1])
  olink_df_reverse = olink_df %>% mutate(contrast = contrast_tbl$contrast %>% .[2], logFC = -1*logFC)
  de_output_tidy = olink_df %>% bind_rows(olink_df_reverse)

  background_expressed_genes = de_output_tidy$gene %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
  ligands = colnames(ligand_target_matrix)
      
  geneset_vs_ligand_activities = list()
  ligand_activities_targets_geneset = list()
  for(i in seq(length(de_output_tidy$contrast %>% unique()))){
    contrast_oi = de_output_tidy$contrast %>% unique() %>% .[i]
    de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC >= logFC_threshold & pval <= pval_threshold)
    geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
    if(verbose == TRUE){
      print("contrast_oi:")
      print(contrast_oi)
    }
    if(length(geneset_oi) > 0){
      if(verbose == TRUE){
        print("Number of upregulated DE genes (gene set of interest): ")
        print(length(geneset_oi))
      }    
      geneset_id = geneset_oi %>% paste(collapse = ".")  
      if(geneset_id %in% names(geneset_vs_ligand_activities)) {
        ligand_activities = geneset_vs_ligand_activities[[geneset_id]]$ligand_activities_df
      } else {
        ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
        geneset_vs_ligand_activities[[geneset_id]] = list(ligand_activities_df = ligand_activities)
      } 
      ligand_activities = ligand_activities %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = aupr_corrected) %>% dplyr::select(-pearson, -auroc, -aupr)
      ligand_target_df = ligand_activities$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
      ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df)  %>% dplyr::group_by(contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity), direction_regulation = "up")
    } else {
      warning(paste0("In condition ", contrast_oi, " there seem to be no upregulated DE genes - so ligand activities will be NA. Please check the DE output."))
      ligand_activities = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, ligand_target_weight = NA, activity_scaled = NA, direction_regulation = "up")
    }
    de_genes_df = de_tbl_geneset %>% dplyr::mutate(contrast = contrast_oi) 
    ligand_activities_targets_geneset[[i]] = list(ligand_activities = ligand_activities, de_genes_df = de_genes_df)
  }
      
  ligand_activities = ligand_activities_targets_geneset %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()
  de_genes_df = ligand_activities_targets_geneset %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()    
      
  return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
}



