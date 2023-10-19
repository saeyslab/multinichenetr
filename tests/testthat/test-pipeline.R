context("MultiNicheNet pipeline")
organism = "human"
if(organism == "human"){
  options(timeout=1000)
  lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% filter(ligand %in% rownames(sce)) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  gr_network = readRDS(url("https://zenodo.org/record/5884439/files/gr_network_human_21122021.rds")) %>% mutate(from = make.names(from), to = make.names(to))
  # ligand_target_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_target_matrix_nsga2r_final.rds"))
  # ligand_target_matrix = ligand_target_matrix[,intersect(rownames(sce), colnames(ligand_target_matrix))]
  ligand_target_matrix = ligand_target_matrix_test
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
  ligand_tf_matrix = readRDS(url("https://zenodo.org/record/5884439/files/ligand_tf_matrix_nsga2r_final.rds"))
  colnames(ligand_tf_matrix) = colnames(ligand_tf_matrix) %>% make.names()
  rownames(ligand_tf_matrix) = rownames(ligand_tf_matrix) %>% make.names()
  weighted_networks = readRDS(url("https://zenodo.org/record/5884439/files/weighted_networks_nsga2r_final.rds"))
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = make.names(from), to = make.names(to))
  weighted_networks$gr = weighted_networks$gr %>% mutate(from = make.names(from), to = make.names(to))
}
test_that("Pipeline for all-vs-all analysis works & plotting functions work", {
  sce = sce %>% alias_to_symbol_SCE(organism = "human") %>% makenames_SCE()
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  batches = NA
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output_naive = multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    findMarkers = TRUE)
  expect_type(output_naive,"list")
  
  # test edge-case of only 4 samples
  sce_test_cor = sce[,SummarizedExperiment::colData(sce)[,sample_id] %in% c("HN25","HN26","HN17","HN6")]
  output = multi_nichenet_analysis_combined(
    sce = sce_test_cor,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl)
  output = make_lite_output(output)
  expect_type(output,"list")
  
  output = multi_nichenet_analysis_combined(
       sce = sce,
       celltype_id = celltype_id,
       sample_id = sample_id,
       group_id = group_id,
       batches = batches,
       covariates = covariates,
       lr_network = lr_network,
       ligand_target_matrix = ligand_target_matrix,
       contrasts_oi = contrasts_oi,
       contrast_tbl = contrast_tbl)
  output$prioritization_tables$group_prioritization_tbl %>% select(id, scaled_lfc_ligand, scaled_lfc_receptor, scaled_p_val_ligand_adapted, scaled_p_val_receptor_adapted, max_scaled_activity, scaled_pb_ligand, scaled_pb_receptor, fraction_expressing_ligand_receptor,  prioritization_score )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  output = make_lite_output(output)
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  # test edge-case of no cells in one sample
  cells_remove = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c("CAF") & SummarizedExperiment::colData(sce)[,sample_id] %in% c("HN16")] %>% colnames()
  other_cells = colnames(sce) %>% setdiff(cells_remove)
  sce_test = sce[, other_cells]
  output = multi_nichenet_analysis_combined(
    sce = sce_test,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl)
  output$prioritization_tables$group_prioritization_tbl %>% select(id, scaled_lfc_ligand, scaled_lfc_receptor, scaled_p_val_ligand_adapted, scaled_p_val_receptor_adapted, max_scaled_activity, scaled_pb_ligand, scaled_pb_receptor, fraction_expressing_ligand_receptor,  prioritization_score )
  expect_type(output,"list")
  # test plotting functions
  group_oi = "High"
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score)
  lr_prod_plot = make_sample_lr_prod_plots(output$prioritization_tables, prioritized_tbl_oi)
  expect_true("ggplot" %in% class(lr_prod_plot)) 
  
  lr_prod_activity_plot = make_sample_lr_prod_activity_plots(output$prioritization_tables, prioritized_tbl_oi)
  expect_true("ggplot" %in% class(lr_prod_activity_plot))
  
  ligands_oi = output$prioritization_tables$ligand_activities_target_de_tbl %>% inner_join(contrast_tbl) %>% group_by(group, receiver) %>% distinct(ligand, receiver, group, activity) %>% top_n(5, activity) %>% pull(ligand) %>% unique()
  ligand_activity_plot = make_ligand_activity_plots(output$prioritization_tables, ligands_oi, contrast_tbl)
  expect_true("ggplot" %in% class(ligand_activity_plot))
  
  group_oi = "High"
  receiver_oi = "Malignant"
  
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% top_n(25, activity_scaled) %>% arrange(-activity_scaled)
  ligand_activity_target_plot = make_ligand_activity_target_plot(group_oi = group_oi, receiver_oi = receiver_oi, prioritized_tbl_oi = prioritized_tbl_oi, prioritization_tables = output$prioritization_tables, ligand_activities_targets_DEgenes = output$ligand_activities_targets_DEgenes, contrast_tbl, output$grouping_tbl, output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
  expect_true("ggplot" %in% class(ligand_activity_target_plot$combined_plot))
  expect_true("ggplot" %in% class(ligand_activity_target_plot$legends))
  
  targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
  sample_target_plot = make_DEgene_dotplot_pseudobulk(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, grouping_tbl = output$grouping_tbl)
  expect_true("list" %in% class(sample_target_plot))
  sample_target_plot_reversed = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, grouping_tbl = output$grouping_tbl)
  expect_true("list" %in% class(sample_target_plot_reversed))

  group_lfc_exprs_activity_plot = make_group_lfc_exprs_activity_plot(output$prioritization_tables, prioritized_tbl_oi, receiver_oi = receiver_oi)
  expect_true("ggplot" %in% class(group_lfc_exprs_activity_plot))
  
  target_oi = "PTHLH"
  target_violin_plot = make_target_violin_plot(sce_receiver = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id_receiver = celltype_id)
  expect_true("ggplot" %in% class(target_violin_plot))
  target_feature_plot = make_target_feature_plot(sce_receiver = sce, target_oi = target_oi, group_oi = group_oi, group_id = group_id, celltype_id_receiver = celltype_id, receivers_oi = c("Malignant","myofibroblast","CAF"))
  expect_true("ggplot" %in% class(target_feature_plot))
  
  ligand_oi = "APP"
  receptor_oi = "TNFRSF21"
  sender_oi = "Malignant"
  receiver_oi = "Malignant"
  
  ligand_receptor_feature_plot = make_ligand_receptor_feature_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id, senders_oi = c("Malignant","myofibroblast","CAF"), receivers_oi = c("Malignant","myofibroblast","CAF"))
  expect_true("ggplot" %in% class(ligand_receptor_feature_plot))
  
  ligand_receptor_violin_plot = make_ligand_receptor_violin_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, sender_oi = sender_oi, receiver_oi = receiver_oi, sample_id = sample_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id)
  expect_true("ggplot" %in% class(ligand_receptor_violin_plot))
  
  prioritized_tbl_oi_prep = output$prioritization_tables$group_prioritization_tbl %>% 
    distinct(id, sender, receiver, ligand, receptor, group, prioritization_score, ligand_receptor_lfc_avg, fraction_expressing_ligand_receptor) %>% 
    filter(ligand_receptor_lfc_avg > 0 & fraction_expressing_ligand_receptor > 0) %>% top_n(40, prioritization_score) 
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% 
    filter(id %in% prioritized_tbl_oi_prep$id) %>% 
    distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_prep)
  prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
  senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique())
  colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
  expect_type(circos_list,"list")
  
  group_oi = "High"
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(ligand_receptor_lfc_avg > 0 & fraction_expressing_ligand_receptor > 0 & group == group_oi) %>% top_n(25, prioritization_score)
  senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique())
  colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
  circos_list = make_circos_one_group(prioritized_tbl_oi, colors_sender, colors_receiver)
  expect_type(circos_list,"list")
  
  # the correlation plotting functions
  
  group_oi = "High"
  receiver_oi = "Malignant"
  top_n_target = 250
  lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% inner_join(output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.75 | spearman > 0.75))
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.75 | spearman < -0.75))
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
  # lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% dplyr::inner_join(output$ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(group == group_oi) %>% dplyr::ungroup() %>% dplyr::distinct(ligand, target), by = c("ligand", "target")) %>% dplyr::filter(pearson > 0.66 | spearman > 0.66)
  
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% distinct(id, ligand, receptor, sender, receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity_scaled, fraction_expressing_ligand_receptor,  prioritization_score) %>% filter(fraction_expressing_ligand_receptor > 0 & ligand_receptor_lfc_avg > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(250, prioritization_score)
  prioritized_tbl_oi = prioritized_tbl_oi %>% filter(id %in% lr_target_prior_cor_filtered$id)
  prioritized_tbl_oi = prioritized_tbl_oi %>% group_by(ligand, sender, group) %>% top_n(2, prioritization_score)
  
  lr_target_correlation_plot = make_lr_target_correlation_plot(output$prioritization_tables, prioritized_tbl_oi, lr_target_prior_cor_filtered, output$grouping_tbl, output$celltype_info, receiver_oi)
  expect_type(lr_target_correlation_plot,"list")

  graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered, prioritized_tbl_oi = prioritized_tbl_oi, colors = c("blue","red"))
  expect_type(graph_plot,"list")
  
  ligand_oi = "ADAM17"
  receptor_oi = "ITGB1"
  sender_oi = "CAF"
  receiver_oi = "Malignant"
  lr_target_scatter_plot = make_lr_target_scatter_plot(output$prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, output$celltype_info, output$grouping_tbl, lr_target_prior_cor_filtered)
  expect_true("ggplot" %in% class(lr_target_scatter_plot))
  
  lr_target_prior_cor_heatmap = make_lr_target_prior_cor_heatmap(lr_target_prior_cor_filtered)
  expect_true("ggplot" %in% class(lr_target_prior_cor_heatmap))
  
  # check the signaling pathways of these target genes
  targets_all = lr_target_prior_cor_filtered %>% filter(ligand == ligand_oi & receiver == receiver_oi & sender == sender_oi & receptor == receptor_oi)  %>% pull(target) %>% unique()
  
  active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligand_oi, receptors_all = receptor_oi, targets_all = targets_all, weighted_networks = weighted_networks)
  data_source_network = nichenetr::infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
  
  active_signaling_network_min_max = active_signaling_network
  active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  colors = c("ligand" = "darkred", "receptor" = "orange", "target" = "royalblue", "mediator" = "grey60")
  ggraph_signaling_path = suppressWarnings(make_ggraph_signaling_path(active_signaling_network_min_max, colors, ligand_oi, receptor_oi, targets_all))
  expect_type(ggraph_signaling_path,"list")
  
  # for coming calculations: reduce running time by having only one contrast of interest
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  
  # test whether input checks are stringent enough: celltype_id, sample_id, group_id
  SummarizedExperiment::colData(sce)$`Test Celltype` = SummarizedExperiment::colData(sce)$celltype 
  SummarizedExperiment::colData(sce)$`Test Group` = SummarizedExperiment::colData(sce)$pEMT 
  SummarizedExperiment::colData(sce)$`Test Sample` = SummarizedExperiment::colData(sce)$tumor 
  SummarizedExperiment::colData(sce)$false_celltype = SummarizedExperiment::colData(sce)$seurat_clusters %>% as.double()
  
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = "Test Celltype",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = "Test Sample",
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = "Test Group",
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))

  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = "batch",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id =  "batch",
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id =  "batch",
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = "false_celltype",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  # test whether input checks are stringent enough: batches
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = "dataset",
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  # test whether input checks are stringent enough: ligand-target matrix and ligand-receptor network

  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = ligand_target_matrix,
    ligand_target_matrix = lr_network,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network %>% rename(source = ligand, target = receptor),
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  # test whether input checks are stringent enough: contrast definition and contrast tbl
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("'High-Medium'"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("High-Low","Low-High"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  # expect_warning(multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   covariates = covariates,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = c("'High-Low','High-(Low+Low+Low)/3'"),
  #   contrast_tbl = tibble(contrast = c("High-Low","High-(Low+Low+Low)/3"), group = c("High","High"))
  # ))
  # other input checks

  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_sce = "Spatial"
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_pb = "Spatial"
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    de_method_oi = "EdgePython"
  ))
  # expect_warning(multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   min_cells = 0
  # ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    min_cells = "0"
  ))
  # expect_warning(multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   logFC_threshold = -0.25
  # ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    logFC_threshold = "0.25"
  ))
  # expect_warning(multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   p_val_threshold = 2
  # ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    p_val_threshold = "0.1"
  ))
  # expect_warning(multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   fraction_cutoff = 2
  # ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    p_val_adj = "maybe"
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    top_n_target = "0"
  ))
  # check whether parallel computation works
  # contrasts_oi = c("'High-Low'")
  # contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  # output = multi_nichenet_analysis_combined(
  #   sce = sce,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl, 
  #   n.cores = 4)
  # expect_type(output,"list")
  # expect_type(output$prioritization_tables,"list")
})
test_that("Pipeline for separate analysis works", {
  
  celltype_id_receiver = "celltype"
  celltype_id_sender = "celltype"
  
  sce_sender = sce[, SummarizedExperiment::colData(sce)[,celltype_id_sender] == "CAF"]
  sce_receiver = sce[, SummarizedExperiment::colData(sce)[,celltype_id_receiver] == "Malignant"]
  
  sample_id = "tumor"
  group_id = "pEMT"
  batches = NA
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  
  output_naive = multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl, 
    findMarkers = TRUE
  )
  expect_type(output_naive,"list")
  
  output = multi_nichenet_analysis_separate(
       sce_receiver = sce_receiver,
       sce_sender = sce_sender,
       celltype_id_receiver = celltype_id_receiver,
       celltype_id_sender = celltype_id_sender,
       sample_id = sample_id,
       group_id = group_id,
       batches = batches,
       covariates = covariates,
       lr_network = lr_network,
       ligand_target_matrix = ligand_target_matrix,
       contrasts_oi = contrasts_oi,
       contrast_tbl = contrast_tbl
       )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  # for coming calculations: reduce running time by having only one contrast of interest
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  
  # test whether input checks are stringent enough: celltype_id, sample_id, group_id
  SummarizedExperiment::colData(sce)$`Test Celltype` = SummarizedExperiment::colData(sce)$celltype 
  SummarizedExperiment::colData(sce)$`Test Group` = SummarizedExperiment::colData(sce)$pEMT 
  SummarizedExperiment::colData(sce)$`Test Sample` = SummarizedExperiment::colData(sce)$tumor 
  SummarizedExperiment::colData(sce)$false_celltype = SummarizedExperiment::colData(sce)$seurat_clusters %>% as.double()
  
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = "Test Celltype",
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = "Test Celltype",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = "Test Sample",
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = "Test Group",
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = "batch",
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = "batch",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id =  "batch",
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id =  "batch",
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = "false_celltype",
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_sender,
    celltype_id_sender =  "false_celltype",
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  # test whether input checks are stringent enough: batches
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = "dataset",
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  # test whether input checks are stringent enough: ligand-target matrix and ligand-receptor network

  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = ligand_target_matrix,
    ligand_target_matrix = lr_network,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network %>% rename(source = ligand, target = receptor),
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  # test whether input checks are stringent enough: contrast definition and contrast tbl
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("'High-Medium'"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("High-Low","Low-High"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  # expect_warning(multi_nichenet_analysis_separate(
  #   sce_receiver = sce_receiver,
  #   sce_sender = sce_sender,
  #   celltype_id_receiver = celltype_id_receiver,
  #   celltype_id_sender = celltype_id_sender,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   covariates = covariates,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = c("'High-Low','High-(Low+Low+Low)/3'"),
  #   contrast_tbl = tibble(contrast = c("High-Low","High-(Low+Low+Low)/3"), group = c("High","High"))
  # ))
  # other input checks

  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_sce = "Spatial"
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_pb = "Spatial"
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    de_method_oi = "EdgePython"
  ))
  # expect_warning(multi_nichenet_analysis_separate(
  #   sce_receiver = sce_receiver,
  #   sce_sender = sce_sender,
  #   celltype_id_receiver = celltype_id_receiver,
  #   celltype_id_sender = celltype_id_sender,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   min_cells = 0
  # ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    min_cells = "0"
  ))
  # expect_warning(multi_nichenet_analysis_separate(
  #   sce_receiver = sce_receiver,
  #   sce_sender = sce_sender,
  #   celltype_id_receiver = celltype_id_receiver,
  #   celltype_id_sender = celltype_id_sender,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   logFC_threshold = -0.25
  # ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    logFC_threshold = "0.25"
  ))
  # expect_warning(multi_nichenet_analysis_separate(
  #   sce_receiver = sce_receiver,
  #   sce_sender = sce_sender,
  #   celltype_id_receiver = celltype_id_receiver,
  #   celltype_id_sender = celltype_id_sender,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   p_val_threshold = 2
  # ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    p_val_threshold = "0.1"
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    fraction_cutoff = 2
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    p_val_adj = "maybe"
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    top_n_target = "0"
  ))
  
  # # check whether parallel computation works
  # contrasts_oi = c("'High-Low'")
  # contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  # output = multi_nichenet_analysis_separate(
  #   sce_receiver = sce_receiver,
  #   sce_sender = sce_sender,
  #   celltype_id_receiver = celltype_id_receiver,
  #   celltype_id_sender = celltype_id_sender,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   batches = batches,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl, 
  #   n.cores = 4)
  # expect_type(output,"list")
  # expect_type(output$prioritization_tables,"list")
  
})
test_that("Pipeline with wrapper function works", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  batches = NA
  covariates = NA
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    sender_receiver_separate = FALSE
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  top_pairs = get_top_n_lr_pairs(output$prioritization_tables, 50, groups_oi = NULL, senders_oi = NULL, receivers_oi = NULL, rank_per_group = TRUE)
  expect_type(top_pairs,"list")
  top_pairs = get_top_n_lr_pairs(output$prioritization_tables, 50, groups_oi = NULL, senders_oi = NULL, receivers_oi = NULL, rank_per_group = FALSE)
  expect_type(top_pairs,"list")
  
  DE_info = get_DE_info(
     sce = sce,
     sample_id = sample_id,
     celltype_id = celltype_id,
     group_id = group_id,
     batches = batches,
     covariates = covariates,
     contrasts = contrasts_oi)
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  comparison_plots = compare_normal_emp_pvals(DE_info, DE_info_emp)
  

})
test_that("Pipeline with wrapper function works - while correcting for batch effects", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  batches = "batch"
  covariates = NA
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    sender_receiver_separate = FALSE
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  output = output %>% make_lite_output()
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  # test batch effect indicated plots
  ligand_oi = "TNC"
  receptor_oi = "ITGB1"
  group_oi = "High"
  sender_oi = "CAF"
  receiver_oi = "Malignant"
  p_violin = make_ligand_receptor_violin_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, sender_oi = sender_oi, receiver_oi = receiver_oi, sample_id = sample_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id, batch_oi = batches)
  expect_true("ggplot" %in% class(p_violin))
  
  target_oi = "PTHLH"
  target_violin_plot = make_target_violin_plot(sce_receiver = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id_receiver = celltype_id, batch_oi = batches)
  expect_true("ggplot" %in% class(target_violin_plot))
  
  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score)
  p_sample_lr_prod_activity_batch = make_sample_lr_prod_activity_batch_plots(output$prioritization_tables, prioritized_tbl_oi, output$grouping_tbl, batch_oi = batches)
  expect_true("ggplot" %in% class(p_sample_lr_prod_activity_batch))
  
  targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
  sample_target_plot = make_DEgene_dotplot_pseudobulk_batch(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, batch_oi = batches, grouping_tbl = output$grouping_tbl)
  expect_true("list" %in% class(sample_target_plot))
  
  batches = NA
  covariates = "batch"
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    sender_receiver_separate = FALSE
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
})
test_that("Pipeline for separate analysis works - while correcting for batch effects", {
  sample_id = "tumor"
  group_id = "pEMT"

  celltype_id_receiver = "celltype"
  celltype_id_sender = "celltype"
  
  sce_sender = sce[, SummarizedExperiment::colData(sce)[,celltype_id_sender] == "CAF"]
  sce_receiver = sce[, SummarizedExperiment::colData(sce)[,celltype_id_receiver] == "Malignant"]
  
  batches = "batch"
  covariates = NA
  
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output = multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  batches = NA
  covariates = "batch"
  
  output = multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    batches = batches,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
})
test_that("Unsupervised analysis functions work", {

  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  senders_oi = "CAF"
  receivers_oi = "Malignant"
  fraction_cutoff = 0.05
  batches = NA
  lr_prod_mat = calculate_LR_pb_prod_matrix(sce, sample_id, celltype_id, group_id, senders_oi, receivers_oi, fraction_cutoff, lr_network, batches)
  expect_true("matrix" %in% class(lr_prod_mat))
  
  metadata_tbl = SummarizedExperiment::colData(sce)[, c(sample_id, group_id)]  %>% data.frame() %>% tibble::as_tibble() %>% distinct()
  colnames(metadata_tbl) = c("sample_id","group_id")
  metadata_tbl2 = metadata_tbl %>% mutate(group_id = as.double(factor(group_id)))
  output_pca = pca_LR_pb_prod_matrix(lr_prod_mat, metadata_tbl2)
  expect_type(output_pca,"list")
})
