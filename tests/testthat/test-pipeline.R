context("MultiNicheNet pipeline")
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
test_that("Pipeline for all-vs-all analysis works & plotting functions work", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output = multi_nichenet_analysis_combined(
       sce = sce,
       celltype_id = celltype_id,
       sample_id = sample_id,
       group_id = group_id,
       covariates = covariates,
       lr_network = lr_network,
       ligand_target_matrix = ligand_target_matrix,
       contrasts_oi = contrasts_oi,
       contrast_tbl = contrast_tbl)
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  output = make_lite_output(output)
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
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
  ligand_activity_target_plot = make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi, output$ligand_activities_targets_DEgenes, contrast_tbl, output$grouping_tbl, output$celltype_info, plot_legend = FALSE)
  expect_true("ggplot" %in% class(ligand_activity_target_plot$combined_plot))
  expect_true("ggplot" %in% class(ligand_activity_target_plot$legends))
  
  targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
  sample_target_plot = make_sample_target_plots(receiver_info = output$celltype_info, targets_oi, receiver_oi, output$grouping_tbl)
  expect_true("ggplot" %in% class(sample_target_plot))
  sample_target_plot_reversed = make_sample_target_plots_reversed(receiver_info = output$celltype_info, targets_oi, receiver_oi, output$grouping_tbl)
  expect_true("ggplot" %in% class(sample_target_plot_reversed))
  
  group_lfc_exprs_activity_plot = make_group_lfc_exprs_activity_plot(output$prioritization_tables, prioritized_tbl_oi, receiver_oi = receiver_oi)
  expect_true("ggplot" %in% class(group_lfc_exprs_activity_plot))
  
  target_oi = "RAB31"
  target_violin_plot = make_target_violin_plot(sce_receiver = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id_receiver = celltype_id)
  expect_true("ggplot" %in% class(target_violin_plot))
  target_feature_plot = make_target_feature_plot(sce_receiver = sce, target_oi = target_oi, group_oi = group_oi, group_id = group_id, celltype_id_receiver = celltype_id, receivers_oi = c("Malignant","myofibroblast","CAF"))
  expect_true("ggplot" %in% class(target_feature_plot))
  
  ligand_oi = "DLL1"
  receptor_oi = "NOTCH3"
  sender_oi = "Malignant"
  receiver_oi = "myofibroblast"
  
  ligand_receptor_feature_plot = make_ligand_receptor_feature_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id, senders_oi = c("Malignant","myofibroblast","CAF"), receivers_oi = c("Malignant","myofibroblast","CAF"))
  expect_true("ggplot" %in% class(ligand_receptor_feature_plot))
  
  ligand_receptor_violin_plot = make_ligand_receptor_violin_plot(sce_sender = sce, sce_receiver = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, sender_oi = sender_oi, receiver_oi = receiver_oi, sample_id = sample_id, celltype_id_sender = celltype_id, celltype_id_receiver = celltype_id)
  expect_true("ggplot" %in% class(ligand_receptor_violin_plot))
  
  prioritized_tbl_oi_prep = output$prioritization_tables$group_prioritization_tbl %>% 
    distinct(id, sender, receiver, ligand, receptor, group, prioritization_score, ligand_receptor_lfc_avg, fraction_expressing_ligand_receptor) %>% 
    filter(ligand_receptor_lfc_avg > 0 & fraction_expressing_ligand_receptor > 0) %>% top_n(30, prioritization_score) 
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))

  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id =  "Valencia",
    group_id = group_id,
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
    group_id =  "Valencia",
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  # test whether input checks are stringent enough: covariates
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = "dataset",
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_warning(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  expect_warning(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("'High-Low','High-(Low+Low+Low)/3'"),
    contrast_tbl = tibble(contrast = c("High-Low","High-(Low+Low+Low)/3"), group = c("High","High"))
  ))
  # other input checks
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    prioritizing_weights = c("Anil" = 1, "Canalla" = 2, "Fuera De Mestalla" = 3, "Lim Vete Ya" = 4, "Marcelino OE OE" = 5)
  ))
  expect_error(multi_nichenet_analysis_combined(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = "Valencia",
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
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
    celltype_id_sender = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
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
    sample_id =  "Valencia",
    group_id = group_id,
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
    group_id =  "Valencia",
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  # test whether input checks are stringent enough: covariates
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = "dataset",
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
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_warning(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  expect_warning(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("'High-Low','High-(Low+Low+Low)/3'"),
    contrast_tbl = tibble(contrast = c("High-Low","High-(Low+Low+Low)/3"), group = c("High","High"))
  ))
  # other input checks
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    prioritizing_weights = c("Anil" = 1, "Canalla" = 2, "Fuera De Mestalla" = 3, "Lim Vete Ya" = 4, "Marcelino OE OE" = 5)
  ))
  expect_error(multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  #   covariates = covariates,
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
  covariates = NA
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
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
test_that("Pipeline with wrapper function works - while correcting for batch effects", {
  sample_id = "tumor"
  group_id = "pEMT"
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = "batch"
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = multi_nichenet_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
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
  
  covariates = "batch"
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output = multi_nichenet_analysis_separate(
    sce_receiver = sce_receiver,
    sce_sender = sce_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
})
