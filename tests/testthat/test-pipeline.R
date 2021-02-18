context("MultiNicheNet pipeline")
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
test_that("Pipeline for all-vs-all analysis works", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output = ms_mg_nichenet_analysis_combined(
       seurat_obj = seurat_obj,
       celltype_id = celltype_id,
       sample_id = sample_id,
       group_id = group_id,
       covariates = covariates,
       lr_network = lr_network,
       ligand_target_matrix = ligand_target_matrix,
       contrasts_oi = contrasts_oi,
       contrast_tbl = contrast_tbl       )
  expect_type(output,"list")
  expect_type(output$prioritization_tables,"list")
  
  # for coming calculations: reduce running time by having only one contrast of interest
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  
  # test whether input checks are stringent enough: celltype_id, sample_id, group_id
  seurat_obj@meta.data$`Test Celltype` = seurat_obj@meta.data$celltype 
  seurat_obj@meta.data$`Test Group` = seurat_obj@meta.data$pEMT 
  seurat_obj@meta.data$`Test Sample` = seurat_obj@meta.data$tumor 
  seurat_obj@meta.data$false_celltype = seurat_obj@meta.data$seurat_clusters %>% as.double()
  
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = "Test Celltype",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = "Test Sample",
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = "Test Group",
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))

  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id =  "Valencia",
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id =  "Valencia",
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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

  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = ligand_target_matrix,
    ligand_target_matrix = lr_network,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("'High-Medium'"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = c("High-Low","Low-High"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_warning(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  expect_warning(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  # expect_warning(ms_mg_nichenet_analysis_combined(
  #   seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  # expect_warning(ms_mg_nichenet_analysis_combined(
  #   seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  # expect_warning(ms_mg_nichenet_analysis_combined(
  #   seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  # expect_warning(ms_mg_nichenet_analysis_combined(
  #   seurat_obj = seurat_obj,
  #   celltype_id = celltype_id,
  #   sample_id = sample_id,
  #   group_id = group_id,
  #   covariates = covariates,
  #   lr_network = lr_network,
  #   ligand_target_matrix = ligand_target_matrix,
  #   contrasts_oi = contrasts_oi,
  #   contrast_tbl = contrast_tbl,
  #   frac_cutoff = 2
  # ))
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
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
})
test_that("Pipeline for separate analysis works", {
  seurat_obj_receiver = seurat_obj %>% subset(subset = celltype == "Malignant")
  seurat_obj_sender = seurat_obj %>% subset(subset = celltype == "CAF")
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id_receiver = "celltype"
  celltype_id_sender = "celltype"
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  output = ms_mg_nichenet_analysis_separate(
       seurat_obj_receiver = seurat_obj_receiver,
       seurat_obj_sender = seurat_obj_sender,
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
  seurat_obj@meta.data$`Test Celltype` = seurat_obj@meta.data$celltype 
  seurat_obj@meta.data$`Test Group` = seurat_obj@meta.data$pEMT 
  seurat_obj@meta.data$`Test Sample` = seurat_obj@meta.data$tumor 
  seurat_obj@meta.data$false_celltype = seurat_obj@meta.data$seurat_clusters %>% as.double()
  
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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

  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_warning(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_warning(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  # expect_warning(ms_mg_nichenet_analysis_separate(
  #   seurat_obj_receiver = seurat_obj_receiver,
  #   seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  # expect_warning(ms_mg_nichenet_analysis_separate(
  #   seurat_obj_receiver = seurat_obj_receiver,
  #   seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  # expect_warning(ms_mg_nichenet_analysis_separate(
  #   seurat_obj_receiver = seurat_obj_receiver,
  #   seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    frac_cutoff = 2
  ))
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
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
  
  
})
test_that("Pipeline with wrapper function works", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = NA
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))
  output = ms_mg_nichenet_analysis(
    seurat_obj = seurat_obj,
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

