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

  # test whether input checks are stringent enough
  expect_error(ms_mg_nichenet_analysis_combined(
    seurat_obj = seurat_obj,
    celltype_id = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
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
  
  # test whether input checks are stringent enough
  expect_error(ms_mg_nichenet_analysis_separate(
    seurat_obj_receiver = seurat_obj_receiver,
    seurat_obj_sender = seurat_obj_sender,
    celltype_id_receiver = celltype_id_receiver,
    celltype_id_sender = celltype_id_sender,
    sample_id = sample_id,
    group_id = group_id,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
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
  
  # test whether input checks are stringent enough
  expect_error(ms_mg_nichenet_analysis(
    seurat_obj = seurat_obj,
    celltype_id = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    sender_receiver_separate = FALSE
  ))
  
  
})

