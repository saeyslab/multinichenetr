context("Utils tests")
test_that("Average expression per region can be defined: visium", {

  avg_expression = get_avg_region_expression(genes = c("Glul","Cyp2f2"), seurat_obj_visium, assay_oi = "SCT")
  expect_type(avg_expression,"list")
  expect_type(avg_expression$expression,"double")
  expect_type(avg_expression$gene,"character")
  expect_type(avg_expression$region,"integer")
  avg_expression = get_avg_region_expression(genes = NULL, seurat_obj_visium, assay_oi = "SCT")
  expect_type(avg_expression,"list")
  expect_type(avg_expression$expression,"double")
  expect_type(avg_expression$gene,"character")
  expect_type(avg_expression$region,"integer")
  avg_expression = get_avg_region_expression(genes = c("Glul","Cyp2f2"), seurat_obj_visium, assay_oi = "Spatial")
  expect_type(avg_expression,"list")
  expect_type(avg_expression$expression,"double")
  expect_type(avg_expression$gene,"character")
  expect_type(avg_expression$region,"integer")
  # test whether input checks are stringent enough
  expect_error(get_avg_region_expression(genes = c("Glul","Cyp2f2"), seurat_obj_visium, assay_oi = "Integrated"))
  expect_error(get_avg_region_expression(genes = c("Glutamine","Glutamate"), seurat_obj_visium, assay_oi = "Integrated"))
})
test_that("Average expression per region can be defined: scRNAseq object for one celltype", {
  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  seurat_obj_scrnaseq_central_hepatocytes = get_seuratObj_celltype_region(seurat_obj_scrnaseq, celltype_oi = "hepatocyte", region_oi = "central", spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25)
  avg_expression = get_exprs_frac_region(seurat_obj_scrnaseq_central_hepatocytes, exprs_assay_oi = "SCT")

  expect_type(avg_expression,"list")
  expect_type(avg_expression$expression,"double")
  expect_type(avg_expression$gene,"character")
  expect_type(avg_expression$region_oi,"integer")

  ncells_celltype_region = seurat_obj_scrnaseq_central_hepatocytes %>% subset(ident = 1) %>% Cells() %>% length()
  ncells_celltype = seurat_obj_scrnaseq_central_hepatocytes %>% Cells() %>% length()
  frac_cells = ncells_celltype_region/ncells_celltype

  DE_table_central_hepatocytes = get_DE_celltype_region(seurat_obj_scrnaseq_central_hepatocytes, frac_cells)
  expect_type(DE_table_central_hepatocytes,"list")
  expect_type(DE_table_central_hepatocytes$p_val,"double")
  expect_type(DE_table_central_hepatocytes$avg_logFC,"double")
  expect_type(DE_table_central_hepatocytes$gene,"character")

 })
test_that("Average expression per celltype and celltype markers can be defined", {
  assay_oi = "SCT"
  features_oi = rownames(seurat_obj_scrnaseq) %>% head(1000)

  exprs_df = get_exprs_frac_celltype(seurat_obj_scrnaseq, exprs_assay_oi = assay_oi)
  markers = get_markers(celltype = "hepatocyte", seurat_obj_scrnaseq)

  celltype_markers = c("hepatocyte","cholangiocyte") %>% lapply(get_markers, seurat_obj_scrnaseq, features_oi)
  names(celltype_markers) = c("hepatocyte","cholangiocyte")
  celltype_expression = get_exprs_frac_celltype(seurat_obj_scrnaseq, features_oi = features_oi, exprs_assay_oi = assay_oi)

  expect_type(exprs_df,"list")
  expect_type(markers,"character")
  expect_type(celltype_markers,"list")
  expect_type(celltype_expression,"list")


})
