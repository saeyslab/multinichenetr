context("Integrate scRNAseq and Visium data")
test_that("Cell type prediction of spots works", {

  seurat_obj_visium2 = predict_celltypes_of_spots(seurat_obj_visium, seurat_obj_scrnaseq)
  expect_type(seurat_obj_visium2,"S4")

  prediction_matrix = seurat_obj_visium2[["celltypes"]]@data
  seurat_obj_visium3 = add_prediction_matrix(seurat_obj_visium2, prediction_matrix, assay_name = "testCelltypes")
  expect_type(seurat_obj_visium3,"S4")
  seurat_obj_visium3 = add_prediction_matrix(seurat_obj_visium2, prediction_matrix %>% t(), assay_name = "testCelltypes2")
  expect_type(seurat_obj_visium3,"S4")

  # test whether it works with preselected features
  features_oi = c(seurat_obj_visium[["SCT"]] %>% rownames() %>% sample(1000)) %>% intersect(seurat_obj_scrnaseq[["SCT"]] %>% rownames() %>% sample(1000))

  seurat_obj_visium2 = predict_celltypes_of_spots(seurat_obj_visium, seurat_obj_scrnaseq, integration_features = features_oi)
  expect_type(seurat_obj_visium2,"S4")

  # test whether input checks are stringent enough
  expect_error(predict_celltypes_of_spots(seurat_obj_scrnaseq, seurat_obj_visium))
  expect_error(predict_celltypes_of_spots(seurat_obj_visium, seurat_obj_scrnaseq, integration_method = "ppp"))


})
test_that("Region and pseudospace prediction of cells works", {

  seurat_obj_scrnaseq2 = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  expect_type(seurat_obj_scrnaseq2,"S4")

  prediction_matrix = seurat_obj_scrnaseq2[["regions"]]@data
  seurat_obj_scrnaseq3 = add_prediction_matrix(seurat_obj_scrnaseq2, prediction_matrix, assay_name = "testRegions")
  expect_type(seurat_obj_scrnaseq3,"S4")
  seurat_obj_scrnaseq3 = add_prediction_matrix(seurat_obj_scrnaseq2, prediction_matrix %>% t(), assay_name = "testRegions")
  expect_type(seurat_obj_scrnaseq3,"S4")

  expect_error(add_prediction_matrix(seurat_obj_scrnaseq2, c("should","fail")))

  seurat_obj_scrnaseq2 = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  expect_type(seurat_obj_scrnaseq2,"S4")

  # test whether it works with preselected features
  features_oi = c(seurat_obj_visium[["SCT"]] %>% rownames() %>% sample(1000)) %>% intersect(seurat_obj_scrnaseq[["SCT"]] %>% rownames() %>% sample(1000))

  seurat_obj_scrnaseq2 = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium, integration_features = features_oi)
  expect_type(seurat_obj_scrnaseq2,"S4")

  seurat_obj_scrnaseq2 = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium, integration_features = features_oi)
  expect_type(seurat_obj_scrnaseq2,"S4")

  # test whether input checks are stringent enough
  expect_error(predict_regions_of_cells(seurat_obj_visium, seurat_obj_scrnaseq))
  expect_error(predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium, integration_method = "ppp"))

  expect_error(predict_pseudospace_of_cells(seurat_obj_visium, seurat_obj_scrnaseq))
  expect_error(predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium, integration_method = "ppp"))

})
