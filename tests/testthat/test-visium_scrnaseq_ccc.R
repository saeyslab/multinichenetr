context("Test functions to analyze CCC after integration of Visium with scRNAseq data")
test_that("Extracting region-specific cells into Seurat objects work", {

  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)

  seurat_obj_scrnaseq_portal = get_seuratObj_region_seurat(region_oi = "portal", seurat_obj_scrnaseq, spatial_annotation_assay = "regions", label_transfer_cutoff = 0.25)
  seurat_obj_scrnaseq_central = get_seuratObj_region(region_oi = "central", seurat_obj_scrnaseq, spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25)
  seurat_obj_scrnaseq_portal_triad = get_seuratObj_region(region_oi = c("portal","portaltriad"), seurat_obj_scrnaseq, spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25)
  expect_type(seurat_obj_scrnaseq_portal,"S4")
  expect_type(seurat_obj_scrnaseq_central,"S4")
  expect_type(seurat_obj_scrnaseq_portal_triad,"S4")

  expect_error(get_seuratObj_region(region_oi = c("portal","portaltriad"), seurat_obj_scrnaseq, spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 1.25))
  expect_error(get_seuratObj_region(region_oi = c("portal","portaltriad"), seurat_obj_scrnaseq, spatial_annotation_assay = "regions", integration_method = "harmony", label_transfer_cutoff = 0.25))
  expect_error(get_seuratObj_region(region_oi = c("portal","portaltriad"), seurat_obj_scrnaseq, spatial_annotation_assay = "regionss", integration_method = "seurat", label_transfer_cutoff = 0.25))
  expect_error(get_seuratObj_region(region_oi = c("portal","portaltriad","capsule"), seurat_obj_scrnaseq, spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25))


})
test_that("Pseudospace-based Extracting region-specific cells into Seurat objects work", {

  seurat_obj_scrnaseq = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)

  seurat_obj_scrnaseq_oi = get_seuratObj_pseudospace_seurat(pseudospace_range = c(0.2,0.6), seurat_obj_scrnaseq, pseudospace_assay = "pseudospace")
  expect_type(seurat_obj_scrnaseq_oi,"S4")

  expect_error(get_seuratObj_pseudospace(pseudospace_range = c(0.2,0.1), seurat_obj_scrnaseq, pseudospace_assay = "pseudospace", integration_method = "seurat"))
  expect_error(get_seuratObj_pseudospace(pseudospace_range = c(0.2,1.1), seurat_obj_scrnaseq, pseudospace_assay = "pseudospace", integration_method = "seurat"))
  expect_error(get_seuratObj_pseudospace(pseudospace_range = c("portal","portaltriad"), seurat_obj_scrnaseq, pseudospace_assay = "pseudospace", integration_method = "harmony"))
  expect_error(get_seuratObj_pseudospace(pseudospace_range = c(0.2,0.6), seurat_obj_scrnaseq, pseudospace_assay = "pseudospaceX", integration_method = "seurat"))

})
test_that("Extracting and annotating celltypes of interest from regions of interest into Seurat objects work", {

  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  seurat_obj_scrnaseq_central_hepatocytes = get_seuratObj_celltype_region(seurat_obj_scrnaseq, celltype_oi = "hepatocyte", region_oi = "central", spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25)
  expect_type(seurat_obj_scrnaseq_central_hepatocytes,"S4")
  expect_type(seurat_obj_scrnaseq_central_hepatocytes %>% Idents(),"integer")

  expect_error(get_seuratObj_celltype_region(seurat_obj_scrnaseq,region_oi = c("portal"),  celltype_oi = "unicorns", spatial_annotation_assay = "regions", integration_method = "seurat", label_transfer_cutoff = 0.25))
  rm(seurat_obj_scrnaseq)

  seurat_obj_scrnaseq = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  seurat_obj_scrnaseq_range_oi_hepatocytes = get_seuratObj_celltype_pseudospace(seurat_obj_scrnaseq, celltype_oi = "hepatocyte", pseudospace_range = c(0.2, 0.5), pseudospace_assay = "pseudospace", integration_method = "seurat")
  expect_type(seurat_obj_scrnaseq_range_oi_hepatocytes,"S4")
  expect_type(seurat_obj_scrnaseq_range_oi_hepatocytes %>% Idents(),"integer")

  expect_error(get_seuratObj_celltype_pseudospace(seurat_obj_scrnaseq, celltype_oi = "unicorns", pseudospace_range = c(0.66, 1), pseudospace_assay = "pseudospace", integration_method = "seurat"))

})
test_that("Getting all necessary input information for the CCC analysis works + spatially-aware LR network + ligand-target links and activities", {

  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))#' head(lr_network)
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  ccc_info_portal_hepatocyte = get_ccc_info_celltype_region(region_oi = "portal", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "region", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  expect_type(ccc_info_portal_hepatocyte,"list")
  expect_type(ccc_info_portal_hepatocyte$region,"character")
  expect_type(ccc_info_portal_hepatocyte$celltype,"character")
  expect_type(ccc_info_portal_hepatocyte$n_cells,"integer")
  expect_type(ccc_info_portal_hepatocyte$cells,"character")
  expect_type(ccc_info_portal_hepatocyte$frac_cells,"double")
  expect_type(ccc_info_portal_hepatocyte$DE_table,"list")
  expect_type(ccc_info_portal_hepatocyte$de_genes,"character")
  expect_type(ccc_info_portal_hepatocyte$expressed_genes,"character")
  expect_type(ccc_info_portal_hepatocyte$ligands_de,"character")
  expect_type(ccc_info_portal_hepatocyte$receptors_de,"character")
  expect_type(ccc_info_portal_hepatocyte$ligands_expressed,"character")
  expect_type(ccc_info_portal_hepatocyte$receptors_expressed,"character")

  expect_error(get_ccc_info_celltype_region(region_oi = "portal", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "something", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10))

  rm(seurat_obj_scrnaseq)
  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  regions_oi = c("central","portal","portaltriad")
  ccc_info_hepatocyte = get_ccc_info_celltype_regions_oi(regions_oi  = regions_oi, celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "region", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  expect_type(ccc_info_hepatocyte,"list")
  expect_type(ccc_info_hepatocyte[[regions_oi[1]]],"list")
  expect_type(ccc_info_hepatocyte[[regions_oi[2]]]$region,"character")
  expect_type(ccc_info_hepatocyte[[regions_oi[3]]],"list")
  expect_type(ccc_info_hepatocyte[[regions_oi[2]]],"list")

  expect_error( get_ccc_info_celltype_regions_oi(regions_oi  = c("Mestalla",regions_oi), celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "region", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10))
  expect_error( get_ccc_info_celltype_regions_oi(regions_oi  = regions_oi, celltype_oi = "hepatocyteX", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "region", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10))

  rm(seurat_obj_scrnaseq)
  seurat_obj_scrnaseq = predict_regions_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  celltypes_oi = c("hepatocyte","cholangiocyte")
  regions_oi = c("central","portal")
  t1 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "region", spatial_annotation_assay = "regions", regions_oi = regions_oi, celltypes_oi = celltypes_oi, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  t2 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "region", spatial_annotation_assay = "regions", regions_oi = regions_oi, celltypes_oi = NULL, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  t3 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "region", spatial_annotation_assay = "regions", regions_oi = NULL, celltypes_oi = celltypes_oi, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  t4 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "region", spatial_annotation_assay = "regions", regions_oi = NULL, celltypes_oi = NULL, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10)
  expect_type(t1[[celltypes_oi[2]]],"list")
  expect_type(t2[[celltypes_oi[2]]],"list")
  expect_type(t3[[celltypes_oi[2]]],"list")
  expect_type(t4[[celltypes_oi[2]]],"list")
  spatialCCC_info = t4

  lr_network = lr_network %>% mutate(bonafide_lr = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% distinct(from, to, bonafide_lr)

  assay_oi = "SCT"
  features = union(lr_network$from, lr_network$to) %>% intersect(rownames(seurat_obj_scrnaseq))
  celltype_markers = spatialCCC_info %>% names() %>% lapply(get_markers, seurat_obj_scrnaseq, features)
  names(celltype_markers) = spatialCCC_info %>% names()
  celltype_expression = get_exprs_frac_celltype(seurat_obj_scrnaseq, features_oi = features, exprs_assay_oi = assay_oi)
  lr_network_portal = get_region_specific_lr_network(region_oi = "portal", seurat_obj_scrnaseq, spatialCCC_info, lr_network, celltype_markers, celltype_expression)#' }
  lr_network_all = get_spatial_specific_lr_network(seurat_obj_scrnaseq, spatialCCC_info, lr_network)
  expect_type(lr_network_portal,"list")
  expect_type(lr_network_all,"list")
  expect_gt(nrow(lr_network_portal), 0)
  expect_gt(nrow(lr_network_all), 0)

  # now test ligand-activities-targets for one region
  # ligand_target_matrix = readRDS("C:/Users/rbrowaey/work/Research/NicheNet/StaticNicheNet/paper/data_nichenet/ligand_target_matrix.rds")
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
  ligand_activities_targets_portal = get_region_specific_ligand_activities_targets(region_oi = "portal", seurat_obj_scrnaseq, spatialCCC_info, lr_network, ligand_target_matrix, target_n = 100)
  expect_type(ligand_activities_targets_portal,"list")

  ligand_activities_targets_all = get_spatial_specific_ligand_activities_targets(seurat_obj_scrnaseq, spatialCCC_info, lr_network, ligand_target_matrix, target_n = 100)
  expect_type(ligand_activities_targets_all,"list")


})

test_that("Getting all necessary input information for the CCC analysis works: pseudospace  + spatially-aware LR network", {

  seurat_obj_scrnaseq = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))#' head(lr_network)
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

  region2pseudospace = list("central_pseudospace" = c(0,0.3))
  ccc_info_central_pseudotime_hepatocyte = get_ccc_info_celltype_region(region_oi = "central_pseudospace", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "pseudospace", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)

  expect_type(ccc_info_central_pseudotime_hepatocyte,"list")
  expect_type(ccc_info_central_pseudotime_hepatocyte$region,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$celltype,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$n_cells,"integer")
  expect_type(ccc_info_central_pseudotime_hepatocyte$cells,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$frac_cells,"double")
  expect_type(ccc_info_central_pseudotime_hepatocyte$DE_table,"list")
  expect_type(ccc_info_central_pseudotime_hepatocyte$de_genes,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$expressed_genes,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$ligands_de,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$receptors_de,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$ligands_expressed,"character")
  expect_type(ccc_info_central_pseudotime_hepatocyte$receptors_expressed,"character")

  expect_error(get_ccc_info_celltype_region(region_oi = "central_pseudospace", celltype_oi = "hepatocyteX", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "pseudospace", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace))
  expect_error(get_ccc_info_celltype_region(region_oi = "central", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "pseudospace", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace))
  expect_error(get_ccc_info_celltype_region(region_oi = "central_pseudospace", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "region", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace))
  expect_error(get_ccc_info_celltype_region(region_oi = "central_pseudospace", celltype_oi = "hepatocyte", seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "pseudospace", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = NULL))

  rm(seurat_obj_scrnaseq)
  seurat_obj_scrnaseq = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  regions_oi = c("central_pseudospace","portal_pseudospace")
  region2pseudospace = list("central_pseudospace" = c(0,0.3), "portal_pseudospace" = c(0.5,1))
  ccc_info_hepatocyte = get_ccc_info_celltype_regions_oi(celltype_oi = "hepatocyte", regions_oi  = regions_oi, seurat_obj_scrnaseq = seurat_obj_scrnaseq, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, type = "pseudospace", exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)

  expect_type(ccc_info_hepatocyte,"list")
  expect_type(ccc_info_hepatocyte[[regions_oi[1]]],"list")
  expect_type(ccc_info_hepatocyte[[regions_oi[2]]]$region,"character")
  expect_type(ccc_info_hepatocyte[[regions_oi[2]]],"list")
  #
  #
  rm(seurat_obj_scrnaseq)
  seurat_obj_scrnaseq = predict_pseudospace_of_cells(seurat_obj_scrnaseq, seurat_obj_visium)
  celltypes_oi = c("hepatocyte","cholangiocyte")
  regions_oi = c("central_pseudospace","portal_pseudospace")
  region2pseudospace = list("central_pseudospace" = c(0,0.3), "portal_pseudospace" = c(0.5,1))

  #' spatialCCC_info_pseudospace = get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = regions_oi, celltypes_oi = celltypes_oi, lr_network = lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)

  t1 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = regions_oi, celltypes_oi = celltypes_oi, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)
  t2 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = regions_oi, celltypes_oi = NULL, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)
  t3 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = NULL, celltypes_oi = celltypes_oi, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)
  t4 = get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = NULL, celltypes_oi = NULL, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = region2pseudospace)
  expect_type(t1[[celltypes_oi[2]]],"list")
  expect_type(t2[[celltypes_oi[2]]],"list")
  expect_type(t3[[celltypes_oi[2]]],"list")
  expect_type(t4[[celltypes_oi[2]]],"list")
  spatialCCC_info = t4
  expect_error(get_spatialCCC_info(seurat_obj_scrnaseq, type = "pseudospace", spatial_annotation_assay = "pseudospace", regions_oi = NULL, celltypes_oi = NULL, lr_network, cutoff_ncells = 10, cutoff_frac_cells = 0.001, exprs_assay_oi = "SCT", exprs_frac_cutoff = 0.10, region2pseudospace = NULL))

  lr_network = lr_network %>% mutate(bonafide_lr = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% distinct(from, to, bonafide_lr)

  assay_oi = "SCT"
  features = union(lr_network$from, lr_network$to) %>% intersect(rownames(seurat_obj_scrnaseq))
  celltype_markers = spatialCCC_info %>% names() %>% lapply(get_markers, seurat_obj_scrnaseq, features)
  names(celltype_markers) = spatialCCC_info %>% names()
  celltype_expression = get_exprs_frac_celltype(seurat_obj_scrnaseq, features_oi = features, exprs_assay_oi = assay_oi)
  lr_network_portal = get_region_specific_lr_network(region_oi = "portal_pseudospace", seurat_obj_scrnaseq, spatialCCC_info, lr_network, celltype_markers, celltype_expression)#' }
  lr_network_all = get_spatial_specific_lr_network(seurat_obj_scrnaseq, spatialCCC_info, lr_network)
  expect_type(lr_network_portal,"list")
  expect_type(lr_network_all,"list")
  expect_gt(nrow(lr_network_portal), 0)
  expect_gt(nrow(lr_network_all), 0)



})





