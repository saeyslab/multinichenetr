context("Use only Visium for cell-cell communication analysis")
test_that("Finding region-specific ligands, receptors, and ligand-receptor interactions work", {

  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  region_markers = FindAllMarkers(seurat_obj_visium, min.pct = 0.10, logfc.threshold = 0.15,return.thresh = 0.01, only.pos = TRUE, verbose = FALSE)
  region_oi = "central"
  lr_network_central = get_ligands_receptors_of_region(region_oi, seurat_obj_visium, region_markers, lr_network)

  expect_type(lr_network_central,"list")
  expect_type(lr_network_central$ligands_region,"character")
  expect_type(lr_network_central$receptors_region,"character")
  expect_type(lr_network_central$lr_network_region,"list")

  # test whether input checks are stringent enough
  expect_error(get_ligands_receptors_of_region("blabla", seurat_obj_visium, region_markers, lr_network))
  expect_error(get_ligands_receptors_of_region(region_oi, seurat_obj_visium, region_markers %>% rename(Gene = gene), lr_network))
  expect_error(get_ligands_receptors_of_region(region_oi, seurat_obj_visium, region_markers, lr_network %>% rename(ligand = from, receptor = to)))
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  expect_warning(get_ligands_receptors_of_region(region_oi, seurat_obj_visium, region_markers, lr_network))

  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  regions_all = c("central","portal")
  lr_network_central_portal = get_ligands_receptors_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network)
  expect_type(lr_network_central_portal,"list")
  expect_type(lr_network_central_portal$central,"list")
  expect_type(lr_network_central_portal$central$ligands_region,"character")
  expect_type(lr_network_central_portal$central$receptors_region,"character")
  expect_type(lr_network_central_portal$central$lr_network_region,"list")
  expect_type(lr_network_central_portal$portal,"list")
  expect_type(lr_network_central_portal$portal$ligands_region,"character")
  expect_type(lr_network_central_portal$portal$receptors_region,"character")
  expect_type(lr_network_central_portal$portal$lr_network_region,"list")

  # test whether input checks are stringent enough
  expect_error(get_ligands_receptors_of_all_regions(c("Anil", "Canalla","Fuera","De","Mestalla"), seurat_obj_visium, region_markers, lr_network))
  expect_error(get_ligands_receptors_of_all_regions(regions_all, seurat_obj_visium, region_markers %>% rename(Gene = gene), lr_network))
  expect_error(get_ligands_receptors_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network %>% rename(ligand = from, receptor = to)))
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  expect_warning(get_ligands_receptors_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network))

  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  region_markers = FindAllMarkers(seurat_obj_visium, min.pct = 0.10, logfc.threshold = 0.25,return.thresh = 0.01, only.pos = TRUE)
  regions_all = seurat_obj_visium %>% Idents() %>% levels()
  lr_network_all_regions = get_ligands_receptors_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network)

  spatial_de_lr_network = plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "SCT", merge = "product", lr_network_de_strict = TRUE)
  expect_type(spatial_de_lr_network,"list")
  expect_type(spatial_de_lr_network$raw_lr_network,"list")
  expect_type(spatial_de_lr_network$exprs_lr_network,"list")
  expect_type(spatial_de_lr_network$plots,"list")

  expect_error(plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "SCT", merge = "VCF"))

  spatial_de_lr_network = plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "SCT", merge = "product", lr_network_de_strict = FALSE)
  expect_type(spatial_de_lr_network,"list")
  expect_type(spatial_de_lr_network$raw_lr_network,"list")
  expect_type(spatial_de_lr_network$exprs_lr_network,"list")
  expect_type(spatial_de_lr_network$plots,"list")
  spatial_de_lr_network = plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "SCT", merge = "average", lr_network_de_strict = TRUE)
  expect_type(spatial_de_lr_network,"list")
  expect_type(spatial_de_lr_network$raw_lr_network,"list")
  expect_type(spatial_de_lr_network$exprs_lr_network,"list")
  expect_type(spatial_de_lr_network$plots,"list")
  spatial_de_lr_network = plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "SCT", merge = "average", lr_network_de_strict = FALSE)
  expect_type(spatial_de_lr_network,"list")
  expect_type(spatial_de_lr_network$raw_lr_network,"list")
  expect_type(spatial_de_lr_network$exprs_lr_network,"list")
  expect_type(spatial_de_lr_network$plots,"list")
  spatial_de_lr_network = plot_spatial_de_lr_network(lr_network_all_regions, seurat_obj_visium, assay_oi = "Spatial", merge = "product", lr_network_de_strict = TRUE)
  expect_type(spatial_de_lr_network,"list")
  expect_type(spatial_de_lr_network$raw_lr_network,"list")
  expect_type(spatial_de_lr_network$exprs_lr_network,"list")
  expect_type(spatial_de_lr_network$plots,"list")


})
test_that("Finding region-specific ligand activities and targets work", {

  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  region_markers = FindAllMarkers(seurat_obj_visium, min.pct = 0.10, logfc.threshold = 0.15,return.thresh = 0.01, only.pos = TRUE, verbose = FALSE)
  region_oi = "central"

  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  # ligand_target_matrix = readRDS("C:/Users/rbrowaey/work/Research/NicheNet/StaticNicheNet/paper/data_nichenet/ligand_target_matrix.rds")

  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

  ligand_activities_targets_central = get_ligand_activities_of_region(region_oi, seurat_obj_visium, region_markers, lr_network, ligand_target_matrix, target_n = 100)

  expect_type(ligand_activities_targets_central,"list")
  expect_type(ligand_activities_targets_central$ligand_activities,"list")
  expect_type(ligand_activities_targets_central$ligand_target_region,"list")
  expect_type(ligand_activities_targets_central$ligand_activities_region,"list")

  # test whether input checks are stringent enough
  expect_error(get_ligand_activities_of_region("blabla", seurat_obj_visium, region_markers, lr_network, ligand_target_matrix, target_n = 100))

  regions_all = c("central","portal")
  ligand_activities_targets_central_portal = get_ligand_activities_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network, ligand_target_matrix, target_n = 100)
  expect_type(ligand_activities_targets_central_portal,"list")
  expect_type(ligand_activities_targets_central_portal$central,"list")
  expect_type(ligand_activities_targets_central_portal$central$ligand_activities,"list")
  expect_type(ligand_activities_targets_central_portal$central$ligand_activities_region,"list")
  expect_type(ligand_activities_targets_central_portal$central$ligand_target_region,"list")
  expect_type(ligand_activities_targets_central_portal$portal,"list")
  expect_type(ligand_activities_targets_central_portal$portal$ligand_activities,"list")
  expect_type(ligand_activities_targets_central_portal$portal$ligand_activities_region,"list")
  expect_type(ligand_activities_targets_central_portal$portal$ligand_target_region,"list")

  # test whether input checks are stringent enough
  expect_error(get_ligand_activities_of_all_regions(regions_all, seurat_obj_visium, region_markers, lr_network, ligand_target_matrix, target_n = c("3",100)))


})
