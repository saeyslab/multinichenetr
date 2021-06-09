#' @title get_ligand_activities_targets_DEgenes
#'
#' @description \code{get_ligand_activities_targets_DEgenes}  Predict NicheNet ligand activities and ligand-target links for the receiver cell types. Uses `nichenetr::predict_ligand_activities()` and `nichenetr::get_weighted_ligand_target_links()` under the hood.
#' @usage get_ligand_activities_targets_DEgenes(receiver_de, receivers_oi,  ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_target = 250, verbose = FALSE, n.cores = 1)
#'
#' @inheritParams multi_nichenet_analysis_separate
#' @inheritParams combine_sender_receiver_info_ic
#' @inheritParams combine_sender_receiver_de
#'
#' @return List with two data frames: one data frame containing the ligand activities and ligand-target links, one containing the DE gene information.
#'
#' @import dplyr
#' @import muscat
#' @import nichenetr
#' @import tibble
#' @import tidyr
#' @importFrom purrr map
#' @importFrom parallel makeCluster parLapply clusterExport stopCluster clusterEvalQ
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' receivers_oi = Idents(seurat_obj) %>% unique()
#' celltype_info = get_avg_frac_exprs_abund(seurat_obj = seurat_obj, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' celltype_de = perform_muscat_de_analysis(
#'    seurat_obj = seurat_obj,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#' receiver_de = celltype_de$de_output_tidy
#' ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
#'    receiver_de = receiver_de,
#'    receivers_oi = receivers_oi,
#'    ligand_target_matrix = ligand_target_matrix)
#' }
#'
#' @export
#'
get_ligand_activities_targets_DEgenes = function(receiver_de, receivers_oi,  ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_target = 250, verbose = FALSE, n.cores = 1){
  
  requireNamespace("dplyr")
  
  n.cores_oi = min(n.cores, length(receivers_oi))
  
  if(n.cores_oi > 1){
    if( Sys.info()[['sysname']] == "Windows"){
      clust = parallel::makeCluster(n.cores_oi)
      parallel::clusterExport(clust, c("receivers_oi","receiver_de","ligand_target_matrix","verbose","logFC_threshold","p_val_threshold", "p_val_adj","top_n_target"), envir = environment())
      parallel::clusterEvalQ(clust, library(dplyr))
      parallel::clusterEvalQ(clust, library(nichenetr))
      parallel::clusterEvalQ(clust, library(muscat))
      parallel::clusterEvalQ(clust, library(tidyr))
      parallel::clusterEvalQ(clust, library(tibble))
      parallel::clusterEvalQ(clust, library(purrr))
    } else {
      clust = parallel::makeCluster(n.cores_oi, type="FORK") 
    }
    
    ligand_activities_targets_geneset_ALL =  parallel::parLapply(clust, receivers_oi,function(receiver_oi, receiver_de,   verbose, ligand_target_matrix, logFC_threshold, p_val_threshold, p_val_adj, top_n_target){
      
      requireNamespace("dplyr")
      
      if(verbose == TRUE){
        print("receiver_oi:")
        print(receiver_oi %>% as.character())
      }
      
      # de_output_tidyXXX = muscat::resDS(receiver_de$sce, receiver_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()
      de_output_tidy = receiver_de
      
      de_output_tidy = de_output_tidy %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast)

      background_expressed_genes = de_output_tidy$gene %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
      ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
      ligands = colnames(ligand_target_matrix)
      
      ligand_activities_targets_geneset = de_output_tidy$contrast %>% unique() %>%
        lapply(function(contrast_oi,de_output_tidy){
          
          if(p_val_adj == TRUE){
            de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_adj < p_val_threshold)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          } else {
            de_tbl_geneset = de_output_tidy  %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_val < p_val_threshold)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          }
          
          if(verbose == TRUE){
            print("contrast_oi:")
            print(contrast_oi)
            print("Number of DE genes (gene set of interest): ")
            print(length(geneset_oi))
          }
          
          if(length(geneset_oi) > 0){
            ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
            ligand_activities = ligand_activities %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = pearson) %>% dplyr::select(-aupr, -auroc)
            
            ligand_target_df = ligand_activities$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
            ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi) %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
          } else {
            warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no DE genes - so ligand activities will be NA. Please check the DE output."))
            ligand_activities = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
          }
          
          de_genes_df = de_tbl_geneset %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(receiver = cluster_id)
          
          return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
        }, de_output_tidy)
      
      ligand_activities = ligand_activities_targets_geneset %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()
      de_genes_df = ligand_activities_targets_geneset %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()
      
      return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
      
    },
    receiver_de, verbose, ligand_target_matrix,  logFC_threshold, p_val_threshold, p_val_adj, top_n_target)
    
    parallel::stopCluster(clust)
    
  } else {
    ligand_activities_targets_geneset_ALL = lapply(receivers_oi,function(receiver_oi, receiver_de, verbose, ligand_target_matrix, logFC_threshold, p_val_threshold, p_val_adj, top_n_target){
      
      requireNamespace("dplyr")
      
      if(verbose == TRUE){
        print("receiver_oi:")
        print(receiver_oi %>% as.character())
      }
      
      de_output_tidy = receiver_de
      de_output_tidy = de_output_tidy %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast)

      background_expressed_genes = de_output_tidy$gene %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
      ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
      ligands = colnames(ligand_target_matrix)
      
      ligand_activities_targets_geneset = de_output_tidy$contrast %>% unique() %>%
        lapply(function(contrast_oi,de_output_tidy){
          
          if(p_val_adj == TRUE){
            de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_adj < p_val_threshold)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          } else {
            de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_val < p_val_threshold)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          }
          
          if(verbose == TRUE){
            print("contrast_oi:")
            print(contrast_oi)
            print("Number of DE genes (gene set of interest): ")
            print(length(geneset_oi))
          }
          
          if(length(geneset_oi) > 0){
            ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
            ligand_activities = ligand_activities %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = pearson) %>% dplyr::select(-aupr, -auroc)
            
            ligand_target_df = ligand_activities$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
            ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi) %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
          } else {
            warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no DE genes - so ligand activities will be NA. Please check the DE output."))
            ligand_activities = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
          }
          
          de_genes_df = de_tbl_geneset %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(receiver = cluster_id)
          
          return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
        }, de_output_tidy)
      
      ligand_activities = ligand_activities_targets_geneset %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()
      de_genes_df = ligand_activities_targets_geneset %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()
      
      return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
      
    },receiver_de,   verbose, ligand_target_matrix,  logFC_threshold, p_val_threshold, p_val_adj, top_n_target)
  }
  
  
  ligand_activities = ligand_activities_targets_geneset_ALL %>% purrr::map("ligand_activities") %>% dplyr::bind_rows() 
  de_genes_df = ligand_activities_targets_geneset_ALL %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()
  
  
  return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
  
}
