#' @title get_ligand_activities_targets_DEgenes
#'
#' @description \code{get_ligand_activities_targets_DEgenes}  Predict NicheNet ligand activities and ligand-target links for the receiver cell types. Uses `nichenetr::predict_ligand_activities()` and `nichenetr::get_weighted_ligand_target_links()` under the hood.
#' @usage get_ligand_activities_targets_DEgenes(receiver_de, receivers_oi,  ligand_target_matrix, logFC_threshold = 0.50, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_target = 250, verbose = FALSE, n.cores = 1)
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
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' batches = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() 
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
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
get_ligand_activities_targets_DEgenes = function(receiver_de, receivers_oi,  ligand_target_matrix, logFC_threshold = 0.50, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_target = 250, verbose = FALSE, n.cores = 1){
  
  requireNamespace("dplyr")
  receivers_oi = receiver_de$cluster_id %>% unique() %>% sort()
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
      
      de_output_tidy = receiver_de
      
      de_output_tidy = de_output_tidy %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast)
      
      background_expressed_genes = de_output_tidy$gene %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
      ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
      ligands = colnames(ligand_target_matrix)
      
      geneset_vs_ligand_activities = list()
      ligand_activities_targets_geneset = list()
      for(i in seq(length(de_output_tidy$contrast %>% unique()))){
        contrast_oi = de_output_tidy$contrast %>% unique() %>% .[i]
        if(p_val_adj == TRUE){
          de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC >= logFC_threshold & p_adj <= p_val_threshold)
          geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
          de_tbl_geneset_down = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC <= -1*logFC_threshold & p_adj <= p_val_threshold)
          geneset_oi_down = de_tbl_geneset_down %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
        } else {
          de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC >= logFC_threshold & p_val <= p_val_threshold)
          geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
          de_tbl_geneset_down = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC <= -1*logFC_threshold & p_val <= p_val_threshold)
          geneset_oi_down = de_tbl_geneset_down %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
        }
        
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
          ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi, direction_regulation = "up") %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
        } else {
          warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no upregulated DE genes - so ligand activities will be NA. Please check the DE output."))
          ligand_activities = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, direction_regulation = "up", ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
        }
        if(length(geneset_oi_down) > 0){
          if(verbose == TRUE){
            print("Number of downregulated DE genes (gene set of interest): ")
            print(length(geneset_oi_down))
          }
          
          geneset_id = geneset_oi_down %>% paste(collapse = ".")
          if(geneset_id %in% names(geneset_vs_ligand_activities)) {
            ligand_activities_down = geneset_vs_ligand_activities[[geneset_id]]$ligand_activities_df
            
          } else {
            ligand_activities_down = nichenetr::predict_ligand_activities(geneset = geneset_oi_down, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
            geneset_vs_ligand_activities[[geneset_id]] = list(ligand_activities_df = ligand_activities_down)
          }
          
          ligand_activities_down = ligand_activities_down %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = aupr_corrected) %>% dplyr::select(-pearson, -auroc, -aupr)
          
          ligand_target_df = ligand_activities_down$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi_down, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
          ligand_activities_down = ligand_activities_down %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi, direction_regulation = "down") %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
        } else {
          warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no downregulated DE genes - so ligand activities will be NA. Please check the DE output."))
          ligand_activities_down = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, direction_regulation = "down", ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
        }
        ligand_activities = ligand_activities %>% bind_rows(ligand_activities_down)
        de_genes_df = de_tbl_geneset %>% bind_rows(de_tbl_geneset_down) %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(receiver = cluster_id)
        
        ligand_activities_targets_geneset[[i]] = list(ligand_activities = ligand_activities, de_genes_df = de_genes_df)
      }
      
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
      
      geneset_vs_ligand_activities = list()
      ligand_activities_targets_geneset = list()
      for(i in seq(length(de_output_tidy$contrast %>% unique()))){
        contrast_oi = de_output_tidy$contrast %>% unique() %>% .[i]
        if(p_val_adj == TRUE){
          de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC >= logFC_threshold & p_adj <= p_val_threshold)
          geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
          de_tbl_geneset_down = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC <= -1*logFC_threshold & p_adj <= p_val_threshold)
          geneset_oi_down = de_tbl_geneset_down %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
        } else {
          de_tbl_geneset = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC >= logFC_threshold & p_val <= p_val_threshold)
          geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          
          de_tbl_geneset_down = de_output_tidy %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC <= -1*logFC_threshold & p_val <= p_val_threshold)
          geneset_oi_down = de_tbl_geneset_down %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
        }
        
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
          ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi, direction_regulation = "up") %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
        } else {
          warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no upregulated DE genes - so ligand activities will be NA. Please check the DE output."))
          ligand_activities = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, direction_regulation = "up", ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
        }
        if(length(geneset_oi_down) > 0){
          if(verbose == TRUE){
            print("Number of downregulated DE genes (gene set of interest): ")
            print(length(geneset_oi_down))
          }
          
          geneset_id = geneset_oi_down %>% paste(collapse = ".")
          if(geneset_id %in% names(geneset_vs_ligand_activities)) {
            ligand_activities_down = geneset_vs_ligand_activities[[geneset_id]]$ligand_activities_df
            
          } else {
            ligand_activities_down = nichenetr::predict_ligand_activities(geneset = geneset_oi_down, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
            geneset_vs_ligand_activities[[geneset_id]] = list(ligand_activities_df = ligand_activities_down)
          }
          
          ligand_activities_down = ligand_activities_down %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = aupr_corrected) %>% dplyr::select(-pearson, -auroc, -aupr)
          
          ligand_target_df = ligand_activities_down$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi_down, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
          ligand_activities_down = ligand_activities_down %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi, direction_regulation = "down") %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
        } else {
          warning(paste0("For Celltype ", receiver_oi, " in condition ", contrast_oi, " there seem to be no downregulated DE genes - so ligand activities will be NA. Please check the DE output."))
          ligand_activities_down = tibble(ligand = ligands, activity = NA, contrast = contrast_oi, target = NA, direction_regulation = "down", ligand_target_weight = NA, receiver = receiver_oi, activity_scaled = NA)
        }
        ligand_activities = ligand_activities %>% bind_rows(ligand_activities_down)
        de_genes_df = de_tbl_geneset %>% bind_rows(de_tbl_geneset_down) %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(receiver = cluster_id)
        
        ligand_activities_targets_geneset[[i]] = list(ligand_activities = ligand_activities, de_genes_df = de_genes_df)
      }
      
      ligand_activities = ligand_activities_targets_geneset %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()
      de_genes_df = ligand_activities_targets_geneset %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()
      
      return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
      
    },receiver_de,   verbose, ligand_target_matrix,  logFC_threshold, p_val_threshold, p_val_adj, top_n_target)
  }
  
  
  ligand_activities = ligand_activities_targets_geneset_ALL %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()  %>% dplyr::mutate(direction_regulation = factor(direction_regulation, levels = c("up","down")))
  de_genes_df = ligand_activities_targets_geneset_ALL %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()
  
  return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))  

}
