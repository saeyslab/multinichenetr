#' @title get_ligand_activities_targets_DEgenes
#'
#' @description \code{get_ligand_activities_targets_DEgenes}  XXXX
#' @usage get_ligand_activities_targets_DEgenes(receiver_de, receivers_oi, receiver_frq_df_group, ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, frac_cutoff = 0.05, p_val_adj = FALSE, top_n_target = 250)
#'
#' @inheritParams ms_mg_nichenet_analysis_separate
#' @inheritParams combine_sender_receiver_info_ic
#' @inheritParams combine_sender_receiver_de
#' @param receiver_frq_df_group XXX
#'
#' @return XXXX
#'
#' @import Seurat
#' @import dplyr
#' @import muscat
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
get_ligand_activities_targets_DEgenes = function(receiver_de, receivers_oi, receiver_frq_df_group, ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, frac_cutoff = 0.05, p_val_adj = FALSE, top_n_target = 250){

  # consider other dplyr::filtering of the target genes? based on fraction of samples in a group that should have enough expression?

  ligand_activities_targets_geneset_ALL = receivers_oi %>%
    lapply(function(receiver_oi, receiver_de, receiver_frq_df_group, frac_cutoff){
      print("receiver_oi:")
      print(receiver_oi %>% as.character())

      de_output_tidy = muscat::resDS(receiver_de$sce, receiver_de$de_output, bind = "row", cpm = FALSE, frq = FALSE) %>% tibble::as_tibble()
      de_output_tidy = de_output_tidy %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj.loc, contrast)

      background_expressed_genes = de_output_tidy$gene %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
      ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
      ligands = colnames(ligand_target_matrix)

      frq_tbl = receiver_frq_df_group %>% dplyr::group_by(gene, celltype) %>% dplyr::summarise(max_frac = max(fraction_group))  %>% dplyr::mutate(present = max_frac > frac_cutoff) %>% dplyr::select(gene, celltype, present, max_frac) %>% dplyr::rename(cluster_id = celltype)

      ligand_activities_targets_geneset = de_output_tidy$contrast %>% unique() %>%
        lapply(function(contrast_oi,de_output_tidy){
          print("contrast_oi:")
          print(contrast_oi)
          if(p_val_adj == TRUE){
            de_tbl_geneset = de_output_tidy %>% dplyr::inner_join(frq_tbl) %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_adj.loc < p_val_threshold & present)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          } else {
            de_tbl_geneset = de_output_tidy %>% dplyr::inner_join(frq_tbl) %>% dplyr::filter(contrast == contrast_oi) %>% dplyr::filter(logFC > logFC_threshold & p_val < p_val_threshold & present)
            geneset_oi = de_tbl_geneset %>% dplyr::pull(gene) %>% unique() %>% dplyr::intersect(rownames(ligand_target_matrix))
          }

          ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
          ligand_activities = ligand_activities %>% dplyr::mutate(contrast = contrast_oi) %>% tidyr::drop_na() %>% dplyr::rename(ligand = test_ligand, activity = pearson) %>% dplyr::select(-aupr, -auroc)

          ligand_target_df = ligand_activities$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows() %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(ligand_target_weight = weight)
          ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df) %>% dplyr::mutate(receiver = receiver_oi)

          de_genes_df = de_tbl_geneset %>% dplyr::mutate(contrast = contrast_oi) %>% dplyr::rename(receiver = cluster_id)

          return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))
        }, de_output_tidy)

      ligand_activities = ligand_activities_targets_geneset %>% purrr::map("ligand_activities") %>% dplyr::bind_rows()
      de_genes_df = ligand_activities_targets_geneset %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()

      return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))

    },receiver_de, receiver_frq_df_group, frac_cutoff)


  ligand_activities = ligand_activities_targets_geneset_ALL %>% purrr::map("ligand_activities") %>% dplyr::bind_rows() %>% dplyr::group_by(receiver, contrast) %>% dplyr::mutate(activity_scaled = nichenetr::scaling_zscore(activity))
  de_genes_df = ligand_activities_targets_geneset_ALL %>% purrr::map("de_genes_df") %>% dplyr::bind_rows()


  return(list(ligand_activities = ligand_activities, de_genes_df = de_genes_df))

}
