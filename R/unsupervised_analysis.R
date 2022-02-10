#' @title calculate_LR_pb_prod_matrix
#'
#' @description \code{calculate_LR_pb_prod_matrix} Calculate a matrix giving the pseudobulk expression product of each ligand-receptor-sender-receiver pair (that is expressed in at least one sample)
#' @usage calculate_LR_pb_prod_matrix(sce, sample_id, celltype_id, group_id, senders_oi, receivers_oi, fraction_cutoff, lr_network, batches)
#' 
#' @inheritParams multi_nichenet_analysis_separate
#' @inheritParams multi_nichenet_analysis_combined
#' @param group_id Metadata column indicating group/condition. If not available: NULL.
#' @param senders_oi character vector of the sender cell types of interest
#' @param receivers_oi character vector of the receiver cell types of interest
#'
#' @return lr_prod_mat: matrix giving the pseudobulk expression product of each ligand-receptor-sender-receiver pair 
#'  
#' @import dplyr
#' @importFrom generics setdiff intersect union
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' senders_oi = "CAF"
#' receivers_oi = "Malignant"
#' fraction_cutoff = 0.05
#' batches = NA
#' lr_prod_mat = calculate_LR_pb_prod_matrix(sce, sample_id, celltype_id, group_id, senders_oi, receivers_oi, fraction_cutoff, lr_network, batches)
#' }
#'
#' @export
#'
calculate_LR_pb_prod_matrix = function(sce, sample_id, celltype_id, group_id, senders_oi, receivers_oi, fraction_cutoff, lr_network, batches){
  
  requireNamespace("dplyr")
  
  ### Cell type Info
  if(!is.null(group_id)){
    celltype_info = suppressMessages(get_avg_frac_exprs_abund(
      sce = sce,
      sample_id = sample_id,
      celltype_id =  celltype_id,
      group_id = group_id, 
      batches = batches))
    
    receiver_info_ic = suppressMessages(process_info_to_ic(
      info_object = celltype_info,
      ic_type = "receiver",
      lr_network = lr_network))
    
    sender_info_ic = suppressMessages(process_info_to_ic(
      info_object = celltype_info,
      ic_type = "sender",
      lr_network = lr_network))
    
    sender_receiver_info = suppressMessages(combine_sender_receiver_info_ic(
      sender_info = sender_info_ic,
      receiver_info = receiver_info_ic,
      senders_oi = senders_oi,
      receivers_oi = receivers_oi,
      lr_network = lr_network))
  } else {
    group_id = "mock_group"
    SummarizedExperiment::colData(sce)[,group_id] = c("mock_group")
    
    celltype_info = suppressMessages(get_avg_frac_exprs_abund(
      sce = sce,
      sample_id = sample_id,
      celltype_id =  celltype_id,
      group_id = group_id, 
      batches = batches))
    
    receiver_info_ic = suppressMessages(process_info_to_ic(
      info_object = celltype_info,
      ic_type = "receiver",
      lr_network = lr_network))
    
    sender_info_ic = suppressMessages(process_info_to_ic(
      info_object = celltype_info,
      ic_type = "sender",
      lr_network = lr_network))
    
    sender_receiver_info = suppressMessages(combine_sender_receiver_info_ic(
      sender_info = sender_info_ic,
      receiver_info = receiver_info_ic,
      senders_oi = senders_oi,
      receivers_oi = receivers_oi,
      lr_network = lr_network))
  }
  
  fraction_expressing_ligand_receptor_df = sender_receiver_info$frq_df %>% dplyr::ungroup() %>% dplyr::select(sample, sender, receiver, ligand, receptor, fraction_ligand, fraction_receptor) %>% dplyr::distinct()  %>%
    dplyr::group_by(ligand, receptor, sender, receiver) %>%
    dplyr::summarise(n_samples = n(), n_expressing = sum(fraction_ligand > fraction_cutoff & fraction_receptor > fraction_cutoff)) %>%
    dplyr::mutate(fraction_expressing_ligand_receptor = n_expressing/n_samples) %>% dplyr::arrange(-fraction_expressing_ligand_receptor) %>% dplyr::select(-n_samples, -n_expressing)  %>% dplyr::ungroup()
  fraction_expressing_ligand_receptor_df = fraction_expressing_ligand_receptor_df %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_"))
  
  ids_oi = fraction_expressing_ligand_receptor_df %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% dplyr::pull(id) %>% unique()
  
  lr_prod_df = sender_receiver_info$pb_df %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sample, id, ligand_receptor_pb_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(id, ligand_receptor_pb_prod)
  
  lr_prod_mat = lr_prod_df %>% dplyr::select(-sample) %>% data.frame() %>% as.matrix()
  rownames(lr_prod_mat) = lr_prod_df$sample
  
  col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
  row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0)) %>% .[. == 0] %>% names()
  
  lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(col_remove),colnames(.) %>% generics::setdiff(col_remove)]
  
  return(lr_prod_mat)
  
}
#' @title pca_LR_pb_prod_matrix
#'
#' @description \code{pca_LR_pb_prod_matrix} Perform PCA and downstream analysis to analyze sample-to-sample heterogeneity in ligand-receptor pair expression. 
#' @usage pca_LR_pb_prod_matrix(lr_prod_mat, metadata_tbl)
#' 
#' @param lr_prod_mat Output of `calculate_LR_pb_prod_matrix`
#' @param metadata_tbl metadata tibble, with at least the column "sample_id"
#'
#'
#' @return list with following elements: \cr
#' `pca_lr_prod_mat`: PCA output \cr
#' `LR_pair_contribution_tbl`: loadings of each ligand-receptor-sender-receiver pair in contributing to each principal component (PC) \cr
#' `coordinate_tbl`: coordinates of each PC \cr
#' `metadata_correlation_tbl`: correlation of each PC with each metadata column of interest (in metadata_tbl) \cr
#'  
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom factoextra get_pca_ind get_pca_var
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' senders_oi = "CAF"
#' receivers_oi = "Malignant"
#' fraction_cutoff = 0.05
#' batches = NA
#' lr_prod_mat = calculate_LR_pb_prod_matrix(sce, sample_id, celltype_id, group_id, senders_oi, receivers_oi, fraction_cutoff, lr_network, batches)
#' metadata_tbl = SummarizedExperiment::colData(sce)[, c(sample_id, group_id)]  %>% data.frame() %>% tibble::as_tibble() %>% distinct()
#' colnames(metadata_tbl) = c("sample_id","group_id")
#' metadata_tbl2 = metadata_tbl %>% mutate(group_id = as.double(factor(group_id)))
#' output_pca = pca_LR_pb_prod_matrix(lr_prod_mat, metadata_tbl)
#' }
#'
#' @export
#'
pca_LR_pb_prod_matrix = function(lr_prod_mat, metadata_tbl){
  
  requireNamespace("dplyr")
  
  lr_prod_mat.pca = prcomp(lr_prod_mat , center = TRUE,scale. = TRUE)
  res_ind = factoextra::get_pca_ind(lr_prod_mat.pca)
  coordinate_tbl = res_ind$coord %>% data.frame() %>% tibble::rownames_to_column("sample_id") %>% tibble::as_tibble()
  var_info = factoextra::get_pca_var(lr_prod_mat.pca)
  var_info_tbl = var_info$contrib %>% data.frame() %>% tibble::rownames_to_column("LR_id") %>% tibble::as_tibble() %>% dplyr::arrange(-Dim.1)
  metadata_df = data.frame(metadata_tbl)
  rownames(metadata_df) = metadata_df[,"sample_id"]
  metadata_df = metadata_df[rownames(res_ind$coord),-1]
  cor_tbl = cor(metadata_df, res_ind$coord) %>% data.frame() %>% tibble::rownames_to_column("metadata_column") %>% tibble::as_tibble() %>% tidyr::gather(dimension, correlation, -metadata_column) %>% dplyr::group_by(metadata_column) %>% dplyr::arrange(-abs(correlation)) %>% dplyr::ungroup()
  return(list(pca_lr_prod_mat = lr_prod_mat.pca, LR_pair_contribution_tbl = var_info_tbl, coordinate_tbl = coordinate_tbl, metadata_correlation_tbl = cor_tbl))
}
