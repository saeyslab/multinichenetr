#' @title lr_target_prior_cor_inference
#'
#' @description \code{lr_target_prior_cor_inference} Calculate the pearson and spearman expression correlation between ligand-receptor pair pseudobulk expression products and DE gene expression product. Add the NicheNet ligand-target regulatory potential score as well as prior knowledge support for the LR-->Target link.
#' @usage lr_target_prior_cor_inference(receivers_oi, abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_LR = 2500)
#'
#' @param receivers_oi Character vector with the names of the receiver cell types of interest
#' @param abundance_expression_info Output of `get_abundance_expression_info_separate` or `get_abundance_expression_info`
#' @param celltype_de Output of `perform_muscat_de_analysis`
#' @param  top_n_LR top nr of LR pairs for which correlation with target genes will be calculated. Is 2500 by default. If you want to calculate correlation for all LR pairs, set this argument to NA.
#' @inheritParams make_sample_lr_prod_plots
#' @inheritParams get_ligand_activities_targets_DEgenes
#' @inheritParams generate_prioritization_tables
#'
#' @return Tibble with expression correlation and prior knowledge support measures for ligand-receptor to target gene links
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom generics intersect
#' @importFrom Hmisc rcorr
#' @importFrom nichenetr scale_quantile
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
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' 
#' metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] 
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' abundance_data_receiver = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' abundance_data_sender = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#' 
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' 
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#' 
#' abundance_expression_info = list(abund_plot_sample = abund_plot, abund_plot_group = abund_plot_boxplot, abundance_data_receiver = abundance_data_receiver, abundance_data_sender = abundance_data_sender, celltype_info = celltype_info, receiver_info_ic = receiver_info_ic, sender_info_ic = sender_info_ic, sender_receiver_info = sender_receiver_info)
#' 
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    contrasts = contrasts_oi)
#'    
#' sender_receiver_de = combine_sender_receiver_de(
#'  sender_de = celltype_de,
#'  receiver_de = celltype_de,
#'  senders_oi = senders_oi,
#'  receivers_oi = receivers_oi,
#'  lr_network = lr_network)
#'  
#' ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
#'    receiver_de = celltype_de,
#'    receivers_oi = receivers_oi,
#'    receiver_frq_df_group = celltype_info$frq_df_group,
#'    ligand_target_matrix = ligand_target_matrix)
#' 
#' 
#' sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
#' metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() 
#' grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
#' colnames(grouping_tbl) = c("sample","group") 
#' 
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 1,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0)
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     prioritizing_weights = prioritizing_weights,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#'
#' receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique()
#' lr_target_prior_cor = lr_target_prior_cor_inference(receivers_oi, abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix)
#' 
#' }
#'
#' @export
#'
lr_target_prior_cor_inference = function(receivers_oi, abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = 0.25, p_val_threshold = 0.05, p_val_adj = FALSE, top_n_LR = 2500){
  
  requireNamespace("dplyr")
  
  # add sender-receiver presence to grouping_tbl
  grouping_tbl = grouping_tbl %>% dplyr::inner_join(prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, sender, receiver, keep_receiver, keep_sender), by = "sample")
  
  # Step1: calculate LR prod matrix for the LR ids of interest
  if(is.na(top_n_LR)){
    ids_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0)  %>% dplyr::pull(id) %>% unique()
  } else {
    ids_oi = get_top_n_lr_pairs(prioritization_tables, top_n_LR, rank_per_group = FALSE)  %>% dplyr::pull(id) %>% unique() ## prioritization-based subset of IDs! -- only those interesting for correlation analyses with target genes
  }
  
  # make ligand-receptor-id mapping
  lig_rec_send_rec_mapping = abundance_expression_info$sender_receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample","sender","receiver")) %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(sender, receiver, ligand, receptor, id) %>% dplyr::distinct()
  
  # make receiver-id mapping to filter later on
  receiver_lr_id_mapping = abundance_expression_info$sender_receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample","sender","receiver")) %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::select(receiver, sample, id, ligand_receptor_pb_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct(receiver, id) 
  
  # make LR prod matrix
  lr_prod_df = abundance_expression_info$sender_receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample","sender","receiver")) %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1) %>% dplyr::select(sample, id, ligand_receptor_pb_prod) %>% dplyr::filter(id %in% ids_oi) %>% dplyr::distinct() %>% tidyr::spread(sample, ligand_receptor_pb_prod)
  
  lr_prod_mat = lr_prod_df %>% dplyr::select(-id) %>% data.frame() %>% as.matrix()
  rownames(lr_prod_mat) = lr_prod_df$id
  
  col_remove = lr_prod_mat %>% apply(2,function(x)sum(x != 0, na.rm = TRUE)) %>% .[. == 0] %>% names()
  row_remove = lr_prod_mat %>% apply(1,function(x)sum(x != 0, na.rm = TRUE)) %>% .[. == 0] %>% names()
  
  lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(row_remove) %>% generics::intersect(receiver_lr_id_mapping$id), colnames(.) %>% generics::setdiff(col_remove)]
  # lr_prod_mat = lr_prod_mat %>% .[rownames(.) %>% generics::setdiff(row_remove), colnames(.) %>% generics::setdiff(col_remove)]
  # lr_prod_mat = lr_prod_mat[receiver_lr_id_mapping$id, ]
  # 
  # Step2: per receiver: subset lr_prod_mat, get DE genes, calculate correlation between LR and Target expression
  lr_target_cor = receivers_oi %>% lapply(function(receiver_oi){
    
    # subset lr_prod_mat
    lr_prod_mat_oi = lr_prod_mat[receiver_lr_id_mapping %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::pull(id) %>% generics::intersect(rownames(lr_prod_mat)),]
    # print(dim(lr_prod_mat_oi))
    # get DE genes
    if(p_val_adj == FALSE){
      targets_oi = celltype_de %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::filter(p_val <= p_val_threshold & abs(logFC) >= logFC_threshold) %>% dplyr::pull(gene)

    } else {
      targets_oi = celltype_de %>% dplyr::filter(cluster_id == receiver_oi) %>% dplyr::filter(p_adj <= p_val_threshold & abs(logFC) >= logFC_threshold) %>% dplyr::pull(gene)
    }

    if("celltype_info" %in% names(abundance_expression_info)){
      pb_df =  abundance_expression_info$celltype_info$pb_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi) %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1)
    }
    if("receiver_info" %in% names(abundance_expression_info)){
      pb_df =  abundance_expression_info$receiver_info$pb_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi) %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1)
    }
    
    ##### target genes correlation #####
    
    target_df = pb_df %>% dplyr::select(sample, gene, pb_sample) %>% dplyr::distinct() %>% tidyr::spread(sample, pb_sample)
    target_mat = target_df %>% dplyr::select(-gene) %>% data.frame() %>% as.matrix()
    # print(dim(target_mat))
    
    rownames(target_mat) = target_df$gene
    
    col_remove = target_mat %>% apply(2,function(x)sum(x != 0, na.rm = TRUE)) %>% .[. == 0] %>% names()
    row_remove = target_mat %>% apply(1,function(x)sum(x != 0, na.rm = TRUE)) %>% .[. == 0] %>% names()
    
    target_mat = target_mat %>% .[rownames(.) %>% generics::setdiff(row_remove) %>% generics::intersect(targets_oi), colnames(.) %>% generics::setdiff(col_remove)]
    # target_mat = target_mat %>% .[rownames(.) %>% generics::setdiff(row_remove) %>% generics::intersect(targets_oi),]
    # target_mat = target_mat[targets_oi, ]
    # print(dim(target_mat))
    
    #  calculate correlation between LR and Target expression
    # make sure the dimensions of both matrices are the same
    common_samples = intersect(colnames(lr_prod_mat_oi), colnames(target_mat))
    if(length(common_samples) < 5){
      warning(paste0("not enough samples for a correlation analysis for the celltype ",receiver_oi))
      cor_df = NULL
    } else {
      lr_prod_mat_oi = lr_prod_mat_oi[,common_samples]
      target_mat = target_mat[,common_samples]
      
      # pearson
      cor_mat = Hmisc::rcorr(lr_prod_mat_oi %>% t(), target_mat %>% t())
      
      cor_df_pearson = cor_mat$r %>% .[,rownames(target_mat)] %>% data.frame() %>% tibble::rownames_to_column("id") %>% tidyr::gather(target, pearson, -id) %>% tibble::as_tibble()
      cor_df_pearson_pval = cor_mat$P %>% .[,rownames(target_mat)] %>% data.frame() %>% tibble::rownames_to_column("id") %>% tidyr::gather(target, pearson_pval, -id) %>% tibble::as_tibble()
      
      # spearman
      cor_mat = Hmisc::rcorr(lr_prod_mat_oi %>% t(), target_mat %>% t(), type = "spearman")
      
      cor_df_spearman = cor_mat$r %>% .[,rownames(target_mat)] %>% data.frame() %>% tibble::rownames_to_column("id") %>% tidyr::gather(target, spearman, -id) %>% tibble::as_tibble()
      cor_df_spearman_pval = cor_mat$P %>% .[,rownames(target_mat)] %>% data.frame() %>% tibble::rownames_to_column("id") %>% tidyr::gather(target, spearman_pval, -id) %>% tibble::as_tibble()
      
      cor_df = lig_rec_send_rec_mapping %>% dplyr::inner_join(cor_df_pearson, by = "id") %>% dplyr::inner_join(cor_df_pearson_pval, by = c("id", "target")) %>% dplyr::inner_join(cor_df_spearman, by = c("id", "target")) %>% dplyr::inner_join(cor_df_spearman_pval, by = c("id", "target")) 
    }
    return(cor_df)
    # print(cor_df)
    # scaling of the correlation metric -- Don't do this for now
    #   cor_df = cor_df %>% dplyr::ungroup() %>% dplyr::mutate(scaled_pearson = nichenetr::scale_quantile(pearson, 0.05), scaled_spearman = nichenetr::scale_quantile(spearman, 0.05))  # is this  scaling necessary? 
  }) %>% bind_rows()
  # print(lr_target_cor %>% head()) ## TOREMOVE
  
  # Step3: Scale the ligand-target prior information scores 
  ligand_target_df = ligand_target_matrix %>% data.frame() %>% tibble::rownames_to_column("target") %>% tidyr::gather(ligand, prior_score, -target) %>% tibble::as_tibble()
  ligand_target_df = ligand_target_df %>% dplyr::group_by(ligand) %>% dplyr::mutate(rank_of_target = rank(desc(prior_score), ties.method = "min")) %>% dplyr::group_by(target) %>% dplyr::mutate(rank_of_ligand = rank(desc(prior_score), ties.method = "min"))

  # ligand_target_df = ligand_target_df %>% dplyr::group_by(ligand) %>% dplyr::mutate(scaled_ligand = nichenetr::scale_quantile(prior_score, 0.0005)) %>% dplyr::group_by(target) %>% dplyr::mutate(scaled_target = nichenetr::scale_quantile(prior_score, 0.0005))
  # ligand_target_df = ligand_target_df %>% dplyr::mutate(scaled_prior_score = 0.5*(scaled_ligand + scaled_target))
  
  # Step4: Combine the ligand-target prior information scores and combine with the correlation based ones!
  if(nrow(lr_target_cor) > 0){
    cor_prior_df = lr_target_cor %>% dplyr::inner_join(ligand_target_df, by = c("ligand", "target")) %>% dplyr::mutate(id_target = paste(id, target, sep = "_")) %>% dplyr::ungroup() 
  } else {
    cor_prior_df = NULL
    print("For no celltypes, sufficient samples (>= 5) were available for a correlation analysis. lr_target_prior_cor, the output of this function, will be NULL. As a result, not all types of downstream visualizations can be created.")
  }
  
  
  # combine prior and correlation information
  # cor_prior_df = cor_prior_df %>% mutate(scaled_cor_score = 0.5*(scaled_pearson + scaled_target_pearson ))
  # cor_prior_df = cor_prior_df %>% mutate(final_score = 0.50*(2*scaled_prior_score + scaled_cor_score)) %>% arrange(-final_score)
  # cor_prior_df = cor_prior_df %>% dplyr::mutate(final_score = (scaled_prior_score + 0.50*scaled_pearson + 0.50*scaled_spearman)/2) %>% dplyr::arrange(-final_score) # I could give more weight to the prior information? - but maybe do not do this for the moment?
  return(cor_prior_df)
  
}

