#' @title get_abundance_expression_info
#'
#' @description \code{get_abundance_expression_info} Visualize cell type abundances. Calculate the average and fraction of expression of each gene per sample and per group. Calculate relative abundances of cell types as well. Under the hood, the following functions are used: `get_avg_frac_exprs_abund`, `process_info_to_ic`, `combine_sender_receiver_info_ic`
#' @usage get_abundance_expression_info(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, lr_network, covariates = NA)
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @inheritParams combine_sender_receiver_info_ic
#' 
#' @return List containing cell type abundance plots, and data frames with average and fraction of expression per sample and per group, and relative cell type abundances as well.
#'
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importFrom tidyr gather
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
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()  
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() 
#' abundance_celltype_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id =  celltype_id, min_cells = 10, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network)
#' }
#'
#' @export
#'
get_abundance_expression_info = function(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, lr_network, covariates = NA){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  ### Receiver abundance plots

  metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] %>% tibble::as_tibble()
  colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
  abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
  abundance_data_receiver = abundance_data %>% process_abund_info("receiver")
  abundance_data_sender = abundance_data %>% process_abund_info("sender")
  
  if(is.na(covariates)){
    abund_plot = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id ~ group_id, scales = "free", space = "free_x") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    
  } else {
    covariate_oi = covariates[1]
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(covariate_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","covariate_oi")
    metadata_abundance = metadata_abundance %>% dplyr::inner_join(extra_metadata) %>% mutate(group_covariate_id = paste(group_id, covariate_oi, sep = "_"))
    
    abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_covariate_id), by = "sample_id")
    abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
    
    for(celltype_oi in abundance_data$celltype_id %>% unique()){
      n_group_covariate_id = abundance_data %>% dplyr::filter(keep == TRUE & celltype_id == celltype_oi) %>% pull(group_covariate_id) %>% unique() %>% length()
      n_groups = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, covariate_oi)) %>% dplyr::filter(keep == TRUE) %>% pull(group_id) %>% unique() %>% length()
      n_covariates = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, covariate_oi)) %>% dplyr::filter(keep == TRUE) %>% pull(covariate_oi) %>% unique() %>% length()
      
      if(n_group_covariate_id < n_groups*n_covariates){
        warning(paste("For celltype",celltype_oi,"not all group-covariate combinations exist - this will likely lead to errors downstream in batch correction and DE analysis"))
      }
      
    }
    
    abund_plot = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id ~ group_covariate_id, scales = "free", space = "free_x") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_covariate_id, n, group = group_covariate_id, color = group_covariate_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    

  }
  
  ### Cell type Info
  
  celltype_info = suppressMessages(get_avg_frac_exprs_abund(
    sce = sce,
    sample_id = sample_id,
    celltype_id =  celltype_id,
    group_id = group_id, 
    covariates = covariates))
  
  
  ### Link LR network to Cell type info
  
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
  
  
  return(list(abund_plot_sample = abund_plot, abund_plot_group = abund_plot_boxplot, abundance_data_receiver = abundance_data_receiver, abundance_data_sender = abundance_data_sender, celltype_info = celltype_info, receiver_info_ic = receiver_info_ic, sender_info_ic = sender_info_ic, sender_receiver_info = sender_receiver_info))
  
}
#' @title get_abundance_expression_info_separate
#'
#' @description \code{get_abundance_expression_info_separate}.Similar to `get_abundance_expression_info`, but now for sender and receiver object separately. Visualize cell type abundances. Calculate the average and fraction of expression of each gene per sample and per group. Calculate relative abundances of cell types as well. Under the hood, the following functions are used: `get_avg_frac_exprs_abund`, `process_info_to_ic`, `combine_sender_receiver_info_ic`
#' @usage get_abundance_expression_info_separate(sce_receiver, sce_sender, sample_id, group_id, celltype_id_receiver, celltype_id_sender, min_cells, senders_oi, receivers_oi, lr_network, covariates = NA)
#'
#' @inheritParams multi_nichenet_analysis_separate
#' @inheritParams combine_sender_receiver_info_ic
#' 
#' @return List containing cell type abundance plots, and data frames with average and fraction of expression per sample and per group, and relative cell type abundances as well.
#'
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sce_receiver = sce[, sce$celltype == "Malignant"]
#' sce_sender = sce[, sce$celltype == "CAF"]
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id_receiver = "celltype"
#' celltype_id_sender = "celltype"
#' senders_oi =   SummarizedExperiment::colData(sce)[,celltype_id_sender] %>% unique() 
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id_receiver] %>% unique() 
#' abundance_celltype_info = get_abundance_expression_info_separate(sce_receiver = sce_receiver, sce_sender = sce_sender, sample_id = sample_id, group_id = group_id, celltype_id_receiver = celltype_id_receiver, celltype_id_sender = celltype_id_sender, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network)
#' }
#'
#' @export
#'
get_abundance_expression_info_separate = function(sce_receiver, sce_sender, sample_id, group_id, celltype_id_receiver, celltype_id_sender, min_cells, senders_oi, receivers_oi, lr_network, covariates = NA){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  ### Receiver plots and info
  
  metadata_abundance = SummarizedExperiment::colData(sce_receiver)[,c(sample_id, group_id, celltype_id_receiver)]
  colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id_receiver")
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id_receiver) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
  abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
  abund_plot_receiver = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id_receiver ~ group_id, scales = "free", space = "free_x") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination")
  
  
  abund_plot_boxplot_receiver = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id_receiver, scales = "free") + theme_bw() + 
    scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
  
  abundance_data_receiver = abundance_data %>% process_abund_info("receiver")
  
  
  receiver_info = suppressMessages(get_avg_frac_exprs_abund(
    sce = sce_receiver,
    sample_id = sample_id,
    celltype_id =  celltype_id_receiver,
    group_id = group_id, 
    covariates = covariates))
  receiver_info_ic = suppressMessages(process_info_to_ic(
    info_object = receiver_info,
    ic_type = "receiver",
    lr_network = lr_network))
  
  ### Sender plots and info
  
  metadata_abundance = SummarizedExperiment::colData(sce_sender)[,c(sample_id, group_id, celltype_id_sender)] 
  colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id_sender")
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id_sender) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
  abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
  abund_plot_sender = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id_sender ~ group_id, scales = "free", space = "free_x") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination")
  
  
  abund_plot_boxplot_sender = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id_sender, scales = "free") + theme_bw() + 
    scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
  
  abundance_data_sender = abundance_data %>% process_abund_info("sender")
  
  
  sender_info = suppressMessages(get_avg_frac_exprs_abund(
    sce = sce_sender,
    sample_id = sample_id,
    celltype_id =  celltype_id_sender,
    group_id = group_id, 
    covariates = covariates))
  sender_info_ic = suppressMessages(process_info_to_ic(
    info_object = sender_info,
    ic_type = "sender",
    lr_network = lr_network))
  
  
  sender_receiver_info = suppressMessages(combine_sender_receiver_info_ic(
    sender_info = sender_info_ic,
    receiver_info = receiver_info_ic,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network))
  
  
  return(list(abund_plot_sample_receiver = abund_plot_receiver, abund_plot_group_receiver = abund_plot_boxplot_receiver, 
              abund_plot_sample_sender = abund_plot_sender, abund_plot_group_sender = abund_plot_boxplot_sender, 
              abundance_data_receiver = abundance_data_receiver, abundance_data_sender = abundance_data_sender, 
              receiver_info = receiver_info, sender_info = sender_info, 
              receiver_info_ic = receiver_info_ic, sender_info_ic = sender_info_ic, sender_receiver_info = sender_receiver_info))
  
}
#' @title get_DE_info
#'
#' @description \code{get_DE_info} Perform differential expression analysis via Muscat - Pseudobulking approach. Also visualize the p-value distribution. Under the hood, the following function is used: `perform_muscat_de_analysis`.
#' @usage get_DE_info(sce, sample_id, group_id, celltype_id, covariates, contrasts_oi, min_cells = 10, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", findMarkers = FALSE)
#'
#' @inheritParams multi_nichenet_analysis_combined
#' @inheritParams perform_muscat_de_analysis
#' @param contrast_tbl see explanation in multi_nichenet_analysis function -- here: only required to give as input if findMarkers = TRUE.
#' 
#' @return List with output of the differential expression analysis in 1) default format(`muscat::pbDS()`), and 2) in a tidy table format (`muscat::resDS()`) (both in the `celltype_de` slot); Histogram plot of the p-values is also returned.
#'
#' @import dplyr
#' @import muscat
#' @import ggplot2
#' @importFrom scran findMarkers
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' DE_info = get_DE_info(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#'}
#'
#' @export
#'
#'
get_DE_info = function(sce, sample_id, group_id, celltype_id, covariates, contrasts_oi, min_cells = 10, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", findMarkers = FALSE, contrast_tbl = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  celltype_de = perform_muscat_de_analysis(
    sce = sce,
    sample_id = sample_id,
    celltype_id = celltype_id,
    group_id = group_id,
    covariates = covariates,
    contrasts = contrasts_oi,
    assay_oi_pb = assay_oi_pb,
    fun_oi_pb = fun_oi_pb,
    de_method_oi = de_method_oi,
    min_cells = min_cells)
  
  hist_pvals = celltype_de$de_output_tidy %>% dplyr::inner_join(celltype_de$de_output_tidy %>% dplyr::group_by(contrast,cluster_id) %>% dplyr::count(), by = c("cluster_id","contrast")) %>% 
    dplyr::mutate(cluster_id = paste0(cluster_id, "\nnr of genes: ", n)) %>% dplyr::mutate(`p-value <= 0.05` = p_val <= 0.05) %>% 
    ggplot(aes(x = p_val, fill = `p-value <= 0.05`)) + 
    geom_histogram(binwidth = 0.05,boundary=0, color = "grey35") + scale_fill_manual(values = c("grey90", "lightsteelblue1")) + 
    facet_grid(contrast~cluster_id) + ggtitle("P-value histograms") + theme_bw() 
  
  if(findMarkers == TRUE){

    genes_filtered = celltype_de$de_output_tidy %>% dplyr::pull(gene) %>% unique()
    celltypes = celltype_de$de_output_tidy %>% dplyr::pull(cluster_id) %>% unique()
    
    celltype_de_findmarkers = celltypes %>% lapply(function(celltype_oi, sce){
      sce_oi = sce[genes_filtered, SummarizedExperiment::colData(sce)[,celltype_id] == celltype_oi]
      DE_tables_list = scran::findMarkers(sce_oi, test.type="t", groups = SummarizedExperiment::colData(sce_oi)[,group_id])
      
      conditions = names(DE_tables_list)
      DE_tables_df = conditions %>% lapply(function(condition_oi, DE_tables_list){
        DE_table_oi = DE_tables_list[[condition_oi]]
        DE_table_oi = DE_table_oi %>% data.frame() %>% tibble::rownames_to_column("gene") %>% tibble::as_tibble() %>% dplyr::mutate(cluster_id = celltype_oi, group = condition_oi) %>% dplyr::select(gene, p.value, FDR, summary.logFC, cluster_id, group)  
      }, DE_tables_list) %>% dplyr::bind_rows()
    }, sce) %>% dplyr::bind_rows() %>% dplyr::rename(logFC = summary.logFC, p_val = p.value, p_adj = FDR) %>% dplyr::inner_join(contrast_tbl, by = "group") %>% dplyr::select(gene, cluster_id, logFC, p_val, p_adj, contrast)
    
    hist_pvals_findmarkers = celltype_de_findmarkers %>% dplyr::inner_join(celltype_de_findmarkers %>% dplyr::group_by(contrast,cluster_id) %>% dplyr::count(), by = c("cluster_id","contrast")) %>% 
      dplyr::mutate(cluster_id = paste0(cluster_id, "\nnr of genes: ", n)) %>% dplyr::mutate(`p-value <= 0.05` = p_val <= 0.05) %>% 
      ggplot(aes(x = p_val, fill = `p-value <= 0.05`)) + 
      geom_histogram(binwidth = 0.05,boundary=0, color = "grey35") + scale_fill_manual(values = c("grey90", "lightsteelblue1")) + 
      facet_grid(contrast~cluster_id) + ggtitle("findMarker P-value histograms") + theme_bw() 
    
    
  } else {
    celltype_de_findmarkers = NA
    hist_pvals_findmarkers = NA
    
  }
  return(list(celltype_de = celltype_de, hist_pvals = hist_pvals, celltype_de_findmarkers = celltype_de_findmarkers, hist_pvals_findmarkers = hist_pvals_findmarkers))
  
}
#' @title get_empirical_pvals
#'
#' @description \code{get_empirical_pvals} Calculate empirical p-values based on a DE output. Show p-value distribution histograms. Under the hood, the following functions are used: `add_empirical_pval_fdr` and `get_FDR_empirical_plots_all`
#' @usage get_empirical_pvals(de_output_tidy)
#'
#' @param de_output_tidy Differential expression analysis output for the sender cell types. `de_output_tidy` slot of the output of `perform_muscat_de_analysis`.
#' 
#' @return `de_output_tidy`, but now 2 columns added with the empirical pvalues (normal and adjusted for multiple testing); Histogram plot of the empirical p-values is also returned.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' DE_info = get_DE_info(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#' DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
#' }
#' 
#' @export
#'
#'
get_empirical_pvals = function(de_output_tidy){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  de_output_tidy_emp = add_empirical_pval_fdr(de_output_tidy, plot = FALSE)
  z_distr_plots_emp_pval = get_FDR_empirical_plots_all(de_output_tidy)
  
  hist_pvals_emp = de_output_tidy_emp %>% inner_join(de_output_tidy_emp %>% group_by(contrast,cluster_id) %>% count(), by = c("cluster_id","contrast")) %>% 
    mutate(cluster_id = paste0(cluster_id, "\nnr of genes: ", n)) %>% mutate(`p-value <= 0.05` = p_emp <= 0.05) %>% 
    ggplot(aes(x = p_emp, fill = `p-value <= 0.05`)) + 
    geom_histogram(binwidth = 0.05,boundary=0, color = "grey35") + scale_fill_manual(values = c("grey90", "lightsteelblue1")) + 
    facet_grid(contrast~cluster_id) + ggtitle("Empirical P-value histograms") + theme_bw() 
  return(list(de_output_tidy_emp = de_output_tidy_emp, z_distr_plots_emp_pval = z_distr_plots_emp_pval, hist_pvals_emp = hist_pvals_emp))
}
#' @title make_lite_output
#'
#' @description \code{make_lite_output} Reduce the size of the MultiNicheNet output object (for memory efficiency), by only keeping expression information for present ligands, receptors, and genes DE in at least one probed condition.
#' @usage make_lite_output(multinichenet_output)
#'
#' @param  multinichenet_output Output of a MultiNicheNet analysis (result of `multi_nichenet_analysis()`).
#' 
#' @return multinichenet output list (= result of `multi_nichenet_analysis()`), but now filtered such that expression information is only returned for present ligands, receptors, and genes DE in at least one probed condition.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' multinichenet_output = multinichenet_output %>% make_lite_output()
#' 
#'}
#'
#' @export
#'
#'
make_lite_output = function(multinichenet_output){
  
  requireNamespace("dplyr")
  
  if("celltype_info" %in% names(multinichenet_output)){
    gene_subset = generics::union(multinichenet_output$prioritization_tables$group_prioritization_tbl$ligand, multinichenet_output$prioritization_tables$group_prioritization_tbl$receptor) %>% generics::union(multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::pull(gene) %>% unique())
    multinichenet_output$celltype_info$avg_df = multinichenet_output$celltype_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$frq_df = multinichenet_output$celltype_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$pb_df = multinichenet_output$celltype_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
    
    multinichenet_output$celltype_info$avg_df_group = multinichenet_output$celltype_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$frq_df_group = multinichenet_output$celltype_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$pb_df_group = multinichenet_output$celltype_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
    
    multinichenet_output$celltype_de = multinichenet_output$celltype_de %>% dplyr::filter(gene %in% gene_subset)
    
    ## maybe also a subset of LR-Sender-Receiver pairs?
    LR_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor  > 0) %>% dplyr::distinct(ligand, receptor, sender, receiver)
    
    multinichenet_output$sender_receiver_info$avg_df = multinichenet_output$sender_receiver_info$avg_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$sender_receiver_info$frq_df = multinichenet_output$sender_receiver_info$frq_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$sender_receiver_info$pb_df = multinichenet_output$sender_receiver_info$pb_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    multinichenet_output$sender_receiver_info$avg_df_group = multinichenet_output$sender_receiver_info$avg_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$sender_receiver_info$frq_df_group = multinichenet_output$sender_receiver_info$frq_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$sender_receiver_info$pb_df_group = multinichenet_output$sender_receiver_info$pb_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    multinichenet_output$sender_receiver_de = multinichenet_output$sender_receiver_de %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    multinichenet_output$prioritization_tables$group_prioritization_tbl = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$prioritization_tables$sample_prioritization_tbl = multinichenet_output$prioritization_tables$sample_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    multinichenet_output$lr_target_prior_cor = multinichenet_output$lr_target_prior_cor %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor")) %>% dplyr::filter(target %in% gene_subset)
    
  } else {
    if("receiver_info" %in% names(multinichenet_output)) {
      # sender
      gene_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl$ligand %>% unique()
      multinichenet_output$sender_info$avg_df = multinichenet_output$sender_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$frq_df = multinichenet_output$sender_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$pb_df = multinichenet_output$sender_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$sender_info$avg_df_group = multinichenet_output$sender_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$frq_df_group = multinichenet_output$sender_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$pb_df_group = multinichenet_output$sender_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$sender_de = multinichenet_output$sender_de %>% dplyr::filter(gene %in% gene_subset)
      
      # receiver
      gene_subset = unique(multinichenet_output$prioritization_tables$group_prioritization_tbl$receptor) %>% generics::union(multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::pull(gene) %>% unique())
      multinichenet_output$receiver_info$avg_df = multinichenet_output$receiver_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$frq_df = multinichenet_output$receiver_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$pb_df = multinichenet_output$receiver_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$receiver_info$avg_df_group = multinichenet_output$receiver_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$frq_df_group = multinichenet_output$receiver_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$pb_df_group = multinichenet_output$receiver_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$receiver_de = multinichenet_output$receiver_de %>% dplyr::filter(gene %in% gene_subset)
      
      ## maybe also a subset of LR-Sender-Receiver pairs?
      LR_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor  > 0) %>% dplyr::distinct(ligand, receptor, sender, receiver)
      
      multinichenet_output$sender_receiver_info$avg_df = multinichenet_output$sender_receiver_info$avg_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$sender_receiver_info$frq_df = multinichenet_output$sender_receiver_info$frq_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$sender_receiver_info$pb_df = multinichenet_output$sender_receiver_info$pb_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      multinichenet_output$sender_receiver_info$avg_df_group = multinichenet_output$sender_receiver_info$avg_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$sender_receiver_info$frq_df_group = multinichenet_output$sender_receiver_info$frq_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$sender_receiver_info$pb_df_group = multinichenet_output$sender_receiver_info$pb_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      multinichenet_output$sender_receiver_de = multinichenet_output$sender_receiver_de %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      multinichenet_output$prioritization_tables$group_prioritization_tbl = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$prioritization_tables$sample_prioritization_tbl = multinichenet_output$prioritization_tables$sample_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      multinichenet_output$lr_target_prior_cor = multinichenet_output$lr_target_prior_cor %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor")) %>% dplyr::filter(target %in% gene_subset)
      
    }

  }
  

  
  return(multinichenet_output)
}
