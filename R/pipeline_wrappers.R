#' @title make_abundance_plots
#'
#' @description \code{make_abundance_plots} Visualize cell type abundances.
#' @usage make_abundance_plots(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, batches = NA)
#'
#' @inheritParams multi_nichenet_analysis
#' @inheritParams combine_sender_receiver_info_ic
#' 
#' @return List containing cell type abundance plots and abundance_data data frame.
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
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()  
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() 
#' abundance_celltype_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id =  celltype_id, min_cells = 10, senders_oi = senders_oi, receivers_oi = receivers_oi)
#' }
#'
#' @export
#'
make_abundance_plots = function(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, batches = NA){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce)[,celltype_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,celltype_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,celltype_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,celltype_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,celltype_id]))))){
      stop("All the cell type labels in SummarizedExperiment::colData(sce)[,celltype_id] should be syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,group_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,group_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,group_id]))))){
      stop("All the group/condition labels in SummarizedExperiment::colData(sce)[,group_id] should be syntactically valid R names - see make.names")
    }
  }
  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,sample_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,sample_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,sample_id]))))){
      stop("All the sample_id labels in SummarizedExperiment::colData(sce)[,sample_id] should be syntactically valid R names - see make.names")
    }
  }
  
  ### Receiver abundance plots
  
  metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] %>% tibble::as_tibble()
  colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
  abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
  abundance_data_receiver = abundance_data %>% process_abund_info("receiver")
  abundance_data_sender = abundance_data %>% process_abund_info("sender")
  
  if(is.na(batches)){
    ## barplots
    # celltype proportion per sample
    abund_barplot = metadata_abundance %>% mutate(celltype_id = factor(celltype_id)) %>% ggplot() +
      aes(x = sample_id, fill = celltype_id) +
      geom_bar(position = "fill") +
      facet_grid(. ~ group_id, scales = "free", space = "free_x") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + ggtitle("Cell type proportions per sample") + ylab("proportion") + xlab("sample")
    
    
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
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination") + xlab("")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    
  } else {
    batch_oi = batches[1]
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    metadata_abundance = metadata_abundance %>% dplyr::inner_join(extra_metadata, by = "sample_id") %>% mutate(group_batch_id = paste(group_id, batch_oi, sep = "_"))
    
    abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_batch_id), by = "sample_id")
    abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
    abundance_data = abundance_data%>% dplyr::inner_join(metadata_abundance %>% distinct(sample_id, group_id, batch_oi), by = "sample_id")
    
    for(celltype_oi in abundance_data$celltype_id %>% unique()){
      n_group_batch_id = abundance_data %>% dplyr::filter(keep == TRUE & celltype_id == celltype_oi) %>% pull(group_batch_id) %>% unique() %>% length()
      n_groups = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, batch_oi), by = c("sample_id","group_id","batch_oi")) %>% dplyr::filter(keep == TRUE) %>% dplyr::pull(group_id) %>% unique() %>% length()
      n_batches = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, batch_oi), by = c("sample_id","group_id","batch_oi")) %>% dplyr::filter(keep == TRUE) %>% dplyr::pull(batch_oi) %>% unique() %>% length()
      
      if(n_group_batch_id < n_groups*n_batches){
        warning(paste("For celltype",celltype_oi,"not all group-batch combinations exist - this will likely lead to errors downstream in batch correction and DE analysis"))
      }
      
    }
    ## barplots
    # celltype proportion per sample
    abund_barplot = metadata_abundance %>% mutate(celltype_id = factor(celltype_id)) %>% ggplot() +
      aes(x = sample_id, fill = celltype_id) +
      geom_bar(position = "fill") +
      facet_grid(. ~ group_batch_id, scales = "free", space = "free_x") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + ggtitle("Cell type proportions per sample") + ylab("proportion") + xlab("sample")
    
    abund_plot = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id ~ group_batch_id, scales = "free", space = "free_x") +
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
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination") + xlab("")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_batch_id, n, group = group_batch_id, color = group_batch_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    
    
  }
  
  return(list(abund_plot_sample = abund_plot, abund_plot_group = abund_plot_boxplot, abund_barplot = abund_barplot,  abundance_data = abundance_data))
  
}
#' @title get_abundance_expression_info
#'
#' @description \code{get_abundance_expression_info} Visualize cell type abundances. Calculate the average and fraction of expression of each gene per sample and per group. Calculate relative abundances of cell types as well. Under the hood, the following functions are used: `get_avg_frac_exprs_abund`, `process_info_to_ic`, `combine_sender_receiver_info_ic`
#' @usage get_abundance_expression_info(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, lr_network, batches = NA)
#'
#' @inheritParams multi_nichenet_analysis
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
get_abundance_expression_info = function(sce, sample_id, group_id, celltype_id, min_cells, senders_oi, receivers_oi, lr_network, batches = NA){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce)[,celltype_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,celltype_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,celltype_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,celltype_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,celltype_id]))))){
      stop("All the cell type labels in SummarizedExperiment::colData(sce)[,celltype_id] should be syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,group_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,group_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,group_id]))))){
      stop("All the group/condition labels in SummarizedExperiment::colData(sce)[,group_id] should be syntactically valid R names - see make.names")
    }
  }
  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,sample_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,sample_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,sample_id]))))){
      stop("All the sample_id labels in SummarizedExperiment::colData(sce)[,sample_id] should be syntactically valid R names - see make.names")
    }
  }

  ### Receiver abundance plots

  metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)] %>% tibble::as_tibble()
  colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
  
  abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
  abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
  abundance_data_receiver = abundance_data %>% process_abund_info("receiver")
  abundance_data_sender = abundance_data %>% process_abund_info("sender")
  
  if(is.na(batches)){
    ## barplots
    # celltype proportion per sample
    abund_barplot = metadata_abundance %>% mutate(celltype_id = factor(celltype_id)) %>% ggplot() +
      aes(x = sample_id, fill = celltype_id) +
      geom_bar(position = "fill") +
      facet_grid(. ~ group_id, scales = "free", space = "free_x") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + ggtitle("Cell type proportions per sample") + ylab("proportion") + xlab("sample")
    
    
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
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination") + xlab("")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    
  } else {
    batch_oi = batches[1]
    extra_metadata = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    metadata_abundance = metadata_abundance %>% dplyr::inner_join(extra_metadata, by = "sample_id") %>% mutate(group_batch_id = paste(group_id, batch_oi, sep = "_"))
    
    abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_batch_id), by = "sample_id")
    abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
    abundance_data = abundance_data%>% dplyr::inner_join(metadata_abundance %>% distinct(sample_id, group_id, batch_oi), by = "sample_id")
    
    for(celltype_oi in abundance_data$celltype_id %>% unique()){
      n_group_batch_id = abundance_data %>% dplyr::filter(keep == TRUE & celltype_id == celltype_oi) %>% pull(group_batch_id) %>% unique() %>% length()
      n_groups = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, batch_oi), by = c("sample_id","group_id","batch_oi")) %>% dplyr::filter(keep == TRUE) %>% dplyr::pull(group_id) %>% unique() %>% length()
      n_batches = abundance_data %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id, group_id, batch_oi), by = c("sample_id","group_id","batch_oi")) %>% dplyr::filter(keep == TRUE) %>% dplyr::pull(batch_oi) %>% unique() %>% length()
      
      if(n_group_batch_id < n_groups*n_batches){
        warning(paste("For celltype",celltype_oi,"not all group-batch combinations exist - this will likely lead to errors downstream in batch correction and DE analysis"))
      }
      
    }
    ## barplots
    # celltype proportion per sample
    abund_barplot = metadata_abundance %>% mutate(celltype_id = factor(celltype_id)) %>% ggplot() +
      aes(x = sample_id, fill = celltype_id) +
      geom_bar(position = "fill") +
      facet_grid(. ~ group_batch_id, scales = "free", space = "free_x") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + ggtitle("Cell type proportions per sample") + ylab("proportion") + xlab("sample")
    
    abund_plot = abundance_data %>% ggplot(aes(sample_id, n, fill = keep)) + geom_bar(stat="identity") + scale_fill_manual(values = c("royalblue", "lightcoral")) + facet_grid(celltype_id ~ group_batch_id, scales = "free", space = "free_x") +
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
      ) + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash")  + ggtitle("Cell type abundances per sample") + ylab("# cells per sample-celltype combination") + xlab("")
    
    
    abund_plot_boxplot = abundance_data %>% ggplot(aes(group_batch_id, n, group = group_batch_id, color = group_batch_id)) + 
      geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
      scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
    

  }
  
  ### Cell type Info
  celltype_info = suppressMessages(get_avg_frac_exprs_abund(
    sce = sce,
    sample_id = sample_id,
    celltype_id =  celltype_id,
    group_id = group_id, 
    batches = batches))
  
  
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
  
  
  return(list(abund_plot_sample = abund_plot, abund_plot_group = abund_plot_boxplot, abund_barplot = abund_barplot, abundance_data_receiver = abundance_data_receiver, abundance_data_sender = abundance_data_sender, celltype_info = celltype_info, receiver_info_ic = receiver_info_ic, sender_info_ic = sender_info_ic, sender_receiver_info = sender_receiver_info))
  
}
#' @title get_DE_info
#'
#' @description \code{get_DE_info} Perform differential expression analysis via Muscat - Pseudobulking approach. Also visualize the p-value distribution. Under the hood, the following function is used: `perform_muscat_de_analysis`.
#' @usage get_DE_info(sce, sample_id, group_id, celltype_id, batches, covariates, contrasts_oi, min_cells = 10, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", findMarkers = FALSE, contrast_tbl = NULL, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7)
#'
#' @inheritParams multi_nichenet_analysis
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
#' batches = NA
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' avg_frac_exprs_abund = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id)
#' DE_info = get_DE_info(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#'}
#'
#' @export
#'
#'
get_DE_info = function(sce, sample_id, group_id, celltype_id, batches, covariates, contrasts_oi, min_cells = 10, assay_oi_pb = "counts", fun_oi_pb = "sum", de_method_oi = "edgeR", findMarkers = FALSE, contrast_tbl = NULL, filterByExpr.min.count = 7, filterByExpr.min.total.count = 15, filterByExpr.large.n = 4, filterByExpr.min.prop = 0.7){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if (class(sce) != "SingleCellExperiment") {
    stop("sce should be a SingleCellExperiment object")
  }
  if (!celltype_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("celltype_id should be a column name in the metadata dataframe of sce")
  }
  if (celltype_id != make.names(celltype_id)) {
    stop("celltype_id should be a syntactically valid R name - check make.names")
  }
  if (!sample_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("sample_id should be a column name in the metadata dataframe of sce")
  }
  if (sample_id != make.names(sample_id)) {
    stop("sample_id should be a syntactically valid R name - check make.names")
  }
  if (!group_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("group_id should be a column name in the metadata dataframe of sce")
  }
  if (group_id != make.names(group_id)) {
    stop("group_id should be a syntactically valid R name - check make.names")
  }
  
  if(is.double(SummarizedExperiment::colData(sce)[,celltype_id])){
    stop("SummarizedExperiment::colData(sce)[,celltype_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,group_id])){
    stop("SummarizedExperiment::colData(sce)[,group_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,sample_id])){
    stop("SummarizedExperiment::colData(sce)[,sample_id] should be a character vector or a factor")
  }
  
  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce)[,celltype_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,celltype_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,celltype_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,celltype_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,celltype_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,celltype_id]))))){
      stop("All the cell type labels in SummarizedExperiment::colData(sce)[,celltype_id] should be syntactically valid R names - see make.names")
    }
  }
  
  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,group_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,group_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,group_id]))))){
      stop("All the group/condition labels in SummarizedExperiment::colData(sce)[,group_id] should be syntactically valid R names - see make.names")
    }
  }
  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  } else{
    is_make_names = unique(sort(SummarizedExperiment::colData(sce)[,sample_id])) == make.names(unique(sort(SummarizedExperiment::colData(sce)[,sample_id])))
    if(sum(is_make_names) != length(unique(sort((SummarizedExperiment::colData(sce)[,sample_id]))))){
      stop("All the sample_id labels in SummarizedExperiment::colData(sce)[,sample_id] should be syntactically valid R names - see make.names")
    }
  }
  
  if(!is.character(contrasts_oi)){
    stop("contrasts should be a character vector")
  }
  
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = SummarizedExperiment::colData(sce)[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    # stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",","," ,", ", ")) %>% unlist() %>% unique()
  conditions_oi = conditions_oi[is.na(suppressWarnings(as.numeric(conditions_oi)))]
  
  if(length(contrasts_oi) != 1 | !is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }
  
  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()
  
  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  
  if(!is.na(batches)){
    if (sum(batches %in% colnames(SummarizedExperiment::colData(sce))) != length(batches) ) {
      stop("batches should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  }
  
  if(length(covariates) > 1){
    covariates_present = TRUE
    if (sum(covariates %in% colnames(SummarizedExperiment::colData(sce))) != length(covariates) ) {
      stop("covariates should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  } else {
    if(!is.na(covariates)){
      covariates_present = TRUE
      if (sum(covariates %in% colnames(SummarizedExperiment::colData(sce))) != length(covariates) ) {
        stop("covariates should be NA or all present as column name(s) in the metadata dataframe of sce")
      }
    } else {
      covariates_present = FALSE
    }
  }
  
  if(!is.character(assay_oi_pb)){
    stop("assay_oi_pb should be a character vector")
  } else {
    if(assay_oi_pb != "counts"){
      warning("are you sure you don't want to use the counts assay?")
    }
  }
  if(!is.character(fun_oi_pb)){
    stop("fun_oi_pb should be a character vector")
  }
  if(!is.character(de_method_oi)){
    stop("de_method_oi should be a character vector")
  }
  
  if(!is.double(min_cells)){
    stop("min_cells should be numeric")
  } else {
    if(min_cells <= 0) {
      warning("min_cells is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(findMarkers == TRUE){
    if(is.null(contrast_tbl)){
      stop("Please provide an input to the argument `contrast_tbl` -- see documentation")
    }
  }

  celltypes = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
  
  DE_list = celltypes %>% lapply(function(celltype_oi, sce){
    sce_oi = sce[, SummarizedExperiment::colData(sce)[,celltype_id] == celltype_oi]
    DE_result = tryCatch(
      {perform_muscat_de_analysis(sce = sce_oi, 
                                  sample_id = sample_id, 
                                  celltype_id = celltype_id, 
                                  group_id = group_id, 
                                  batches = batches, 
                                  covariates = covariates, 
                                  contrasts = contrasts_oi, 
                                  assay_oi_pb = assay_oi_pb, 
                                  fun_oi_pb = fun_oi_pb, 
                                  de_method_oi = de_method_oi, 
                                  min_cells = min_cells, 
                                  filterByExpr.min.count = filterByExpr.min.count, 
                                  filterByExpr.min.total.count = filterByExpr.min.total.count, 
                                  filterByExpr.large.n = filterByExpr.large.n, 
                                  filterByExpr.min.prop = filterByExpr.min.prop)
        }, 
      error = function(cond){
        message(paste0("perform_muscat_de_analysis errored for celltype: ", celltype_oi))
        message("Here's the original error message:")
        message(cond)
        message("")
        print(cond)
        message(paste0("perform_muscat_de_analysis errored for celltype: ", celltype_oi))
        message("")
        print("In case: Error in x[[1]]: subscript out of bounds: this likely means that there are not enough samples per group with sufficient cells of this cell type. This cell type will thus be ignored for further analyses, other cell types will still be considered.")
        return(NA) # occurs when not enough samples per group with sufficient cells in most cases, can also be due to other error messages
      })
  }, sce)
  
  celltype_de = list(
    de_output = c(DE_list %>% purrr::map("de_output")),
    de_output_tidy = DE_list %>% purrr::map("de_output_tidy") %>% bind_rows()
    )
  
  print("DE analysis is done:")
  
  print("included cell types are:")
  included_celltypes = celltypes %>% generics::intersect(celltype_de$de_output_tidy$cluster_id) %>% unique()
  print(included_celltypes)
  
  excluded_celltypes = celltypes %>% generics::setdiff(celltype_de$de_output_tidy$cluster_id) %>% unique()
  if (length(excluded_celltypes) > 0) {
    print("excluded cell types are:")
    print(excluded_celltypes)
    print("These celltypes are not considered in the analysis. After removing samples that contain less cells than the required minimal, some groups don't have 2 or more samples anymore (also relevant for groups not included in your contrasts!). As a result the analysis cannot be run. To solve this: decrease the number of min_cells or change your group_id and pool all samples that belong to groups that are not of interest! ")
  }
  if (length(excluded_celltypes) == length(celltypes)) {
    print("DE analysis did error for all cell types. This might be because of several reasons - check the original error message for this. Here are 2 common reasons in case no cell type past the filtering criteria: 1) no cell type has enough cells in >=2 samples per group. 2) problem in batch definition: not all levels of your batch are in each group - Also for groups not included in your contrasts!")
  }
  
  hist_pvals = celltype_de$de_output_tidy %>% dplyr::inner_join(celltype_de$de_output_tidy %>% dplyr::group_by(contrast,cluster_id) %>% dplyr::count(), by = c("cluster_id","contrast")) %>% 
    dplyr::mutate(cluster_id = paste0(cluster_id, "\nnr of genes: ", n)) %>% dplyr::mutate(`p-value <= 0.05` = p_val <= 0.05) %>% 
    ggplot(aes(x = p_val, fill = `p-value <= 0.05`)) + 
    geom_histogram(binwidth = 0.05,boundary=0, color = "grey35") + scale_fill_manual(values = c("grey90", "lightsteelblue1")) + 
    facet_grid(contrast~cluster_id) + ggtitle("P-value histograms") + theme_bw() 
  
  if(findMarkers == TRUE){
    celltypes = celltype_de$de_output_tidy %>% dplyr::pull(cluster_id) %>% unique()
    
    celltype_de_findmarkers = celltypes %>% lapply(function(celltype_oi, sce){
      genes_expressed = rownames(sce) ## change later if necessary for having a more decent filtering
      sce_oi = sce[intersect(rownames(sce), genes_expressed), SummarizedExperiment::colData(sce)[,celltype_id] == celltype_oi]
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
#' batches = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' DE_info = get_DE_info(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
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
#' @usage make_lite_output(multinichenet_output, top_n_LR = 2500)
#'
#' @param  multinichenet_output Output of a MultiNicheNet analysis (result of `multi_nichenet_analysis()`).
#' @param  top_n_LR top nr of LR pairs for which correlation with target genes will be calculated. Is 2500 by default. If you want to calculate correlation for all LR pairs, set this argument to NA.
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
make_lite_output = function(multinichenet_output, top_n_LR = 2500){
  
  requireNamespace("dplyr")
  
  if("celltype_info" %in% names(multinichenet_output)){
    gene_subset = generics::union( ## to filter the output, keep only genes that are expressed ligands, receptors and/or DE genes
      multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0) %>% .$ligand, 
      multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0) %>% .$receptor) %>% 
      generics::union(multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::pull(gene) %>% unique())
    
    multinichenet_output$celltype_info$avg_df = multinichenet_output$celltype_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$frq_df = multinichenet_output$celltype_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$pb_df = multinichenet_output$celltype_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
    
    multinichenet_output$celltype_info$avg_df_group = multinichenet_output$celltype_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$frq_df_group = multinichenet_output$celltype_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
    multinichenet_output$celltype_info$pb_df_group = multinichenet_output$celltype_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
    
    multinichenet_output$celltype_de = multinichenet_output$celltype_de %>% dplyr::filter(gene %in% gene_subset)
    
    ## maybe also a subset of LR-Sender-Receiver pairs?
    LR_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor  > 0) %>% dplyr::distinct(ligand, receptor, sender, receiver)
    
    if(is.na(top_n_LR)){
      LR_subset_cor = LR_subset
    } else {
      LR_subset_cor = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(group == top_group & fraction_expressing_ligand_receptor > 0) %>% dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% ungroup() %>% top_n(top_n_LR, prioritization_score) %>% dplyr::distinct(ligand, receptor, sender, receiver)  
    }
    
    # multinichenet_output$sender_receiver_info$avg_df = multinichenet_output$sender_receiver_info$avg_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # multinichenet_output$sender_receiver_info$frq_df = multinichenet_output$sender_receiver_info$frq_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # multinichenet_output$sender_receiver_info$pb_df = multinichenet_output$sender_receiver_info$pb_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # 
    # multinichenet_output$sender_receiver_info$avg_df_group = multinichenet_output$sender_receiver_info$avg_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # multinichenet_output$sender_receiver_info$frq_df_group = multinichenet_output$sender_receiver_info$frq_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # multinichenet_output$sender_receiver_info$pb_df_group = multinichenet_output$sender_receiver_info$pb_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    # 
    # multinichenet_output$sender_receiver_de = multinichenet_output$sender_receiver_de %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    multinichenet_output$prioritization_tables$group_prioritization_tbl = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$prioritization_tables$group_prioritization_table_source = multinichenet_output$prioritization_tables$group_prioritization_table_source %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    multinichenet_output$prioritization_tables$sample_prioritization_tbl = multinichenet_output$prioritization_tables$sample_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
    
    if(nrow(multinichenet_output$lr_target_prior_cor) > 0){
      multinichenet_output$lr_target_prior_cor = multinichenet_output$lr_target_prior_cor %>% dplyr::inner_join(LR_subset_cor, by = c("sender", "receiver", "ligand", "receptor")) %>% dplyr::filter(target %in% gene_subset)
    } else{
      multinichenet_output$lr_target_prior_cor = tibble()
    }
    
  } else {
    if("receiver_info" %in% names(multinichenet_output)) {
      # sender
      gene_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0) %>% .$ligand %>% unique()
      
      multinichenet_output$sender_info$avg_df = multinichenet_output$sender_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$frq_df = multinichenet_output$sender_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$pb_df = multinichenet_output$sender_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$sender_info$avg_df_group = multinichenet_output$sender_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$frq_df_group = multinichenet_output$sender_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$sender_info$pb_df_group = multinichenet_output$sender_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$sender_de = multinichenet_output$sender_de %>% dplyr::filter(gene %in% gene_subset)
      
      # receiver
      gene_subset = generics::union(
        multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% dplyr::pull(gene), 
        multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor > 0) %>% .$receptor)  
      
      multinichenet_output$receiver_info$avg_df = multinichenet_output$receiver_info$avg_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$frq_df = multinichenet_output$receiver_info$frq_df %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$pb_df = multinichenet_output$receiver_info$pb_df %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$receiver_info$avg_df_group = multinichenet_output$receiver_info$avg_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$frq_df_group = multinichenet_output$receiver_info$frq_df_group %>% dplyr::filter(gene %in% gene_subset)
      multinichenet_output$receiver_info$pb_df_group = multinichenet_output$receiver_info$pb_df_group %>% dplyr::filter(gene %in% gene_subset)
      
      multinichenet_output$receiver_de = multinichenet_output$receiver_de %>% dplyr::filter(gene %in% gene_subset)
      
      ## maybe also a subset of LR-Sender-Receiver pairs?
      LR_subset = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(fraction_expressing_ligand_receptor  > 0) %>% dplyr::distinct(ligand, receptor, sender, receiver)
      
      if(is.na(top_n_LR)){
        LR_subset_cor = LR_subset
      } else {
        LR_subset_cor = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::filter(group == top_group & fraction_expressing_ligand_receptor > 0) %>% dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score) %>% ungroup() %>% top_n(top_n_LR, prioritization_score) %>% dplyr::distinct(ligand, receptor, sender, receiver)  
      }      
      # multinichenet_output$sender_receiver_info$avg_df = multinichenet_output$sender_receiver_info$avg_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # multinichenet_output$sender_receiver_info$frq_df = multinichenet_output$sender_receiver_info$frq_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # multinichenet_output$sender_receiver_info$pb_df = multinichenet_output$sender_receiver_info$pb_df %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # 
      # multinichenet_output$sender_receiver_info$avg_df_group = multinichenet_output$sender_receiver_info$avg_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # multinichenet_output$sender_receiver_info$frq_df_group = multinichenet_output$sender_receiver_info$frq_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # multinichenet_output$sender_receiver_info$pb_df_group = multinichenet_output$sender_receiver_info$pb_df_group %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      # 
      # multinichenet_output$sender_receiver_de = multinichenet_output$sender_receiver_de %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      multinichenet_output$prioritization_tables$group_prioritization_tbl = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$prioritization_tables$group_prioritization_table_source = multinichenet_output$prioritization_tables$group_prioritization_table_source %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      multinichenet_output$prioritization_tables$sample_prioritization_tbl = multinichenet_output$prioritization_tables$sample_prioritization_tbl %>% dplyr::inner_join(LR_subset, by = c("sender", "receiver", "ligand", "receptor"))
      
      if(nrow(multinichenet_output$lr_target_prior_cor) > 0){
        multinichenet_output$lr_target_prior_cor = multinichenet_output$lr_target_prior_cor %>% dplyr::inner_join(LR_subset_cor, by = c("sender", "receiver", "ligand", "receptor")) %>% dplyr::filter(target %in% gene_subset)
      } else{
        multinichenet_output$lr_target_prior_cor = tibble()
      }
    }

  }

  return(multinichenet_output)
}
#' @title Convert aliases to official gene symbols in a SingleCellExperiment Object
#'
#' @description \code{alias_to_symbol_SCE}Convert aliases to official gene symbols in a SingleCellExperiment Object. Makes use of `nichenetr::convert_alias_to_symbols`
#' @usage alias_to_symbol_SCE(sce, organism)
#'
#' @inheritParams multi_nichenet_analysis
#' @param organism Is sce data from "mouse" or "human"
#' 
#' @return SingleCellExperiment Object
#'
#' @import nichenetr
#'
#' @examples
#' \dontrun{
#' sce = sce %>% alias_to_symbol_SCE("human")
#' }
#'
#' @export
#'
alias_to_symbol_SCE = function(sce, organism) {
  
  requireNamespace("dplyr")
  requireNamespace("nichenetr")
  
  RNA = sce@assays@data
  
  newnames = convert_alias_to_symbols(rownames(sce), organism = organism)
  
  # sometimes: there are doubles:
  doubles =  newnames %>% table() %>% .[. > 1] %>% names()
  genes_remove = (names(newnames[newnames %in% doubles]) != (newnames[newnames %in% doubles])) %>%  .[. == TRUE] %>% names()
  newnames[genes_remove] = genes_remove # set the doubles back to their old names
  
  rownames(sce) = newnames
  
  if(!is.null(RNA$counts)){
    dim_before = dim(RNA$counts)
    rownames(RNA$counts) = newnames
    RNA$counts = RNA$counts %>% .[!is.na(rownames(.)), ]
    dim_after = dim(RNA$counts)
    if(sum(dim_before != dim_after) > 0){
      print("dim counts assay changed")
      print(paste0("before: ",dim_before))
      print(paste0("after: ",dim_before))
    }
  }
  
  if(!is.null(RNA$logcounts)){
    dim_before = dim(RNA$logcounts)
    rownames(RNA$logcounts) = newnames
    RNA$logcounts = RNA$logcounts %>% .[!is.na(rownames(.)), ]
    dim_after = dim(RNA$logcounts)
    if(sum(dim_before != dim_after) > 0){
      print("dim logcounts assay changed")
      print(paste0("before: ",dim_before))
      print(paste0("after: ",dim_before))
    }
  }
  
  sce@assays@data = RNA
  return(sce)
}
#' @title make.names of all genes in a SingleCellExperiment Object
#'
#' @description \code{makenames_SCE} make.names of all genes in a SingleCellExperiment Object. Avoids missing of genes that are sometimes in original symbol, and sometimes in make.names() format
#' @usage makenames_SCE(sce)
#'
#' @inheritParams multi_nichenet_analysis
#' 
#' @return SingleCellExperiment Object
#'
#' @examples
#' \dontrun{
#' sce = sce %>% makenames_SCE()
#' }
#'
#' @export
#'
makenames_SCE = function(sce) {
  
  requireNamespace("dplyr")
  requireNamespace("nichenetr")
  
  RNA = sce@assays@data
  
  newnames = rownames(sce) %>% make.names()
  rownames(sce) = newnames
  
  if(!is.null(RNA$counts)){
    dim_before = dim(RNA$counts)
    rownames(RNA$counts) = newnames
    RNA$counts = RNA$counts %>% .[!is.na(rownames(.)), ]
    dim_after = dim(RNA$counts)
    if(sum(dim_before != dim_after) > 0){
      print("dim counts assay changed")
      print(paste0("before: ",dim_before))
      print(paste0("after: ",dim_before))
    }
  }
  
  if(!is.null(RNA$logcounts)){
    dim_before = dim(RNA$logcounts)
    rownames(RNA$logcounts) = newnames
    RNA$logcounts = RNA$logcounts %>% .[!is.na(rownames(.)), ]
    dim_after = dim(RNA$logcounts)
    if(sum(dim_before != dim_after) > 0){
      print("dim logcounts assay changed")
      print(paste0("before: ",dim_before))
      print(paste0("after: ",dim_before))
    }
  }
  
  sce@assays@data = RNA
  return(sce)
}