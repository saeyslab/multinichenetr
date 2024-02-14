#' @title make_sample_lr_prod_plots
#'
#' @description \code{make_sample_lr_prod_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups
#' @usage make_sample_lr_prod_plots(prioritization_tables, prioritized_tbl_oi)
#'
#' @param prioritization_tables Output of `generate_prioritization_tables` or sublist in the output of `multi_nichenet_analysis`
#' @param prioritized_tbl_oi Subset of `prioritization_tables$group_prioritization_tbl`: the ligand-receptor interactions shown in this subset will be visualized: recommended to consider the top n LR interactions of a group of interest, based on the prioritization_score (eg n = 50; see vignettes for examples). 
#'
#' @return Ligand-Receptor Expression Product Dotplot/Heatmap
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score) 
#' plot_oi = make_sample_lr_prod_plots(output$prioritization_tables, prioritized_tbl_oi)
#' plot_oi
#'}
#'
#' @export
#'
make_sample_lr_prod_plots = function(prioritization_tables, prioritized_tbl_oi){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  pb_exprs_data_subset = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  pb_exprs_data_subset = pb_exprs_data_subset %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = pb_exprs_data_subset$sender_receiver %>% unique()))
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4.25)
  names(keep_sender_receiver_values) = levels(pb_exprs_data_subset$keep_sender_receiver)
  
  p1 = pb_exprs_data_subset %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + xlab("") + ylab("") +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(pb_exprs_data_subset$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill

  return(p1)

}

#' @title make_sample_lr_prod_activity_plots
#'
#' @description \code{make_sample_lr_prod_activity_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups. In addition, show the NicheNet ligand activities in each receiver-celltype combination.
#' @usage make_sample_lr_prod_activity_plots(prioritization_tables, prioritized_tbl_oi, widths = NULL)
#'
#' @inheritParams make_sample_lr_prod_plots
#' @param widths Vector of 4 elements: Width of the LR exprs product panel,  width of the scaled ligand activity panel, width of the ligand activity panel & width of the cell-type specificity panel. Default NULL: automatically defined based on nr of samples vs nr of group-receiver combinations. If manual change: example format: c(6,2,2,1) 
#'
#' @return Ligand-Receptor Expression Product & Ligand Activities Dotplot/Heatmap 
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#' @importFrom viridis scale_color_viridis
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score) 
#' plot_oi = make_sample_lr_prod_activity_plots(output$prioritization_tables, prioritized_tbl_oi)
#' plot_oi
#' }
#'
#' @export
#'
make_sample_lr_prod_activity_plots = function(prioritization_tables, prioritized_tbl_oi, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))   %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  sample_data = sample_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sample_data$sender_receiver %>% unique()))
  
  group_data = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, activity, activity_scaled, direction_regulation, prioritization_score) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data_celltype_specificity = prioritization_tables$group_prioritization_tbl  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_celltype_specificity = group_data_celltype_specificity %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data_celltype_specificity$sender_receiver %>% unique()))
  
  group_data_frac_expression = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_frac_expression = group_data_frac_expression %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data = group_data %>% inner_join(group_data_celltype_specificity) %>% inner_join(group_data_frac_expression)
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  rm(group_data_celltype_specificity)
  rm(group_data_frac_expression)
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
  
  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  p2 = group_data %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  p3 = group_data %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand\nActivity in Receiver")
  max_activity = (group_data$activity) %>% max()
  min_activity = (group_data$activity) %>% min()
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  
  p3 = p3 + custom_scale_fill
  
  # add the plot visualizing cell-type specificity
  # cs_data = group_data %>% filter(group %in% prioritized_tbl_oi$group) %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data = group_data %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
  frac_data = group_data %>% distinct(sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% tidyr::gather(LR, fraction_expression, fraction_ligand_group:fraction_receptor_group)
  frac_data$LR[frac_data$LR == "fraction_ligand_group"] = "ligand"
  frac_data$LR[frac_data$LR == "fraction_receptor_group"] = "receptor"
  
  cs_data = cs_data %>% inner_join(frac_data)
  
  p_cs = cs_data %>% 
    ggplot(aes(LR , lr_interaction, color = celltype_specificity, size = fraction_expression)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    viridis::scale_color_viridis() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled celltype specificity") + labs(size = "Fraction of expression")
  
  
  if(!is.null(widths)){
    p = patchwork::wrap_plots(
      p1,p2,p3,p_cs,
      nrow = 1,guides = "collect",
      widths = widths
    )
  } else {
    p = patchwork::wrap_plots(
      p1,p2,p3,p_cs,
      nrow = 1,guides = "collect",
      widths = c(sample_data$sample %>% unique() %>% length(), 2*(sample_data$group %>% unique() %>% length()), 2*(sample_data$group %>% unique() %>% length()),2*(sample_data$group %>% unique() %>% length()))
    )
  }
  
  return(p)
  
}
#' @title make_sample_lr_prod_activity_batch_plots
#'
#' @description \code{make_sample_lr_prod_activity_batch_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups. In addition, show the NicheNet ligand activities in each receiver-celltype combination. On top of this summary plot, a heatmap indicates the batch value for each displayed sample.
#' @usage make_sample_lr_prod_activity_batch_plots(prioritization_tables, prioritized_tbl_oi, grouping_tbl, batch_oi, widths = NULL, heights = NULL)
#'
#' @inheritParams make_sample_lr_prod_activity_plots
#' @param grouping_tbl Data frame linking the sample_id, group_id and batch_oi
#' @param batch_oi Name of the batch that needs to be visualized for each sample
#' @param heights Vector of 2 elements: Height of the batch panel and height of the ligand-receptor prod+activity panel. Default NULL: automatically defined based on the nr of Ligand-Receptor pairs. If manual change: example format: c(1,5) 
#'
#' @return Ligand-Receptor Expression Product & Ligand Activities Dotplot/Heatmap, complemented with a heatmap indicating the batch of interest
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
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
#' batches = "batch"
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' batch_oi = "batch"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score) 
#' plot_oi = make_sample_lr_prod_activity_batch_plots(output$prioritization_tables, prioritized_tbl_oi, output$grouping_tbl, batch_oi = batch_oi)
#' plot_oi
#' }
#'
#' @export
#'
make_sample_lr_prod_activity_batch_plots = function(prioritization_tables, prioritized_tbl_oi, grouping_tbl, batch_oi, widths = NULL, heights = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))   %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  sample_data = sample_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sample_data$sender_receiver %>% unique()))
  
  group_data = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, activity, activity_scaled, direction_regulation, prioritization_score) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data_celltype_specificity = prioritization_tables$group_prioritization_tbl  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_celltype_specificity = group_data_celltype_specificity %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data_celltype_specificity$sender_receiver %>% unique()))
  
  group_data_frac_expression = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_frac_expression = group_data_frac_expression %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data = group_data %>% inner_join(group_data_celltype_specificity) %>% inner_join(group_data_frac_expression)
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  rm(group_data_celltype_specificity)
  rm(group_data_frac_expression)
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
  
  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  p2 = group_data %>%
    ggplot(aes(direction_regulation, lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  p3 = group_data %>%
    ggplot(aes(direction_regulation, lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand\nActivity in Receiver")
  max_activity = (group_data$activity) %>% max()
  min_activity = (group_data$activity) %>% min()
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  
  p3 = p3 + custom_scale_fill
  
  grouping_tbl_plot = grouping_tbl %>% mutate(batch_ = paste0(" ",batch_oi," "), mock = "", batch = grouping_tbl[[batch_oi]])
  p_batch = grouping_tbl_plot %>% 
    ggplot(aes(sample, mock, fill = batch)) +
    geom_tile(color = "black") +
    facet_grid(batch_~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + scale_fill_brewer(palette = "Set2") + xlab("")
  
  # add the plot visualizing cell-type specificity
  cs_data = group_data %>% filter(group %in% prioritized_tbl_oi$group) %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
  
  # add the plot visualizing cell-type specificity
  # cs_data = group_data %>% filter(group %in% prioritized_tbl_oi$group) %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data = group_data %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
  frac_data = group_data %>% distinct(sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% tidyr::gather(LR, fraction_expression, fraction_ligand_group:fraction_receptor_group)
  frac_data$LR[frac_data$LR == "fraction_ligand_group"] = "ligand"
  frac_data$LR[frac_data$LR == "fraction_receptor_group"] = "receptor"
  
  cs_data = cs_data %>% inner_join(frac_data)
  
  p_cs = cs_data %>% 
    ggplot(aes(LR , lr_interaction, color = celltype_specificity, size = fraction_expression)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    viridis::scale_color_viridis() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled celltype specificity") + labs(size = "Fraction of expression")
  
  
  if(is.null(widths)){
    widths = c(sample_data$sample %>% unique() %>% length(), 2*(sample_data$group %>% unique() %>% length()), 2*(sample_data$group %>% unique() %>% length()), 2*(sample_data$group %>% unique() %>% length()))
  }
  if(is.null(heights)){
    heights = c(1, group_data$id %>% unique() %>% length())
  }
  design <- "D###
             ABCE"
  p = patchwork::wrap_plots(
    A = p1, B = p2, C= p3, D = p_batch, E = p_cs,
    guides = "collect", design = design,
    widths = widths,
    heights = heights
  ) + patchwork::plot_layout(design = design)
  
  return(p)
  
}

#' @title make_ligand_activity_plots
#'
#' @description \code{make_ligand_activity_plots}  Visualize the ligand activities (normal and scaled) of each group-receiver combination
#' @usage make_ligand_activity_plots(prioritization_tables, ligands_oi, contrast_tbl, widths = NULL)
#'
#' @inheritParams make_sample_lr_prod_plots
#' @param ligands_oi Character vector of ligands for which the activities should be visualized
#' @param contrast_tbl Table to link the contrast definitions to the group ids.
#' @param widths Vector of 2 elements: Width of the scaled ligand activity panel, width of the ligand activity panel. Default NULL: automatically defined based number of group-receiver combinations. If manual change: example format: c(3,2) 
#'
#' @return Heatmap of ligand activities (normal and scaled) of each group-receiver combination
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' 
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' ligands_oi = output$prioritization_tables$ligand_activities_target_de_tbl %>% inner_join(contrast_tbl) %>% group_by(group, receiver) %>% distinct(ligand, receiver, group, activity) %>% top_n(5, activity) %>% pull(ligand) %>% unique()
#' plot_oi = make_ligand_activity_plots(output$prioritization_tables, ligands_oi, contrast_tbl)
#' plot_oi
#' }
#'
#' @export
#'
make_ligand_activity_plots = function(prioritization_tables, ligands_oi, contrast_tbl, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  group_data = prioritization_tables$ligand_activities_target_de_tbl %>% dplyr::inner_join(contrast_tbl) %>% dplyr::distinct(group, ligand, receiver, activity, activity_scaled, direction_regulation) %>% dplyr::filter(ligand %in% ligands_oi) %>% dplyr::mutate(group_direction_regulation = paste(group, direction_regulation,sep = "-"))

  p1 = group_data %>%
    ggplot(aes(receiver, ligand, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group_direction_regulation, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      # axis.text.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))

  p1 = p1 + custom_scale_fill

  p2 = group_data %>%
    ggplot(aes(receiver, ligand, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group_direction_regulation, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand\nActivity in Receiver")
  max_activity = (group_data$activity) %>% max()
  min_activity = (group_data$activity) %>% min()
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  
  p2 = p2 + custom_scale_fill


  if(!is.null(widths)){
    p = patchwork::wrap_plots(
      p1,p2,
      nrow = 1,guides = "collect",
      widths = widths
    )
  } else {
    p = patchwork::wrap_plots(
      p1,p2,
      nrow = 1,guides = "collect",
      widths = c(group_data$receiver %>% unique() %>% length(), group_data$receiver %>% unique() %>% length())
    )
  }


  return(p)


}

#' @title make_DEgene_dotplot_pseudobulk
#'
#' @description \code{make_DEgene_dotplot_pseudobulk}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes in rows, samples in columns
#' @usage make_DEgene_dotplot_pseudobulk(genes_oi, celltype_info, prioritization_tables, celltype_oi, grouping_tbl, groups_oi = NULL)
#'
#' @param genes_oi Character vector with names of genes to visualize
#' @param celltype_info `celltype_info` or `receiver_info` slot of the output of the `multi_nichenet_analysis` function 
#' @param prioritization_tables `prioritization_tables` slot of the output of the `generate_prioritization_tables` or `multi_nichenet_analysis` function
#' @param celltype_oi Character vector with names of celltype of interest
#' @param grouping_tbl `grouping_tbl` slot of the output of the `multi_nichenet_analysis` function 
#' @param groups_oi Which groups to show? Default: NULL -- will show all groups.
#'
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' receiver_oi = "Malignant"
#' targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
#' p_target = make_DEgene_dotplot_pseudobulk(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, output$grouping_tbl)
#'}
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk = function(genes_oi, celltype_info, prioritization_tables, celltype_oi, grouping_tbl, groups_oi = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####
  
  keep_tbl = prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, group, receiver, keep_receiver) %>% dplyr::rename(celltype = receiver) %>% dplyr::mutate(keep_receiver = as.logical(keep_receiver))
  
  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)
  
  plot_data = celltype_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample")) %>% dplyr::inner_join(keep_tbl, by = c("sample","group","celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)
  
  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  
  if(!is.null(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)
    
  }
  
  p1 = plot_data %>%
    ggplot(aes(sample, gene, color = scaled_gene_exprs, size = keep_receiver)) +
    geom_point() +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("Samples") +ylab("Genes")  +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  ####  make the plot that indicates fraction of expression ####
  
  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype","group"))
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)
  
  p2 = plot_data %>%
    # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(sample, gene , color = scaled_gene_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  +
    labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("Samples") +ylab("Genes") 
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p2 = p2 + custom_scale_fill
  
  return(list(pseudobulk_plot = p1, singlecell_plot = p2))
  
}
#' @title make_DEgene_dotplot_pseudobulk_reversed
#'
#' @description \code{make_DEgene_dotplot_pseudobulk_reversed}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes and sample positions are reversed compared to `make_DEgene_dotplot_pseudobulk`: genes in columns, samples in rows.
#' @usage make_DEgene_dotplot_pseudobulk_reversed(genes_oi, celltype_info, prioritization_tables, celltype_oi, grouping_tbl, groups_oi = NULL, target_regulation_df = NULL)
#'
#' @inheritParams make_DEgene_dotplot_pseudobulk
#' @param target_regulation_df NULL or a data frame with the columns gene and direction_regulation: indicates whether genes should be divided in up- and downregulated. Default: NULL -- no division.
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' receiver_oi = "Malignant"
#' targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
#' p_target = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, output$grouping_tbl)
#' }
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk_reversed = function(genes_oi, celltype_info, prioritization_tables, celltype_oi, grouping_tbl, groups_oi = NULL, target_regulation_df = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(target_regulation_df)){
    ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####
    keep_tbl = prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, group, receiver, keep_receiver) %>% dplyr::rename(celltype = receiver) %>% dplyr::mutate(keep_receiver = as.logical(keep_receiver))
    
    keep_sender_receiver_values = c(1, 4)
    names(keep_sender_receiver_values) = c(FALSE, TRUE)
    
    plot_data = celltype_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample")) %>% dplyr::inner_join(keep_tbl, by = c("sample","group","celltype"))
    plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% dplyr::ungroup()
    plot_data$gene = factor(plot_data$gene, levels=genes_oi)
    
    plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
    
    if(!is.null(groups_oi)){
      plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)
      
    }
    p1 = plot_data %>%
      ggplot(aes(gene, sample, color = scaled_gene_exprs, size = keep_receiver)) +
      geom_point() +
      facet_grid(group ~ ., scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
        axis.text.y = element_text(size = 9),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.text.x = element_text(size = 9, color = "black", face = "bold"),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("Genes") +ylab("Samples")  +
      scale_size_manual(values = keep_sender_receiver_values)
    max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
    
    p1 = p1 + custom_scale_fill
    
    ####  make the plot that indicates fraction of expression ####
    
    frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
    plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype","group"))
    plot_data$gene = factor(plot_data$gene, levels=genes_oi)
    
    p2 = plot_data %>%
      # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
      # geom_tile(color = "white") +
      ggplot(aes(gene, sample , color = scaled_gene_exprs, size = fraction_sample)) +
      geom_point() +
      facet_grid(group ~ ., scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
        axis.text.y = element_text(size = 9),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.text.x = element_text(size = 9, color = "black", face = "bold"),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      )  +
      labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("Genes") +ylab("Samples") 
    max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
    
    p2 = p2 + custom_scale_fill
    
  } else {
    ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####
    keep_tbl = prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, group, receiver, keep_receiver) %>% dplyr::rename(celltype = receiver) %>% dplyr::mutate(keep_receiver = as.logical(keep_receiver))
    
    keep_sender_receiver_values = c(1, 4)
    names(keep_sender_receiver_values) = c(FALSE, TRUE)
    
    plot_data = celltype_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample")) %>% dplyr::inner_join(keep_tbl, by = c("sample","group","celltype"))
    plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% dplyr::ungroup()
    plot_data$gene = factor(plot_data$gene, levels=genes_oi)
    
    plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi) %>% dplyr::inner_join(target_regulation_df) %>% dplyr::mutate(gene = factor(gene, levels = genes_oi))
    
    if(!is.null(groups_oi)){
      plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)
      
    }
    p1 = plot_data %>%
      ggplot(aes(gene, sample, color = scaled_gene_exprs, size = keep_receiver)) +
      geom_point() +
      facet_grid(group ~ direction_regulation, scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
        axis.text.y = element_text(size = 9),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.y = element_text(size = 9, color = "black", angle = 0),
        strip.text.x = element_text(size = 9, color = "black"),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("Genes") +ylab("Samples") +
      scale_size_manual(values = keep_sender_receiver_values)
    max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
    
    p1 = p1 + custom_scale_fill
    
    ####  make the plot that indicates fraction of expression ####
    
    frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
    plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype"))
    plot_data$gene = factor(plot_data$gene, levels=genes_oi)
    
    p2 = plot_data %>%
      # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
      # geom_tile(color = "white") +
      ggplot(aes(gene, sample , color = scaled_gene_exprs, size = fraction_sample)) +
      geom_point() +
      facet_grid(group ~ direction_regulation, scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
        axis.text.y = element_text(size = 9),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 9, color = "black", face = "bold"),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      )  +
      labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("Genes") +ylab("Samples")
    max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
    
    p2 = p2 + custom_scale_fill
  }
  
  return(list(pseudobulk_plot = p1, singlecell_plot = p2))
  
}
#' @title make_DEgene_dotplot_pseudobulk_batch
#'
#' @description \code{make_DEgene_dotplot_pseudobulk_batch}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes in rows, samples in columns
#' @usage make_DEgene_dotplot_pseudobulk_batch(genes_oi, celltype_info, prioritization_tables, celltype_oi, batch_oi, grouping_tbl, groups_oi = NULL)
#'
#' @param genes_oi Character vector with names of genes to visualize
#' @param celltype_info `celltype_info` or `receiver_info` slot of the output of the `multi_nichenet_analysis` function 
#' @param prioritization_tables `prioritization_tables` slot of the output of the `generate_prioritization_tables` or `multi_nichenet_analysis` function
#' @param celltype_oi Character vector with names of celltype of interest
#' @param grouping_tbl `grouping_tbl` slot of the output of the `multi_nichenet_analysis` function 
#' @param groups_oi Which groups to show? Default: NULL -- will show all groups.
#' @param batch_oi Name of the batch that needs to be visualized for each sample
#'
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
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
#' batches = "batch"
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' receiver_oi = "Malignant"
#' targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(receiver == receiver_oi) %>% pull(gene) %>% unique()
#' p_target = make_DEgene_dotplot_pseudobulk_batch(genes_oi = targets_oi, celltype_info = output$celltype_info, prioritization_tables = output$prioritization_tables, celltype_oi = receiver_oi, batch_oi = batches, output$grouping_tbl)
#'}
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk_batch = function(genes_oi, celltype_info, prioritization_tables, celltype_oi, batch_oi, grouping_tbl, groups_oi = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####
  
  keep_tbl = prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, group, receiver, keep_receiver) %>% dplyr::rename(celltype = receiver) %>% dplyr::mutate(keep_receiver = as.logical(keep_receiver))
  
  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)
  
  plot_data = celltype_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample")) %>% dplyr::inner_join(keep_tbl, by = c("sample","group","celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)
  
  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  
  if(!is.null(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)
    
  }
  
  p1 = plot_data %>%
    ggplot(aes(sample, gene, color = scaled_gene_exprs, size = keep_receiver)) +
    geom_point() +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("") +ylab("Genes")  +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  ####  make the plot that indicates fraction of expression ####
  
  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype","group"))
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)
  
  p2 = plot_data %>%
    # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(sample, gene , color = scaled_gene_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  +
    labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("") +ylab("Genes") 
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p2 = p2 + custom_scale_fill
  
  ### batch plot
  # grouping_tbl = plot_data %>% distinct(sample, group, batch_oi) %>% rename(batch = batch_oi)
  # grouping_tbl_plot = grouping_tbl %>% mutate(batch_ = paste0(" ",batch_oi," "), mock = "")
  grouping_tbl_plot = grouping_tbl %>% mutate(batch_ = paste0(" ",batch_oi," "), mock = "", batch = grouping_tbl[[batch_oi]])
  p_batch = grouping_tbl_plot %>%
    ggplot(aes(sample, mock, fill = batch)) +
    geom_tile(color = "black") +
    facet_grid(.~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.50, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + scale_fill_brewer(palette = "Set2") + xlab("") + labs(fill=batch_oi)
  
  p1 = patchwork::wrap_plots(
    p_batch, p1,
    guides = "collect", nrow = 2,
    heights = c(1, plot_data$gene %>% unique() %>% length())
  )
  p2 = patchwork::wrap_plots(
    p_batch, p2,
    guides = "collect", nrow = 2,
    heights = c(1, plot_data$gene %>% unique() %>% length())
  )
  
  return(list(pseudobulk_plot = p1, singlecell_plot = p2))
  
}
#' @title make_ligand_receptor_violin_plot
#'
#' @description \code{make_ligand_receptor_violin_plot}  Plot combining a violin plot of of the ligand of interest in the sender cell type of interest, and a violin plot of the receptor of interest in the receiver cell type of interest.
#
#' @usage make_ligand_receptor_violin_plot(sce, ligand_oi, receptor_oi, sender_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL)
#'
#' @param sce SingleCellExperiment object
#' @param ligand_oi Character vector of name of the ligand of interest
#' @param receptor_oi Character vector of name of the receptor of interest
#' @param group_oi Character vector of name of the group of interest
#' @param group_id Name of the meta data column that indicates from which group/condition a cell comes from 
#' @param celltype_id Name of the meta data column that indicates the cell type of a cell 
#' @param background_groups Default NULL: all groups in the group_id metadata column will be chosen. But user can fill in a character vector with the names of all gruops of interest.
#' @param sample_id Name of the colData(sce) column in which the id of the sample is defined
#' @param sender_oi Character vector with the names of the sender cell type of interest
#' @param receiver_oi  Character vector with the names of the receiver cell type of interest
#' @param batch_oi Name of a batch of interest based on which the visualization will be split. Default: NA: no batch.
#'
#' @return Plot combining a violin plot of of the ligand of interest in the sender cell type of interest, and a violin plot of the receptor of interest in the receiver cell type of interest.
#'
#' @import ggplot2
#' @import dplyr

#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom generics setdiff
#' @importFrom muscat prepSCE
#' @importFrom RColorBrewer brewer.pal
#' @importFrom generics setdiff
#' @importFrom patchwork wrap_plots
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#' ligand_oi = "DLL1"
#' receptor_oi = "NOTCH3"
#' group_oi = "High"
#' sender_oi = "Malignant"
#' receiver_oi = "myofibroblast"
#' p_violin = make_ligand_receptor_violin_plot(sce = sce, ligand_oi = ligand_oi, receptor_oi = receptor_oi, group_oi = group_oi, group_id = group_id, sender_oi = sender_oi, receiver_oi = receiver_oi, sample_id = sample_id, celltype_id = celltype_id)
#' 
#' }
#'
#' @export
#'
make_ligand_receptor_violin_plot = function(sce, ligand_oi, receptor_oi, sender_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  # ligand plot
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% sender_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[ligand_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[ligand_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    
    p_sender = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the ligand ",ligand_oi, " in sender cell type ", sender_oi, sep = ""))
    
  } else {
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_sender = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi ~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the ligand ",ligand_oi, " in sender cell type ", sender_oi, sep = ""))
    
  }
  
  # receptor plot
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[receptor_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[receptor_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    p_receiver = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the receptor ",receptor_oi, " in receiver cell type ", receiver_oi, sep = ""))
    
  } else {
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_receiver = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the receptor ",receptor_oi, " in receiver cell type ", receiver_oi, sep = ""))
    
  }
  
  return(patchwork::wrap_plots(p_sender, p_receiver, nrow = 2))

}
#' @title make_target_violin_plot
#'
#' @description \code{make_target_violin_plot}  Violin plot of a target gene of interest: per sample, and samples are grouped per group
#' @usage make_target_violin_plot(sce, target_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL)
#'
#' @param target_oi Name of the gene of interest
#' @inheritParams make_ligand_receptor_violin_plot
#'
#' @return ggplot object: Violin plot of a target gene of interest: per sample, and samples are grouped per group
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom generics setdiff
#' @importFrom muscat prepSCE
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#' receiver_oi = "Malignant"
#' group_oi = "High"
#' target_oi = "RAB31"
#' p_violin_target = make_target_violin_plot(sce = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
#' }
#' @export
#'
make_target_violin_plot = function(sce, target_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]

  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[target_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[target_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the target ",target_oi, " in receiver cell type ", receiver_oi, sep = ""))
  } else {
    
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the target ",target_oi, " in receiver cell type ", receiver_oi, sep = ""))
  }
  
  return(p_violin)

}
#' @title make_lr_target_scatter_plot
#'
#' @description \code{make_lr_target_scatter_plot}  Plot Ligand-Receptor pseudobulk expression product values vs pseudobulk expression of correlated target genes supported by prior information.
#' @usage make_lr_target_scatter_plot(prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, receiver_info, grouping_tbl, lr_target_prior_cor_filtered)
#'
#' @param lr_target_prior_cor_filtered Data frame filtered from `lr_target_prior_cor` (= output of `multi_nichenet_analysis` or `lr_target_prior_cor_inference`). Filter should be done to keep onl LR-->Target links that are both supported by prior knowledge and correlation in terms of expression.
#' @inheritParams make_ligand_receptor_violin_plot
#' @inheritParams make_sample_lr_prod_plots
#' @inheritParams make_DEgene_dotplot_pseudobulk_reversed
#' @param receiver_info `celltype_info` or `receiver_info` slot of the output of the `multi_nichenet_analysis` function 
#' 
#' @return ggplot object with plot of LR expression vs target expression
#'
#' @import ggplot2
#' @import dplyr
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#' ligand_oi ="IL24"
#' receptor_oi = "IL22RA1" 
#' sender_oi = "CAF"
#' receiver_oi ="Malignant"
#' lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% filter(scaled_prior_score > 0.50 & (pearson > 0.66 | spearman > 0.66))
#' lr_target_scatter_plot = make_lr_target_scatter_plot(output$prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, output$celltype_info, output$grouping_tbl, lr_target_prior_cor_filtered) 
#' }
#'
#' @export
#'
make_lr_target_scatter_plot = function(prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, receiver_info, grouping_tbl, lr_target_prior_cor_filtered){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  grouping_tbl = grouping_tbl %>% dplyr::inner_join(prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, keep_receiver, keep_sender))
  
  ligand_receptor_pb_prod_df = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1) %>% dplyr::filter(ligand == ligand_oi, receptor == receptor_oi, sender == sender_oi, receiver == receiver_oi) %>% dplyr::distinct(sample, ligand_receptor_pb_prod, id)
  
  targets_oi = ligand_receptor_pb_prod_df %>% dplyr::inner_join(lr_target_prior_cor_filtered, by = "id") %>% dplyr::pull(target) %>% unique()
  
  target_pb_df = receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1) %>% dplyr::filter(gene %in% targets_oi & celltype == receiver_oi) %>% dplyr::distinct(gene, sample, pb_sample) %>% dplyr::rename(target_pb = pb_sample)
  
  ligand_receptor_target_pb_df = grouping_tbl %>% dplyr::inner_join(ligand_receptor_pb_prod_df, by = "sample") %>% dplyr::inner_join(target_pb_df, by = "sample")
  
  p = ligand_receptor_target_pb_df %>% ggplot(aes(ligand_receptor_pb_prod, target_pb)) + 
    geom_point(aes(color = group), size = 2) + 
    geom_smooth(method = "lm",alpha = 0.10, color = "grey50") + 
    facet_grid(.~gene) +
    theme_bw() + 
    scale_color_brewer(palette = "Set2") + ggtitle(paste(ligand_oi,"-",receptor_oi," expression vs target expression between ", sender_oi," and ", receiver_oi, " cells.", sep = "")) +
    xlab("LR pair pseudobulk expression product") + ylab("Target gene pseudobulk expression")
  return(p)
} 
#' @title make_lr_target_prior_cor_heatmap
#'
#' @description \code{make_lr_target_prior_cor_heatmap}  Plot Ligand-Receptor-->Target gene links that are both supported by prior knowledge and have correlation in expression
#' @usage make_lr_target_prior_cor_heatmap(lr_target_prior_cor_filtered, add_grid = TRUE)
#'
#' @inheritParams make_lr_target_scatter_plot
#' @param add_grid add a ggplot-facet grid to easier link LR pairs to target genes. Default: TRUE.
#' 
#' @return ggplot object with plot of LR-->Target links
#'
#' @import ggplot2
#' @import dplyr
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#' lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% filter(scaled_prior_score > 0.50 & (pearson > 0.66 | spearman > 0.66))
#' lr_target_prior_cor_heatmap = make_lr_target_prior_cor_heatmap(lr_target_prior_cor_filtered) 
#' }
#'
#' @export
#'
make_lr_target_prior_cor_heatmap = function(lr_target_prior_cor_filtered, add_grid = TRUE){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  lr_target_prior_cor_filtered = lr_target_prior_cor_filtered %>% dplyr::mutate(pearson = abs(pearson))
  
  if(add_grid == TRUE){
    p_lr_target = lr_target_prior_cor_filtered %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      ggplot(aes(target, lr_interaction, color = prior_score, size = pearson)) + # alternative to scaled_prior_score: scaled_ligand or prior_score
      geom_point() +
      facet_grid(sender_receiver~target, scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        # axis.title = element_blank(),
        # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.175, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x.top = element_text(size = 0, color = "black", face = "bold", angle = 0),
        strip.text.y.right = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background.y = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
        strip.background.x = element_rect(color="darkgrey", fill="whitesmoke", size=0, linetype="solid")
        ) + labs(color = "Prior knowledge score\nof correlated targets")  + labs(size = "Absolute value\nPearson correlation LR and target\npseudobulk expression")
  } else {
    p_lr_target = lr_target_prior_cor_filtered %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      ggplot(aes(target, lr_interaction, color = prior_score, size = pearson)) + # alternative to scaled_prior_score: scaled_ligand or prior_score
      geom_point() +
      facet_grid(sender_receiver~., scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        # axis.title = element_blank(),
        # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.175, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x.top = element_text(size = 0, color = "black", face = "bold", angle = 0),
        strip.text.y.right = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background.y = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
        strip.background.x = element_rect(color="darkgrey", fill="whitesmoke", size=0, linetype="solid")
      ) + labs(color = "Prior knowledge score\nof correlated targets")  + labs(size = "Absolute value\nPearson correlation LR and target\npseudobulk expression")
  }
  
  # custom_scale_color = scale_color_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "BuPu")) ,values = c(0, 0.25, 0.40, 0.50, 0.65, 0.775, 0.90, 1),  limits = c(0, 1))
  custom_scale_color = scale_color_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.025, 0.075, 0.25, 0.40, 0.55, 1),  limits = c(0, max(lr_target_prior_cor_filtered$prior_score)))
  
  # custom_scale_color = scale_color_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "BuPu")) ,values = c(0, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1),  limits = c(0, 1))
  
  p_lr_target = p_lr_target + custom_scale_color + xlab("Correlated target genes\nsupported by prior knowledge") + ylab("Prioritzed LR pairs")
  p_lr_target
}
#' @title make_lr_target_correlation_plot
#'
#' @description \code{make_lr_target_correlation_plot}  Plot Ligand-Receptor expression, Ligand-Receptor-->Target gene links that are both supported by prior knowledge and have correlation in expression, and Target expression
#' @usage make_lr_target_correlation_plot(prioritization_tables, prioritized_tbl_oi, lr_target_prior_cor_filtered, grouping_tbl, receiver_info, receiver_oi, plot_legend = TRUE, heights = NULL, widths = NULL)
#'
#' @inheritParams make_lr_target_scatter_plot
#' @inheritParams make_sample_lr_prod_plots
#' @inheritParams make_ligand_activity_target_plot
#' @param receiver_info `celltype_info` or `receiver_info` slot of the output of the `multi_nichenet_analysis` function 
#' 
#' @return ggplot object with a combined plot of LR expression vs target expression
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom generics intersect
#' @importFrom patchwork wrap_plots
#' @importFrom ggpubr as_ggplot
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      )
#' group_oi = "High"
#' receiver_oi = "Malignant"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% distinct(id, ligand, receptor, sender, receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity_scaled, fraction_expressing_ligand_receptor,  prioritization_score) %>% filter(fraction_expressing_ligand_receptor > 0 & ligand_receptor_lfc_avg > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(250, prioritization_score) 
#' lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% filter(scaled_prior_score > 0.50 & (pearson > 0.66 | spearman > 0.66))
#' prioritized_tbl_oi = prioritized_tbl_oi %>% filter(id %in% lr_target_prior_cor_filtered$id)
#' prioritized_tbl_oi = prioritized_tbl_oi %>% group_by(ligand, sender, group) %>% top_n(2, prioritization_score)
#' lr_target_correlation_plot = make_lr_target_correlation_plot(output$prioritization_tables, prioritized_tbl_oi, lr_target_prior_cor_filtered, output$grouping_tbl, output$celltype_info, receiver_oi) 
#' }
#'
#' @export
#'
make_lr_target_correlation_plot = function(prioritization_tables, prioritized_tbl_oi, lr_target_prior_cor_filtered, grouping_tbl, receiver_info, receiver_oi, plot_legend = TRUE, heights = NULL, widths = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  lr_target_prior_cor_filtered = lr_target_prior_cor_filtered %>% dplyr::mutate(pearson = abs(pearson), direction_regulation = factor(direction_regulation, levels = c("up","down")))
  
  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - ")) %>%  dplyr::arrange(sender)
  group_data = prioritization_tables$group_prioritization_table_source %>% dplyr::distinct(id, group) %>% dplyr::filter(id %in% sample_data$id)
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
  
  sample_data = sample_data %>% dplyr::mutate(sample = factor(sample)) %>% dplyr::mutate(sample = factor(sample, levels = rev(levels(sample)))) # reverse - to have the same order as used for the targets
  
  ##### LR-->Target plot
  
  lr_target_data = lr_target_prior_cor_filtered %>% dplyr::inner_join(sample_data, by = c("sender", "receiver", "ligand", "receptor", "id", "group")) %>% dplyr::select(id, lr_interaction, sender, receiver, sender_receiver, group, target, prior_score, pearson, direction_regulation) %>% dplyr::distinct() %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "))
    lr_target_data = lr_target_data %>% dplyr::mutate(target = factor(target))
    order_targets = lr_target_data$target %>% levels()
    # lr_target_data_filled = lr_target_data %>% dplyr::distinct(id, target, scaled_prior_score) %>% tidyr::spread(target, scaled_prior_score, fill = 0) %>% tidyr::gather(target, score, -id) # alternative to scaled_prior_score: scaled_ligand or prior_score
    # lr_target_data = lr_target_data_filled %>% dplyr::inner_join(lr_target_data %>% dplyr::select(id, lr_interaction, sender_receiver, group, pearson), by = "id")
    lr_target_data = lr_target_data %>% dplyr::mutate(lr_interaction = factor(lr_interaction))
    order_lr_interaction = lr_target_data$lr_interaction %>% levels()
    
    p2 = lr_target_data %>%
      ggplot(aes(target, lr_interaction, color = prior_score, size = pearson)) +
      geom_point() +
      facet_grid(sender_receiver~direction_regulation, scales = "free", space = "free") +
      scale_x_discrete(position = "top") +
      theme_light() +
      theme(
        axis.ticks = element_blank(),
        # axis.title = element_blank(),
        # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 9),
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.50, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x.top = element_text(size = 9, color = "black",  angle = 0),
        strip.text.y.right = element_text(size = 9, color = "black",  angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + labs(fill = "Prior knowledge score\nof correlated targets", size= "Absolute value\nPearson correlation LR and target\npseudobulk expression") 
    custom_scale_color = scale_color_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.0125, 0.0040, 0.175, 0.35, 0.50, 1),  limits = c(0, max(lr_target_prior_cor_filtered$prior_score)))
    # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "BuPu")) ,values = c(0, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1),  limits = c(0, 1))
    
    p_lr_target = p2 + custom_scale_color + xlab("Correlated target genes\nsupported by prior knowledge") + ylab("")
  
  ##### Target expression plot
  # Target expression
  groups_oi = group_data %>% dplyr::pull(group) %>% unique()
  target_regulation_df = lr_target_prior_cor_filtered %>% dplyr::distinct(target,  direction_regulation)  %>% dplyr::rename(gene = target)
  p_target_exprs = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = order_targets, celltype_info = receiver_info, prioritization_tables = prioritization_tables, celltype_oi = receiver_oi, grouping_tbl = grouping_tbl, groups_oi = groups_oi, target_regulation_df = target_regulation_df) %>% .$pseudobulk_plot + xlab("") + ylab("")
  
  ##### LR product plot

  sample_data = sample_data %>% dplyr::filter(lr_interaction %in% order_lr_interaction) %>% dplyr::mutate(lr_interaction = factor(lr_interaction))
  
  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      # axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\n pseudobulk product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(sample_data$scaled_LR_pb_prod ) %>% max()
  custom_scale_color = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_lr_exprs = p1 + custom_scale_color + xlab("Sample") + ylab("Prioritized LR pairs")

  
  #### now combine all three plots
  # Combine the plots
  n_targets = length(order_targets)
  n_ligands = length(lr_target_data$id %>% unique())
  n_samples = grouping_tbl %>% dplyr::filter(group %in% groups_oi) %>% dplyr::pull(sample) %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_lr_exprs)), ggpubr::as_ggplot(ggpubr::get_legend(p_lr_target)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_target_exprs)))
  
  if(is.null(heights)){
    heights = c(n_ligands, n_samples)
  }
  if(is.null(widths)){
    widths = c(n_samples, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "AB
               #C"
    combined_plot = patchwork::wrap_plots(A = p_lr_exprs + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          B = p_lr_target + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_target_exprs + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "AB
               LC"
    
    combined_plot = patchwork::wrap_plots(A = p_lr_exprs + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          
                                          B = p_lr_target + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_target_exprs + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
  
  
}  
#' @title make_ggraph_ligand_target_links
#'
#' @description \code{make_ggraph_ligand_target_links} Make a network showing the gene regulatory links between ligands from sender cell types to their induced ligands/receptors in receiver cell types. Lins are only drawn if the ligand/receptor in the receiver is a potential downstream target of the ligand based on prior knowledge and sufficient correlation in expression across the different samples.
#' @usage make_ggraph_ligand_target_links(lr_target_prior_cor_filtered, prioritized_tbl_oi, colors)
#'
#' @inheritParams make_lr_target_scatter_plot
#' @inheritParams make_sample_lr_prod_plots
#' @param colors Named vector of colors associated to each sender cell type. Vector = color, names = sender names. 
#' 
#' @return ggplot object with plot of LR-->Target links, the graph object itself, and the prioritized LR links
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr set_rownames
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @import ggraph
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% filter(scaled_prior_score > 0.50 & (pearson > 0.66 | spearman > 0.66))
#'  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% distinct(id, ligand, receptor, sender, receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity_scaled, fraction_expressing_ligand_receptor,  prioritization_score) %>% filter(fraction_expressing_ligand_receptor > 0 & ligand_receptor_lfc_avg > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(250, prioritization_score)
#'  prioritized_tbl_oi = prioritized_tbl_oi %>% filter(id %in% lr_target_prior_cor_filtered$id)
#'  prioritized_tbl_oi = prioritized_tbl_oi %>% group_by(ligand, sender, group) %>% top_n(2, prioritization_score)
#' graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered, prioritized_tbl_oi = prioritized_tbl_oi, colors = c("blue","red"))
#' graph_plot$plot
#' }
#'
#' @export
#'
make_ggraph_ligand_target_links = function(lr_target_prior_cor_filtered, prioritized_tbl_oi, colors){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  
  lr_target_prior_cor_filtered = lr_target_prior_cor_filtered %>% dplyr::inner_join(prioritized_tbl_oi, by = c("sender", "receiver", "ligand", "receptor", "id", "group"))
  
  source_df_lr = prioritized_tbl_oi %>% dplyr::mutate(celltype_ligand = paste(sender, ligand, sep = "_"), celltype_receptor = paste(receiver, receptor, sep = "_")) %>% dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, ligand, receptor) 
  source_df_lrt = lr_target_prior_cor_filtered %>% dplyr::mutate(celltype_ligand = paste(sender, ligand, sep = "_"), celltype_target = paste(receiver, target, sep = "_"), celltype_receptor = paste(receiver, receptor, sep = "_")) %>% dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, celltype_target, ligand, target, receptor, direction_regulation) 
  
  lr_gr_network = dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & !celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand"),
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_receptor  & !celltype_target %in% source_df_lr$celltype_ligand) %>% dplyr::mutate(type_target = "receptor")
  ) %>% dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand/receptor")
  )
  
  ligand_receptor_network = lr_gr_network %>% dplyr::filter(celltype_receptor %in% celltype_target) %>% dplyr::select(celltype_ligand, celltype_receptor, direction_regulation, group) %>% dplyr::distinct() %>% dplyr::rename(sender = celltype_ligand, receiver = celltype_receptor) %>% dplyr::mutate(type = "Ligand-Receptor", weight = 1)
  
  ligand_target_network = lr_gr_network %>% dplyr::select(celltype_ligand, celltype_target, direction_regulation, group) %>% dplyr::distinct() %>% dplyr::rename(sender = celltype_ligand, receiver = celltype_target) %>% dplyr::mutate(type = "Ligand-Target", weight = 1)
  
  links = ligand_target_network %>% dplyr::bind_rows(ligand_receptor_network) 
  nodes = lr_gr_network %>% dplyr::select(celltype_ligand, sender, ligand) %>% dplyr::rename(celltype = sender, node = celltype_ligand, gene = ligand) %>% dplyr::mutate(type_gene = "ligand") %>% dplyr::bind_rows(
    lr_gr_network %>% dplyr::select(celltype_receptor, receiver, receptor) %>% dplyr::rename(celltype = receiver, node = celltype_receptor, gene = receptor) %>% dplyr::mutate(type_gene = "receptor")
  ) %>% dplyr::bind_rows(
    lr_gr_network %>% dplyr::select(celltype_target, receiver, target, type_target) %>% dplyr::rename(celltype = receiver, node = celltype_target, gene = target, type_gene = type_target)
  ) %>% dplyr::distinct() %>% dplyr::filter(node %in% c(links$sender, links$receiver))
  
  double_nodes =  nodes %>% dplyr::group_by(node) %>% dplyr::count() %>% dplyr::filter(n > 1) %>% pull(node)
  nodes = dplyr::bind_rows(
    nodes %>% dplyr::filter(node %in% double_nodes) %>% dplyr::mutate(type_gene = "ligand/receptor") ,
    nodes %>% dplyr::filter(!node %in% double_nodes)
  ) %>% dplyr:: distinct()
  
  nodes = nodes %>% data.frame() %>% magrittr::set_rownames(nodes$node)
  
  colors_regulation = NULL
  colors_regulation["up"] = "indianred1"
  colors_regulation["down"] = "steelblue2"
  
  # create the network object
  network = igraph::graph_from_data_frame(d=links %>% dplyr::filter(type == "Ligand-Target"), vertices = nodes, directed=T)
  graph = tidygraph::as_tbl_graph(network) 
  set.seed(1919)
  plot =  ggraph::ggraph(graph, layout = 'dh') + 
    ggraph::geom_edge_fan(aes(color = direction_regulation), edge_width = 1, arrow = arrow(length = unit(3, 'mm')), end_cap = ggraph::circle(5.5, 'mm'), start_cap = ggraph::circle(3, 'mm')) + 
    # ggraph::geom_edge_loop(aes(color = direction_regulation), edge_width = 1, alpha = 0.70)  + 
    ggraph::geom_node_label(aes(label = gene, fill = celltype), fontface = "bold", size = 3.5, nudge_x = 0, nudge_y = 0, color = "whitesmoke") +
    ggraph::theme_graph(foreground = 'black', fg_text_colour = 'white', base_family = 'Helvetica') + facet_grid(. ~group)  + ggraph::scale_edge_color_manual(values = colors_regulation) + scale_fill_manual(values = colors)
  
  return(list(plot = plot, graph = graph, source_df_lr = source_df_lr, source_df_lt = links %>% dplyr::filter(type == "Ligand-Target"), nodes_df = nodes))
}
#' @title visualize_network
#'
#' @description \code{visualize_network} Visualize a network showing the gene regulatory links between ligands from sender cell types to their induced ligands/receptors in receiver cell types. Links are only drawn if the ligand/receptor in the receiver is a potential downstream target of the ligand (based on prior knowledge, and optionally with sufficient correlation in expression across the different samples).
#' @usage visualize_network(network, colors)
#'
#' @param network Output of `infer_intercellular_regulatory_network`
#' @param colors Named vector of colors associated to each sender cell type. Vector = color, names = sender names. 
#' 
#' @return ggplot object with plot of LR-->Target links, together with the tidygraph and igraph objects themselves
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr set_rownames
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @import ggraph
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' lr_target_prior_cor_filtered = output$lr_target_prior_cor %>% filter(scaled_prior_score > 0.50 & (pearson > 0.66 | spearman > 0.66))
#'  prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% distinct(id, ligand, receptor, sender, receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity_scaled, fraction_expressing_ligand_receptor,  prioritization_score) %>% filter(fraction_expressing_ligand_receptor > 0 & ligand_receptor_lfc_avg > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(250, prioritization_score)
#'  prioritized_tbl_oi = prioritized_tbl_oi %>% filter(id %in% lr_target_prior_cor_filtered$id)
#'  prioritized_tbl_oi = prioritized_tbl_oi %>% group_by(ligand, sender, group) %>% top_n(2, prioritization_score)
#'  lr_target_df = lr_target_prior_cor_filtered  %>% distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 
#'  network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi)
#' graph_plot = visualize_network(network, colors = c("purple","orange"))
#' graph_plot$plot
#' }
#'
#' @export
#'
visualize_network = function(network, colors){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  
  nodes = network$nodes %>% data.frame() %>% magrittr::set_rownames(network$nodes$node)
  
  colors_regulation = NULL
  colors_regulation["up"] = "indianred1"
  colors_regulation["down"] = "steelblue2"
  
  # create the network object
  network_igraph = igraph::graph_from_data_frame(d=network$links %>% dplyr::filter(type == "Ligand-Target"), vertices = nodes, directed=T)
  network_tidygraph = tidygraph::as_tbl_graph(network_igraph) 
  set.seed(1919)
  plot =  ggraph::ggraph(network_tidygraph, layout = 'dh') + 
    ggraph::geom_edge_fan(aes(color = direction_regulation), edge_width = 1, arrow = arrow(length = unit(3, 'mm')), end_cap = ggraph::circle(5.5, 'mm'), start_cap = ggraph::circle(3, 'mm')) + 
    #ggraph::geom_edge_loop(aes(color = direction_regulation), edge_width = 1, alpha = 0.70)  + 
    ggraph::geom_node_label(aes(label = gene, fill = celltype), fontface = "bold", size = 3.5, nudge_x = 0, nudge_y = 0, color = "whitesmoke") +
    ggraph::theme_graph(foreground = 'black', fg_text_colour = 'white', base_family = 'Helvetica') + facet_grid(. ~group)  + ggraph::scale_edge_color_manual(values = colors_regulation) + scale_fill_manual(values = colors)
  
  
  return(list(plot = plot, network_igraph = network_igraph, network_tidygraph = network_tidygraph))
}
#' @title make_ggraph_signaling_path
#'
#' @description \code{make_ggraph_signaling_path} Visualize the Ligand-Receptor to target signaling paths
#' @usage make_ggraph_signaling_path(signaling_graph_list, colors, ligands_all, receptors_all, targets_all)
#'
#' @param signaling_graph_list Output of `nichenetr::get_ligand_signaling_path_with_receptor`
#' @param colors Named vector of colors associated to each node type: Example: colors <- c("ligand" = "indianred2", "receptor" = "orange", "target" = "steelblue2", "mediator" = "grey25"). 
#' @param ligands_all Name of the ligand(s)
#' @param receptors_all Name of the receptor(s)
#' @param targets_all Name of the target(s)
#' 
#' @return ggraph and tidygraph objec of signaling paths between predefined LR-->Target links
#'
#' @import ggplot2
#' @import dplyr
#' @import ggraph
#' @importFrom magrittr set_rownames set_names
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
#' ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
#' gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))

#' ligands_all = "COL1A1" # this can be a list of multiple ligands if required
#' receptors_all = "ITGB1"
#' targets_all = c("S100A1","SERPINE1")

#' active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, receptors_all = receptors_all, targets_all = targets_all, weighted_networks = weighted_networks, top_n_regulators = 2)
#' data_source_network = nichenetr::infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)

#' active_signaling_network_min_max = active_signaling_network
#' active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
#' active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
#' colors = c("ligand" = "indianred2", "receptor" = "orange", "target" = "steelblue2", "mediator" = "grey75")
#' ggraph_signaling_path = make_ggraph_signaling_path(active_signaling_network_min_max, colors)#' colors = c("ligand" = "indianred2", "receptor" = "orange", "target" = "steelblue2", "mediator" = "grey25")
#' 
#' }
#'
#' @export
#'
make_ggraph_signaling_path = function(signaling_graph_list, colors, ligands_all, receptors_all, targets_all){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  edge_colors = c("signaling" = "indianred1", "gene regulatory" = "steelblue2")
  edge_df = dplyr::bind_rows(signaling_graph_list$sig %>% dplyr::mutate(interaction_type = "signaling"), 
                      signaling_graph_list$gr %>% dplyr::mutate(interaction_type = "gene regulatory"))
  
  nodes = union(edge_df$from, edge_df$to)
  node2type = rep("mediator", times = length(nodes)) %>% magrittr::set_names(nodes)
  node2type[which(nodes %in% ligands_all)] = "ligand"
  node2type[which(nodes %in% receptors_all)] = "receptor"
  node2type[which(nodes %in% targets_all)] = "target"
  nodes_df = tibble::tibble(node = nodes) %>% dplyr::mutate(node_type = node2type[node])
  nodes_df = data.frame(nodes_df) %>% magrittr::set_rownames(nodes_df$node)
  network = igraph::graph_from_data_frame(d = edge_df, vertices = nodes_df, directed=T)
  graph = tidygraph::as_tbl_graph(network) 
  set.seed(1919)
  plot = suppressWarnings(ggraph::ggraph(graph,layout = "nicely") + 
    ggraph::geom_edge_link(aes(colour = interaction_type, width = weight), arrow = arrow(length = unit(4, 'mm')), end_cap = ggraph::circle(5.5, 'mm'), start_cap = ggraph::circle(3, 'mm')) +
    ggraph::geom_node_label(aes(label = name, color = node_type), fontface = "bold", size = 4.5, nudge_x = 0, nudge_y = 0, family= "Serif") +
    ggraph::theme_graph(foreground = 'black', fg_text_colour = 'white') + scale_color_manual(values = colors)) + scale_edge_color_manual(values = edge_colors) + scale_edge_width_continuous(range = c(0.5, 1.5))
  return(list(plot = plot, graph = graph))
  
}
#' @title make_ligand_activity_target_plot
#'
#' @description \code{make_ligand_activity_target_plot}  Summary plot showing the activity of prioritized ligands acting on a receiver cell type of interest, together with the predicted target genes and their sample-by-sample expression
#' @usage make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi, prioritization_tables, ligand_activities_targets_DEgenes, contrast_tbl, grouping_tbl, receiver_info, ligand_target_matrix, groups_oi = NULL, plot_legend = TRUE, heights = NULL, widths = NULL)
#'
#' @param receiver_oi Character vector of receiver cell type of interest
#' @param group_oi Character vector: name of the group of interest
#' @param ligand_activities_targets_DEgenes Sublist in the output of `multi_nichenet_analysis`
#' @param plot_legend if TRUE (default): show legend on the same figure, if FALSE (recommended): show legend in separate figure
#' @param heights Vector of 2 elements: height of the ligand-activity-target panel, height of the target expression panel. Default NULL: automatically defined based on nr of ligands and samples. If manual change: example format: c(1,1) 
#' @param widths Vector of 3 elements: Width of the scaled ligand activity panel, width of the ligand activity panel, width of the ligand-target heatmap panel. Default NULL: automatically defined based on nr of target genes and group-receiver combinations. If manual change: example format: c(1,1,10) 
#' @inheritParams make_DEgene_dotplot_pseudobulk_reversed
#' @inheritParams make_ligand_activity_plots
#' @inheritParams make_DEgene_dotplot_pseudobulk
#' @inheritParams make_sample_lr_prod_plots
#' @inheritParams multi_nichenet_analysis
#' @param receiver_info `celltype_info` or `receiver_info` slot of the output of the `multi_nichenet_analysis` function 
#' 
#' @return Summary plot showing the activity of prioritized ligands acting on a receiver cell type of interest, together with the predicted target genes and their sample-by-sample expression
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nichenetr prepare_ligand_target_visualization make_heatmap_ggplot
#' @importFrom generics intersect
#' @importFrom patchwork wrap_plots
#' @importFrom ggpubr as_ggplot
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' group_oi = "High"
#' receiver_oi = "Malignant"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi & receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% top_n(25, activity_scaled) %>% arrange(-activity_scaled)
#' combined_plot = make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi, output$prioritization_tables, output$ligand_activities_targets_DEgenes, contrast_tbl, output$grouping_tbl, output$celltype_info,ligand_target_matrix, plot_legend = FALSE)
#' 
#' }
#'
#' @export
#'
make_ligand_activity_target_plot = function(group_oi, receiver_oi, prioritized_tbl_oi, prioritization_tables, ligand_activities_targets_DEgenes, contrast_tbl, grouping_tbl, receiver_info, ligand_target_matrix, groups_oi = NULL, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(groups_oi)){
    groups_oi = contrast_tbl %>% dplyr::pull(group) %>% unique() 
  }
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
  
  # Ligand-Target heatmap
  active_ligand_target_links_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::select(ligand, target, ligand_target_weight, direction_regulation) %>% dplyr::rename(weight = ligand_target_weight )
  
  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.2
  }
  
  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)
  
  order_ligands_ = generics::intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets_ = active_ligand_target_links_df$target %>% unique() %>% generics::intersect(rownames(active_ligand_target_links))
  
  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }
  
  vis_ligand_target_df = vis_ligand_target %>% data.frame() %>% tibble::rownames_to_column("ligand") %>% tidyr::gather("target","score", -ligand) %>% tibble::as_tibble() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands))  %>% dplyr::inner_join(active_ligand_target_links_df %>% distinct(target, direction_regulation)) %>% dplyr::mutate(target = factor(target, levels = order_targets))
  
  p_ligand_target_network = vis_ligand_target_df %>% ggplot(aes(target,ligand,fill = score)) + 
    geom_tile(color = "whitesmoke", size = 0.5) + 
    facet_grid(.~direction_regulation, scales = "free", space = "free") +
    scale_fill_gradient2(low = "white", mid = "purple", high = "darkred", midpoint = 0.14) + theme_light() +
    scale_x_discrete(position = "top") + 
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.50, "lines"),
      strip.text.x = element_text(size = 9, color = "black"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Regulatory Potential") + xlab("Predicted target genes") + ylab("Prioritized ligands")
  
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", "plum1", "orchid2","orchid4","violetred"),values = c(0, 0.05, 0.50, 0.80, 1),  limits = c(0, max(ligand_target_matrix)))
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.025, 0.075, 0.25, 0.40, 0.55, 1),  limits = c(0, max(ligand_target_matrix)))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.04, 0.12, 0.30, 0.40, 0.55, 1),  limits = c(0, max(ligand_target_matrix)))
  
  p_ligand_target_network = p_ligand_target_network + custom_scale_fill
  
  # Ligand-Activity-Scaled -----
  ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% order_ligands_ & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(group %in% groups_oi) %>% dplyr::select(ligand, group, direction_regulation, activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands)) 
  
  p_ligand_activity_scaled = ligand_activity_df %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(direction_regulation , ligand, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 10, color = "black"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver") + ylab("Prioritized ligands") + xlab("Scaled ligand activity")
  max_activity = abs(ligand_activity_df$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))
  p_ligand_activity_scaled = p_ligand_activity_scaled + custom_scale_fill
  
  # Ligand-Activity -----
  ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% order_ligands_ & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(group %in% groups_oi) %>% dplyr::select(ligand, group, direction_regulation, activity) %>% dplyr::distinct() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands)) 
  
  p_ligand_activity = ligand_activity_df %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(direction_regulation , ligand, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.50, "lines"),
      panel.spacing.y = unit(0.50, "lines"),
      strip.text.x = element_text(size = 10, color = "black"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand Activity\nin Receiver") + ylab("Prioritized ligands") + xlab("Ligand activity")
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  p_ligand_activity = p_ligand_activity + custom_scale_fill
  
  # Target expression
  target_regulation_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::distinct(target,  direction_regulation)  %>% dplyr::rename(gene = target)
  
  p_targets = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = order_targets_, celltype_info = receiver_info, prioritization_tables = prioritization_tables, celltype_oi = receiver_oi, grouping_tbl = grouping_tbl, groups_oi = groups_oi, target_regulation_df = target_regulation_df)
  
  # Combine the plots -----
  n_groups = ligand_activity_df$group %>% unique() %>% length()
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_samples = grouping_tbl %>% dplyr::filter(group %in% groups_oi) %>% dplyr::pull(sample) %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets$pseudobulk_plot)))
  
  if(is.null(heights)){
    heights = c(n_ligands + 3, n_samples)
  }
  if(is.null(widths)){
    widths = c(n_groups*2 + 0.75, n_groups*2, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "AaB
               ##C"
    combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_activity + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets$pseudobulk_plot + theme(legend.position = "none") + xlab(""),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "AaB
               L#C"
    
    combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_activity + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets$pseudobulk_plot + theme(legend.position = "none") + xlab(""),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
}  
#' @title make_circos_group_comparison
#'
#' @description \code{make_circos_group_comparison}  Make a circos plot with top prioritized ligand-receptor interactions for each group of interest. In each circos, all the possible LR pairs will be shown, but arrows will only be drawn between the ones belonging to the group of interest.
#' @usage make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#'
#' @inheritParams make_sample_lr_prod_plots
#' @param colors_sender Named vector of colors associated to each sender cell type. Vector = color, names = sender names. If sender and receiver cell types are the same, recommended that this vector is the same as `colors_receiver`.
#' @param colors_receiver Named vector of colors associated to each receiver cell type. Vector = color, names = sender names. Vector = color, names = sender names. If sender and receiver cell types are the same, recommended that this vector is the same as `colors_receiver`.
#'
#' @return a list with a circos plot for each group of interest, and a legend showing the color corresponding to each sender/receiver cell type.
#'
#' @import ggplot2
#' @import dplyr
#' @import circlize
#' @importFrom tibble tibble
#' @importFrom generics intersect
#' @importFrom ComplexHeatmap draw Legend
#' @importFrom magrittr set_names
#' @importFrom grDevices recordPlot
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
# prioritized_tbl_oi_prep = output$prioritization_tables$group_prioritization_tbl %>% 
#   distinct(id, sender, receiver, ligand, receptor, group, prioritization_score, ligand_receptor_lfc_avg, fraction_expressing_ligand_receptor) %>% 
#   filter(ligand_receptor_lfc_avg > 0 & fraction_expressing_ligand_receptor > 0) %>% top_n(30, prioritization_score) 
# prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% 
#   filter(id %in% prioritized_tbl_oi_prep$id) %>% 
#   distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_prep)
# prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
# senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique())
# colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
# colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
# circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#'}
#' @export
#'
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)

    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
      gap = width_different_cell
      return(c(sector,gap))
    }) %>% unlist()
    sender_gaps = sender_gaps[-length(sender_gaps)]
    
    receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
      gap = width_different_cell
      return(c(sector,gap))
    }) %>% unlist()
    receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    
    gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    
    if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
      warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#' @title make_circos_one_group
#'
#' @description \code{make_circos_one_group}  Make a circos plot with top prioritized ligand-receptor interactions for one group of interest. 
#' @usage make_circos_one_group(prioritized_tbl_oi, colors_sender, colors_receiver)
#'
#' @inheritParams make_circos_group_comparison
#'
#' @return a list with a circos plot for one group of interest, and a legend showing the color corresponding to each sender/receiver cell type.
#'
#' @import ggplot2
#' @import dplyr
#' @import circlize
#' @importFrom tibble tibble
#' @importFrom generics intersect
#' @importFrom ComplexHeatmap draw Legend
#' @importFrom magrittr set_names
#' @importFrom grDevices recordPlot
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
#' output = multi_nichenet_analysis(
#'      sce = sce, 
#'      celltype_id = celltype_id, 
#'      sample_id = sample_id, 
#'      group_id = group_id,
#'      batches = batches,
#'      lr_network = lr_network, 
#'      ligand_target_matrix = ligand_target_matrix, 
#'      contrasts_oi = contrasts_oi, 
#'      contrast_tbl = contrast_tbl
#'      
#'      )
#' 
# group_oi = "High"
# prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(ligand_receptor_lfc_avg > 0 & fraction_expressing_ligand_receptor > 0 & group == group_oi) %>% top_n(25, prioritization_score)
# senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique())
# colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
# colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
# circos_list = make_circos_one_group(prioritized_tbl_oi, colors_sender, colors_receiver)
#' }
#' @export
#'
make_circos_one_group = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
      gap = width_different_cell
      return(c(sector,gap))
    }) %>% unlist()
    sender_gaps = sender_gaps[-length(sender_gaps)]
    
    receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
      gap = width_different_cell
      return(c(sector,gap))
    }) %>% unlist()
    receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    
    gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    
    # print(length(gaps))
    # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
    if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
      warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    }
    
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 # transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}
#' @title compare_normal_emp_pvals
#'
#' @description \code{compare_normal_emp_pvals}  Compare nr and rank of DE genes between normal p-values and empirical p-values
#' @usage compare_normal_emp_pvals(DE_info, DE_info_emp, adj_pval = FALSE)
#'
#' @param DE_info Output of `get_DE_info`
#' @param DE_info_emp Output of `get_empirical_pvals`
#' @param adj_pval Should the adjusted p-values be compared (TRUE) or the non-adjusted ones (FALSE)? Defautl: FALSE
#'
#' @return a list òf plots for each celltype-contrast pair: an upset plot and line plot are drawn.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom UpSetR upset
#' @importFrom magrittr set_rownames
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
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' DE_info = get_DE_info(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    covariates = covariates,
#'    contrasts = contrasts_oi)
#' DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
#' comparison_plots = compare_normal_emp_pvals(DE_info, DE_info_emp)
#' }
#' @export
#'
compare_normal_emp_pvals = function(DE_info, DE_info_emp, adj_pval = FALSE){
  
  requireNamespace("dplyr")
  requireNamespace(("ggplot2"))
  
  comparison_plots = DE_info$celltype_de$de_output_tidy$cluster_id %>% unique() %>% lapply(function(celltype_oi, adjusted = FALSE){
    if(adjusted == TRUE){
      de_genes_normal = DE_info$celltype_de$de_output_tidy %>% dplyr::filter(cluster_id == celltype_oi) %>% dplyr::filter(p_adj <= 0.05) %>% dplyr::pull(gene) %>% unique()
      de_genes_emp = DE_info_emp$de_output_tidy_emp %>% dplyr::filter(cluster_id == celltype_oi) %>% dplyr::filter(p_adj_emp <= 0.05) %>% dplyr::pull(gene) %>% unique()
      
    } else {
      de_genes_normal = DE_info$celltype_de$de_output_tidy %>% dplyr::filter(cluster_id == celltype_oi) %>% dplyr::filter(p_val <= 0.05) %>% dplyr::pull(gene) %>% unique()
      de_genes_emp = DE_info_emp$de_output_tidy_emp %>% dplyr::filter(cluster_id == celltype_oi) %>% dplyr::filter(p_emp <= 0.05) %>% dplyr::pull(gene) %>% unique()
      
    }
    
    upset_df = tibble::tibble(gene = union(de_genes_normal, de_genes_emp), normal = as.double(gene %in% de_genes_normal), empirical = as.double(gene %in% de_genes_emp)) %>% data.frame() %>% magrittr::set_rownames(.$gene) %>% dplyr::select(-gene)
    colnames(upset_df) = paste(colnames(upset_df), celltype_oi, sep = "-")
    p_upset = UpSetR::upset(upset_df, sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on") 
    
    p_ranking = DE_info_emp$de_output_tidy_emp %>% dplyr::filter(gene %in% union(de_genes_normal, de_genes_emp) & cluster_id == celltype_oi) %>% dplyr::group_by(cluster_id, contrast) %>% mutate(normal = rank(p_val), empirical = rank(p_emp)) %>% dplyr::filter(normal != empirical) %>% dplyr::mutate(empirical_lower = empirical < normal) %>% tidyr::gather(rank_type, rank, normal:empirical) %>% dplyr::select(cluster_id, contrast,  gene, rank_type, rank, empirical_lower)  %>% 
      ggplot(aes(rank_type, rank, group = gene, color = empirical_lower)) + geom_line(aes(group = gene)) + facet_grid(cluster_id ~ contrast) + theme_bw()
    
    return(list(p_upset, p_ranking))
    
  }, adjusted = adj_pval) 
  return(comparison_plots)
}