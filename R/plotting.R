#' @title make_sample_lr_prod_plots
#'
#' @description \code{make_sample_lr_prod_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups
#' @usage make_sample_lr_prod_plots(prioritization_tables, prioritized_tbl_oi)
#'
#' @param prioritization_tables XXX
#' @param prioritized_tbl_oi XXX
#'
#' @return Ligand-Receptor Expression Product Dotplot/Heatmap
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_sample_lr_prod_plots = function(prioritization_tables, prioritized_tbl_oi){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  filtered_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>%  dplyr::arrange(sender) %>% dplyr::group_by(sender) %>%  dplyr::arrange(receiver)
  
  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4.25)
  names(keep_sender_receiver_values) = levels(filtered_data$keep_sender_receiver)
  
  
  p1 = filtered_data %>%
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
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(filtered_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill

  return(p1)

}

#' @title make_sample_lr_prod_activity_plots
#'
#' @description \code{make_sample_lr_prod_activity_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups. In addition, show the NicheNet ligand activities in each receiver-celltype combination.
#' @usage make_sample_lr_prod_activity_plots(prioritization_tables, prioritized_tbl_oi, widths = NULL)
#'
#' @param prioritization_tables XXX
#' @param prioritized_tbl_oi XXX
#' @param widths Vector of 3 elements: Width of the LR exprs product panel,  width of the scaled ligand activity panel, width of the ligand activity panel. Default NULL: automatically defined based on nr of samples vs nr of group-receiver combinations. If manual change: example format: c(5,1,1) 
#'
#' @return Ligand-Receptor Expression Product & Ligand Activities Dotplot/Heatmap 
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_sample_lr_prod_activity_plots = function(prioritization_tables, prioritized_tbl_oi, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>%  dplyr::arrange(sender) %>% dplyr::group_by(sender) %>%  dplyr::arrange(receiver)

  group_data = prioritization_tables$group_prioritization_tbl %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity, activity_scaled, fraction_ligand_group, prioritization_score, scaled_avg_exprs_ligand) %>% dplyr::filter(id %in% sample_data$id)

  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
  
  
  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand-Receptor expression in samples\n\n") +
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
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(receiver, lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
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
  max_activity = abs(group_data$activity_scaled) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.40, 0.50, 0.60, 0.70, 0.825, 1),  limits = c(-1*max_activity, max_activity))

  p2 = p2 + custom_scale_fill

  p3 = group_data %>%
    ggplot(aes(receiver, lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
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
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))

  p3 = p3 + custom_scale_fill


  if(!is.null(widths)){
    p = patchwork::wrap_plots(
      p1,p2,p3,
      nrow = 1,guides = "collect",
      widths = widths
    )
  } else {
    p = patchwork::wrap_plots(
      p1,p2,p3,
      nrow = 1,guides = "collect",
      widths = c(sample_data$sample %>% unique() %>% length(), sample_data$receiver %>% unique() %>% length(), sample_data$receiver %>% unique() %>% length())
    )
  }


  return(p)


}


#' @title make_sample_lr_prod_activity_covariate_plots
#'
#' @description \code{make_sample_lr_prod_activity_covariate_plots}  XXXX EXPERIMENTAL FEATURE!!!!!!!!!!!!!!!!!!!!NOT WORKINGN RIGHT NOW
#' @usage make_sample_lr_prod_activity_covariate_plots(prioritization_tables, prioritized_tbl_oi, grouping_tbl, covariate_oi, widths = NULL, heights = NULL)
#'
#' @param prioritization_tables XXX
#' @param prioritized_tbl_oi XXX
#' @param grouping_tbl Data frame linking the sample_id, group_id and covariate_oi
#' @param covariate_oi Name of the covariate/batch that needs to be visualized for each sample
#' @param widths Vector of 3 elements: Width of the LR exprs product panel,  width of the scaled ligand activity panel, width of the ligand activity panel. Default NULL: automatically defined based on nr of samples vs nr of group-receiver combinations. If manual change: example format: c(5,1,1) 
#' @param heights Vector of 2 elements: Height of the covariate panel and height of the ligand-receptor prod+activity panel. Default NULL: automatically defined based on the nr of Ligand-Receptor pairs. If manual change: example format: c(1,5) 
#'
#' @return XXXX
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_sample_lr_prod_activity_covariate_plots = function(prioritization_tables, prioritized_tbl_oi, grouping_tbl, covariate_oi, widths = NULL, heights = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>%  dplyr::arrange(sender) %>% dplyr::group_by(sender) %>%  dplyr::arrange(receiver)
  
  group_data = prioritization_tables$group_prioritization_tbl %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity, activity_scaled, fraction_ligand_group, prioritization_score, scaled_avg_exprs_ligand) %>% dplyr::filter(id %in% sample_data$id)
  
  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = scaled_LR_frac)) +
    geom_point() +
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand-Receptor expression in samples\n\n") +
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
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Scaled L-R\navg exprs fraction product")
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  
  
  p2 = group_data %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(receiver, lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
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
  max_activity = abs(group_data$activity_scaled) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.40, 0.50, 0.60, 0.70, 0.825, 1),  limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  p3 = group_data %>%
    ggplot(aes(receiver, lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
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
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))
  
  p3 = p3 + custom_scale_fill
  
  grouping_tbl_plot = grouping_tbl %>% mutate(covariate_ = paste0(" ",covariate_oi," "), mock = "", covariate = grouping_tbl[[covariate_oi]])
  p_covariate = grouping_tbl_plot %>% 
    ggplot(aes(sample, mock, fill = covariate)) +
    geom_tile(color = "black") +
    facet_grid(covariate_~group, scales = "free", space = "free", switch = "y") +
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
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + scale_fill_brewer(palette = "Set2")
  
  if(!is.null(widths)){
    design <- "D##
               ABC"
    p = patchwork::wrap_plots(
      A = p1, B = p2, C= p3, D = p_covariate,
      guides = "collect", design = design,
      widths = widths,
      heights = heights
    ) + patchwork::plot_layout(design = design)
  } else {
    design <- "D##
               ABC"
    p = patchwork::wrap_plots(
      A = p1, B = p2, C= p3, D = p_covariate,
      guides = "collect", design = design,
      widths = c(sample_data$sample %>% unique() %>% length(), sample_data$receiver %>% unique() %>% length(), sample_data$receiver %>% unique() %>% length()),
      heights = c(1, group_data$lr_interaction %>% unique() %>% length(), )
    ) + patchwork::plot_layout(design = design)
  }
  
  
  # p = cowplot::plot_grid(p_covariate, NULL, NULL, p1, p2, p3, align = "vh", rel_widths = widths, rel_heights =  heights, axis = "tblr")
  
  
  return(p)
  
  
}

#' @title make_ligand_activity_plots
#'
#' @description \code{make_ligand_activity_plots}  Visualize the ligand activities (normal and scaled) of each group-receiver combination
#' @usage make_ligand_activity_plots(prioritization_tables, ligands_oi, contrast_tbl, widths = NULL)
#'
#' @param prioritization_tables XXX
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
#' print("XXXX")
#' }
#'
#' @export
#'
make_ligand_activity_plots = function(prioritization_tables, ligands_oi, contrast_tbl, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  group_data = prioritization_tables$ligand_activities_target_de_tbl %>% dplyr::inner_join(contrast_tbl) %>% dplyr::distinct(group, ligand, receiver, activity, activity_scaled) %>% dplyr::filter(ligand %in% ligands_oi)


  p1 = group_data %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(receiver, ligand, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
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
  max_activity = abs(group_data$activity_scaled) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.40, 0.50, 0.60, 0.70, 0.825, 1),  limits = c(-1*max_activity, max_activity))

  p1 = p1 + custom_scale_fill

  p2 = group_data %>%
    ggplot(aes(receiver, ligand, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(.~group, scales = "free", space = "free") +
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
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))

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

#' @title make_sample_target_plots
#'
#' @description \code{make_sample_target_plots}  Heatmap/Dotplot of scaled target gene expression per sample, compared between the groups of interest. (samples in columns, genes in rows)
#' @usage make_sample_target_plots(receiver_info, targets_oi, receiver_oi, grouping_tbl)
#'
#' @param receiver_info XXX
#' @param targets_oi Character vector of genes of which the expression should be visualized
#' @param receiver_oi Name of the receiver cell type of interest for which expression of genes should be visualized
#' @param grouping_tbl Data frame linking the sample_id, group_id and covariate_oi 
#'
#' @return Heatmap/Dotplot of scaled target gene expression (samples in columns, genes in rows)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nichenetr scaling_zscore
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_sample_target_plots = function(receiver_info, targets_oi, receiver_oi, grouping_tbl){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  avg_df =  receiver_info$avg_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi)
  frq_df =  receiver_info$frq_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi)

  filtered_data = avg_df %>% dplyr::inner_join(frq_df) %>% dplyr::inner_join(grouping_tbl)

  filtered_data = filtered_data %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled_target_exprs = nichenetr::scaling_zscore(average_sample), scaled_target_frac = nichenetr::scaling_zscore(fraction_sample)) %>% dplyr::ungroup()
  filtered_data$gene = factor(filtered_data$gene, levels=targets_oi)

  p1 = filtered_data %>%
    ggplot(aes(sample, gene, color = scaled_target_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(.~group, scales = "free", space = "free") +
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
    ) + labs(color = "Scaled target\navg expression", size= "Target\n exprs fraction")
  max_lfc = abs(filtered_data$scaled_target_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  p1 = p1 + custom_scale_fill
  return(p1)
}

#' @title make_sample_target_plots_reversed
#'
#' @description \code{make_sample_target_plots_reversed}  Heatmap/Dotplot of scaled target gene expression per sample, compared between the groups of interest. (genes in columns, samples in rows)
#' @usage make_sample_target_plots_reversed(receiver_info, targets_oi, receiver_oi, grouping_tbl)
#'
#' @inheritParams make_sample_target_plots
#'
#' @return Heatmap/Dotplot of scaled target gene expression (samples in rows, genes in columns)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nichenetr scaling_zscore
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_sample_target_plots_reversed = function(receiver_info, targets_oi, receiver_oi, grouping_tbl){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  avg_df =  receiver_info$avg_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi)
  frq_df =  receiver_info$frq_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% receiver_oi)

  filtered_data = avg_df %>% dplyr::inner_join(frq_df) %>% dplyr::inner_join(grouping_tbl)

  filtered_data = filtered_data %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled_target_exprs = nichenetr::scaling_zscore(average_sample), scaled_target_frac = nichenetr::scaling_zscore(fraction_sample)) %>% dplyr::ungroup()

  filtered_data$gene = factor(filtered_data$gene, levels=targets_oi)

  p1 = filtered_data %>%
    # ggplot(aes(gene, sample , fill = scaled_target_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(gene, sample , color = scaled_target_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(group~., scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
      axis.text.y = element_text(size = 9),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 11, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) +
    # labs(color = "Scaled target\navg expression")
    labs(color = "Scaled target\navg expression", size= "Target\n exprs fraction")
  max_lfc = abs(filtered_data$scaled_target_exprs) %>% max()
  # custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  plot = p1 + custom_scale_fill


  return( plot  )
}

#' @title make_group_lfc_exprs_activity_plot
#'
#' @description \code{make_group_lfc_exprs_activity_plot}  XXXX
#' @usage make_group_lfc_exprs_activity_plot(prioritization_tables, prioritized_tbl_oi, receiver_oi, heights = NULL)
#'
#' @param prioritization_tables XXX
#' @param prioritized_tbl_oi XXX
#' @param receiver_oi XXX
#' @param heights XXX
#'
#' @return XXXX
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork plot_layout wrap_plots
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_group_lfc_exprs_activity_plot = function(prioritization_tables, prioritized_tbl_oi, receiver_oi, heights = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  filtered_data = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id & receiver == receiver_oi) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>%  dplyr::arrange(sender) %>% dplyr::group_by(sender) %>%  dplyr::arrange(receiver)
  plot_data = prioritization_tables$group_prioritization_tbl %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " to ")) %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, ligand_receptor_lfc_avg, activity, activity_scaled, fraction_ligand_group, prioritization_score, scaled_avg_exprs_ligand) %>% dplyr::filter(lr_interaction %in% filtered_data$lr_interaction & receiver == receiver_oi)

  p_exprs = plot_data %>%
    ggplot(aes(lr_interaction, group, color = ligand_receptor_lfc_avg, size = scaled_avg_exprs_ligand)) +
    geom_point() +
    facet_grid(sender_receiver~., scales = "free", space = "free") +
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
    ) + labs(color = "Average of\nLogFC ligand & \nLogFC receptor", size= "Scaled average of\nligand expression")
  max_lfc = abs(plot_data$ligand_receptor_lfc_avg) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  p_exprs = p_exprs + custom_scale_fill

  p_activity_scaled = plot_data %>%
    ggplot(aes(lr_interaction, group, fill = activity_scaled)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "purple", mid = "whitesmoke", high = "orange") +
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
    ) + labs(fill = "Scaled NicheNet\nligand activity")
  max_activity = abs(plot_data$activity_scaled) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.40, 0.50, 0.60, 0.70, 0.825, 1),  limits = c(-1*max_activity, max_activity))
  p_activity_scaled = p_activity_scaled + custom_scale_fill

  p_activity = plot_data %>%
    ggplot(aes(lr_interaction, group, fill = activity)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "white", mid = "whitesmoke", high = "orange") +
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
    ) + labs(fill = "NicheNet\nligand activity")
  max_activity = (plot_data$activity) %>% max()
  min_activity = (plot_data$activity) %>% min()
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "OrRd") %>% .[-7]),values = c(0, 0.25, 0.40, 0.55, 0.70, 0.825, 1),  limits = c(min_activity-0.01, max_activity))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))

  p_activity = p_activity + custom_scale_fill

  if(!is.null(heights)){
    # heights should be a vector of lenght 2
    p = patchwork::wrap_plots(
      p_exprs ,
      p_activity + theme(axis.text.x = element_blank()),
      p_activity_scaled + theme(axis.text.x = element_blank()),
      nrow = 3,guides = "collect"
    ) + patchwork::plot_layout(heights = heights)
  } else {
    p = patchwork::wrap_plots(
      p_exprs ,
      p_activity + theme(axis.text.x = element_blank()),
      p_activity_scaled + theme(axis.text.x = element_blank()),
      nrow = 3,guides = "collect"
    ) + patchwork::plot_layout(heights = c(plot_data$sender %>% unique() %>% length(), 1, 1))
  }

  return(p)
}

#' @title make_circos_group_comparison
#'
#' @description \code{make_circos_group_comparison}  XXXX
#' @usage make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
#'
#' @param prioritized_tbl_oi XXX
#' @param colors_sender XXX
#' @param colors_receiver XXX
#'
#' @return XXXX
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
#' print("XXXX")
#' }
#'
#' @export
#'
make_circos_group_comparison = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  # Link each cell type to a color
  grid_col_ligand = colors_sender

  grid_col_receptor = colors_receiver

  grid_col_tbl_ligand = tibble::tibble(sender = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_receptor = tibble::tibble(receiver = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

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
    #circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score) %>% dplyr::mutate(ligand = paste(sender, ligand, sep = "_"), receptor = paste(receptor, receiver, sep = "_"))

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

    # print(intersecting_ligands_receptors)

    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(receptor, " ", sep = ""))
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
  # grid_col_all = c(grid_col_receptor, grid_col_ligand)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_receptor[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))

  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_ligand[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))

  p_legend = grDevices::recordPlot()

  all_plots$legend = p_legend

  return(all_plots)
}
#' @title make_featureplot
#'
#' @description \code{make_featureplot}  Make a panel of UMAP plots showing the expression of a gene of interest in the group of interest vs other groups
#' @usage make_featureplot(sce_subset_oi, sce_subset_bg, title_umap, gene_oi, group_oi, background_groups, group_id)
#'
#' @param sce_subset_oi SingleCellExperiment object containing the cells of the celltype and group of interest
#' @param sce_subset_bg SingleCellExperiment object containing the background cells (cells not of interest)
#' @param title_umap Character vector: title of the umap
#' @param gene_oi Character vector: name of the gene of interest
#' @param group_oi Character vector: name of the group of interest
#' @param background_groups Character vector: names of the background group(s)
#' @param group_id Name of the colData(sce) column in which the id of the group is defined
#'
#' @return UMAP plots showing the expression of gene of interest, compared between a group of interest and background groups
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom scater plotReducedDim
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_featureplot = function(sce_subset_oi, sce_subset_bg, title_umap, gene_oi, group_oi, background_groups, group_id){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  p_dim = scater::plotReducedDim(sce_subset_oi,  "UMAP", colour_by = group_id) + ggtitle(title_umap) + theme(title = element_text(face = "bold"))

  p_oi = scater::plotReducedDim(sce_subset_oi,  "UMAP", colour_by = gene_oi)
  p_bg = scater::plotReducedDim(sce_subset_bg,  "UMAP", colour_by = gene_oi)

  p_oi = p_oi + ggtitle(paste(gene_oi, group_oi, sep = " in ")) + theme(title = element_text(face = "bold"))
  p_bg = p_bg + ggtitle(paste(gene_oi, background_groups %>% paste0(collapse = " & "), sep = " in ")) + theme(title = element_text(face = "bold"))
  
  wrapped_plots = patchwork::wrap_plots(p_dim,
                                        p_oi,
                                        p_bg,
                                        ncol = 3,guides = "collect")
}

#' @title make_ligand_receptor_feature_plot
#'
#' @description \code{make_ligand_receptor_feature_plot}  UMAP showing ligand and receptor expression in groups and celltypes of interest vs other groups and cell types.
#' @usage make_ligand_receptor_feature_plot(sce_sender, sce_receiver, ligand_oi, receptor_oi, group_oi, group_id, celltype_id_sender, celltype_id_receiver, senders_oi, receivers_oi, background_groups = NULL)
#'
#' @param ligand_oi Character vector of name of the ligand of interest
#' @param receptor_oi Character vector of name of the receptor of interest
#' @param group_oi Character vector of name of the group of interest
#' @param senders_oi Character vector with the names of the sender cell types of interest
#' @param receivers_oi Character vector with the names of the receiver cell types of interest
#' @param background_groups Default NULL: all groups in the group_id metadata column will be chosen. But user can fill in a character vector with the names of all gruops of interest.
#' @inheritParams multi_nichenet_analysis_separate
#'
#' @return Plot object: UMAP showing ligand and receptor expression in groups and celltypes of interest vs other groups and cell types.
#'
#' @import ggplot2
#' @import dplyr

#' @importFrom RColorBrewer brewer.pal
#' @importFrom generics setdiff
#' @importFrom patchwork wrap_plots
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_ligand_receptor_feature_plot = function(sce_sender, sce_receiver, ligand_oi, receptor_oi, group_oi, group_id, celltype_id_sender, celltype_id_receiver, senders_oi, receivers_oi, background_groups = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce_sender)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  # subset Sender - Nebulosa
  sce_subset_oi = sce_sender[, SummarizedExperiment::colData(sce_sender)[,group_id] %in% group_oi]
  sce_subset_bg = sce_sender[, SummarizedExperiment::colData(sce_sender)[,group_id] %in% background_groups]

  sce_subset_oi =  sce_subset_oi[, SummarizedExperiment::colData(sce_subset_oi)[,celltype_id_sender] %in% senders_oi]
  sce_subset_bg =  sce_subset_bg[, SummarizedExperiment::colData(sce_subset_bg)[,celltype_id_sender] %in% senders_oi]

  sender_plots_feature = make_featureplot(sce_subset_oi,sce_subset_bg, "Sender UMAP", ligand_oi, group_oi, background_groups, group_id)
  
  # subset Receiver - Nebulosa
  sce_subset_oi = sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,group_id] %in% group_oi]
  sce_subset_bg = sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,group_id] %in% background_groups]
  
  sce_subset_oi =  sce_subset_oi[, SummarizedExperiment::colData(sce_subset_oi)[,celltype_id_receiver] %in% receivers_oi]
  sce_subset_bg =  sce_subset_bg[, SummarizedExperiment::colData(sce_subset_bg)[,celltype_id_receiver] %in% receivers_oi]

  receiver_plots_feature = make_featureplot(sce_subset_oi, sce_subset_bg, "Receiver UMAP", receptor_oi, group_oi, background_groups, group_id)
  
  p_feature = patchwork::wrap_plots(sender_plots_feature, receiver_plots_feature, nrow = 2)

  return(p_feature)
}

#' @title make_ligand_receptor_violin_plot
#'
#' @description \code{make_ligand_receptor_violin_plot}  Plot combining a violin plot of of the ligand of interest in the sender cell type of interest, and a violin plot of the receptor of interest in the receiver cell type of interest.
#
#' @usage make_ligand_receptor_violin_plot(sce_sender, sce_receiver, ligand_oi, receptor_oi, sender_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id_sender, celltype_id_receiver, background_groups = NULL)
#'
#' @param sample_id Name of the colData(sce) column in which the id of the sample is defined
#' @param sender_oi Character vector with the names of the sender cell type of interest
#' @param receiver_oi  Character vector with the names of the receiver cell type of interest
#' @inheritParams make_ligand_receptor_feature_plot
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
#' print("XXXX")
#' }
#'
#' @export
#'
make_ligand_receptor_violin_plot = function(sce_sender, sce_receiver, ligand_oi, receptor_oi, sender_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id_sender, celltype_id_receiver, background_groups = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce_sender)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  # ligand plot
  sce_subset =  sce_sender[, SummarizedExperiment::colData(sce_sender)[,celltype_id_sender] %in% sender_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id_receiver, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[ligand_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[ligand_oi,])) %>% dplyr::inner_join(SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble())
  
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
  
  # receptor plot
  sce_subset =  sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id_receiver, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[receptor_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[receptor_oi,])) %>% dplyr::inner_join(SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble())
  
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
  
  return(patchwork::wrap_plots(p_sender, p_receiver, nrow = 2))

}

#' @title make_target_violin_plot
#'
#' @description \code{make_target_violin_plot}  Violin plot of a target gene of interest: per sample, and samples are grouped per group
#' @usage make_target_violin_plot(sce_receiver, target_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id_receiver, prioritized_tbl_oi, background_groups = NULL)
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
#' print("XXXX")
#' }
#'
#' @export
#'
make_target_violin_plot = function(sce_receiver, target_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id_receiver, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce_receiver)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  sce_subset =  sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,celltype_id_receiver] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]

  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id_receiver, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[target_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[target_oi,])) %>% dplyr::inner_join(SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble())
  
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
  

  return(p_violin)

}

#' @title make_target_feature_plot
#'
#' @description \code{make_target_feature_plot}  UMAP showing target gene expression in groups and celltypes of interest vs other groups and cell types.
#' @usage make_target_feature_plot(sce_receiver, target_oi, group_oi, group_id, celltype_id_receiver, receivers_oi, prioritized_tbl_oi)
#'
#' @param target_oi Character vector of the name of the target gene of interest
#' @inheritParams make_ligand_receptor_feature_plot
#'
#' @return UMAP showing target gene expression in groups and celltypes of interest vs other groups and cell types.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom SummarizedExperiment colData

#' @importFrom generics setdiff
#'
#' @examples
#' \dontrun{
#' print("XXXX")
#' }
#'
#' @export
#'
make_target_feature_plot = function(sce_receiver, target_oi, group_oi, group_id, celltype_id_receiver, receivers_oi, prioritized_tbl_oi, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce_receiver)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }

  # subset Receiver - Nebulosa
  sce_subset_oi = sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,group_id] %in% group_oi]
  sce_subset_bg = sce_receiver[, SummarizedExperiment::colData(sce_receiver)[,group_id] %in% background_groups]
  
  sce_subset_oi =  sce_subset_oi[, SummarizedExperiment::colData(sce_subset_oi)[,celltype_id_receiver] %in% receivers_oi]
  sce_subset_bg =  sce_subset_bg[, SummarizedExperiment::colData(sce_subset_bg)[,celltype_id_receiver] %in% receivers_oi]
  
  receiver_plots_feature = make_featureplot(sce_subset_oi, sce_subset_bg, "Receiver UMAP", target_oi, group_oi, background_groups, group_id)
  
  p_feature = receiver_plots_feature
  
  return(p_feature)
}

#' @title make_ligand_activity_target_plot
#'
#' @description \code{make_ligand_activity_target_plot}  XXXX
#' @usage make_ligand_activity_target_plot(group_oi, receiver_oi, prioritized_tbl_oi, ligand_activities_targets_DEgenes, contrast_tbl, grouping_tbl, receiver_info, plot_legend = TRUE, heights = NULL, widths = NULL)
#'
#' @param group_oi XXX
#' @param receiver_oi XXX
#' @param prioritized_tbl_oi XXX
#' @param ligand_activities_targets_DEgenes XXX
#' @param contrast_tbl XXX
#' @param grouping_tbl XXX
#' @param receiver_info XXX
#' @param plot_legend XXX
#' @param heights XXX
#' @param widths XXX
#'
#' @return XXXX
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
#' print("XXXX")
#' }
#'
#' @export
#'
make_ligand_activity_target_plot = function(group_oi, receiver_oi, prioritized_tbl_oi, ligand_activities_targets_DEgenes, contrast_tbl, grouping_tbl, receiver_info, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()

  # Ligand-Target heatmap
  active_ligand_target_links_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )

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

  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

  # Ligand-Activity-Scaled
  ligand_pearson_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl)  %>% dplyr::select(ligand, group, activity_scaled) %>% dplyr::distinct() %>% tidyr::spread(group, activity_scaled)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), ] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "PuRd"),values = c(0, 0.35, 0.425, 0.525, 0.625, 0.75, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill

  # Ligand-Activity
  ligand_pearson_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl)  %>% dplyr::select(ligand, group, activity) %>% dplyr::distinct() %>% tidyr::spread(group, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), ] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "Oranges"),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill

  # Target expression
  groups_oi = contrast_tbl %>% dplyr::pull(group) %>% unique()
  p_targets = make_sample_target_plots_reversed(receiver_info, order_targets_, receiver_oi, grouping_tbl %>% dplyr::filter(group %in% groups_oi))

  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_samples = grouping_tbl %>% dplyr::filter(group %in% groups_oi) %>% dplyr::pull(sample) %>% length()

  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))

  if(is.null(heights)){
    heights = c(n_ligands + 3, n_samples)
  }
  if(is.null(widths)){
    widths = c(n_groups + 1.33, n_groups, n_targets)
  }

  if(plot_legend == FALSE){
    design <- "AaB
               ##C"
    combined_plot = patchwork::wrap_plots(A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))

  } else {
    design <- "AaB
               L#C"

    combined_plot = patchwork::wrap_plots(A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }

}
