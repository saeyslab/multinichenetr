---
title: "MultiNicheNet analysis: include additional data modalities: serum proteomics in MIS-C"
author: "Robin Browaeys"
package: "multinichenetr 2.0.0"
output: 
  BiocStyle::html_document
output_dir: "/Users/robinb/Work/multinichenetr/vignettes"
vignette: >
  %\VignetteIndexEntry{MultiNicheNet analysis: include additional data modalities: serum proteomics in MIS-C}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8} 
date: 13 May 2024
link-citations: true
---

<style type="text/css">
.smaller {
  font-size: 10px
}
</style>

<!-- github markdown built using
rmarkdown::render("vignettes/add_proteomics_MISC.Rmd", clean = FALSE )
-->



In this vignette, you can learn how to extend the default MultiNicheNet prioritization scheme to **include suitable additional data modalities for prioritizing cell-cell communication patterns**. If you are new to MultiNicheNet, we recommend learning the basics first by reading and running this vignette: [basis_analysis_steps_MISC.knit.md](basis_analysis_steps_MISC.knit.md). In general we would always recommend users to run the analysis step-by-step, but the wrapper function is applied here to reduce the length of the vignette. In this vignette, we will use the MultiNicheNet wrapper function to perform all core steps of the analysis in one line of code. Afterwards, we will demonstrate one possible way to score interactions based on additional data modalities, and how you can incorporate these scorings in the prioritization. 

As example expression data of interacting cells, we will here use scRNAseq data of immune cells in MIS-C patients and healthy siblings from this paper of Hoste et al.: [TIM3+ TRBV11-2 T cells and IFNγ signature in patrolling monocytes and CD16+ NK cells delineate MIS-C](https://rupress.org/jem/article/219/2/e20211381/212918/TIM3-TRBV11-2-T-cells-and-IFN-signature-in) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6362434.svg)](https://doi.org/10.5281/zenodo.6362434). 
MIS-C (multisystem inflammatory syndrome in children) is a novel rare immunodysregulation syndrome that can arise after SARS-CoV-2 infection in children. We will use MultiNicheNet to explore immune cell crosstalk enriched in MIS-C compared to healthy siblings and adult COVID-19 patients. 

In addition to this dataset, we will use publicly available processed OLINK **serum proteomics data** from Diorio et al (https://www.nature.com/articles/s41467-021-27544-6) who profiled MIS-C patients and healthy controls. This data could be informative for prioritizing cell-cell interactions because it may give us information about which ligands are more strongly expressed at the protein level in MIS-C patients versus healthy siblings.   

In this vignette, we will first prepare the MultiNicheNet core analysis, run the MultiNicheNet core analysis, then extend the prioritization with serum proteomics data, and finally interpret that output. All steps before the the incorporation of proteomics data, are exactly the same as described in this vignette: [pairwise_analysis_MISC.knit.md](pairwise_analysis_MISC.knit.md)

# Preparation of the MultiNicheNet core analysis


```r
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
```

## Load NicheNet's ligand-receptor network and ligand-target matrix

The Nichenet v2 networks and matrices for both mouse and human can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7074291.svg)](https://doi.org/10.5281/zenodo.7074291). 

We will read these object in for human because our expression data is of human patients. 
Gene names are here made syntactically valid via `make.names()` to avoid the loss of genes (eg H2-M3) in downstream visualizations.


```r
organism = "human"
```


```r
options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
      )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
    ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
    ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}
```

## Read in SingleCellExperiment Object

In this vignette, we will load in a subset of the scRNAseq data of the MIS-C [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8010790.svg)](https://doi.org/10.5281/zenodo.8010790). For the sake of demonstration, this subset only contains 3 cell types. These celltypes are some of the cell types that were found to be most interesting related to MIS-C according to Hoste et al. 

Because the NicheNet 2.0. networks are in the most recent version of the official gene symbols, we will make sure that the gene symbols used in the expression data are also updated (= converted from their "aliases" to official gene symbols). Afterwards, we will make them again syntactically valid. 


```r
sce = readRDS(url(
  "https://zenodo.org/record/8010790/files/sce_subset_misc.rds"
  ))
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
```

## Prepare the settings of the MultiNicheNet cell-cell communication analysis

In this step, we will formalize our research question into MultiNicheNet input arguments.

### Define in which metadata columns we can find the **group**, **sample** and **cell type** IDs

In this case study, we want to study differences in cell-cell communication patterns between MIS-C patients (M) and their healthy siblings (S). The meta data columns that indicate this disease status(=group/condition of interest) is `MIS.C.AgeTier`. 

Cell type annotations are indicated in the `Annotation_v2.0` column, and the sample is indicated by the `ShortID` column. 
If your cells are annotated in multiple hierarchical levels, we recommend using a relatively high level in the hierarchy. This for 2 reasons: 1) MultiNicheNet focuses on differential expression and not differential abundance, and 2) there should be sufficient cells per sample-celltype combination (see later).


```r
sample_id = "ShortID"
group_id = "MIS.C.AgeTier"
celltype_id = "Annotation_v2.0"
```

__Important__: It is required that each sample-id is uniquely assigned to only one condition/group of interest. See the vignettes about paired and multifactorial analysis to see how to define your analysis input when you have multiple samples (and conditions) per patient.

If you would have batch effects or covariates you can correct for, you can define this here as well. However, this is not applicable to this dataset. Therefore we will use the following NA settings:


```r
covariates = NA
batches = NA
```

__Important__: for categorical covariates and batches, there should be at least one sample for every group-batch combination. If one of your groups/conditions lacks a certain level of your batch, you won't be able to correct for the batch effect because the model is then not able to distinguish batch from group/condition effects.

__Important__: The column names of group, sample, cell type, batches and covariates should be syntactically valid (`make.names`)

__Important__: All group, sample, cell type, batch and covariate names should be syntactically valid as well (`make.names`) (eg through `SummarizedExperiment::colData(sce)$ShortID = SummarizedExperiment::colData(sce)$ShortID %>% make.names()`)

### Define the contrasts of interest.

Here, we want to compare patients in the MIS-C (M) group vs healthy control siblings (S) (= M vs S). We want to know which cell-cell communication patterns are specific for the M vs S group, and the S vs M group. 

To perform this comparison, we need to set the following contrasts:


```r
contrasts_oi = c("'M-S','S-M'") 
```

__Very Important__ Note the format to indicate the contrasts! This formatting should be adhered to very strictly, and white spaces are not allowed! Check `?get_DE_info` for explanation about how to define this well. The most important points are that: 
*each contrast is surrounded by single quotation marks
*contrasts are separated by a comma without any white space 
*all contrasts together are surrounded by double quotation marks. 

If you compare against two groups, you should divide by 2 (as demonstrated here), if you compare against three groups, you should divide by 3 and so on.

For downstream visualizations and linking contrasts to their main condition, we also need to run the following:
This is necessary because we will also calculate cell-type+condition specificity of ligands and receptors. 


```r
contrast_tbl = tibble(
  contrast = c("M-S","S-M"), 
  group = c("M","S")
  )
```

### Define the sender and receiver cell types of interest.

If you want to focus the analysis on specific cell types (e.g. because you know which cell types reside in the same microenvironments based on spatial data), you can define this here. If you have sufficient computational resources and no specific idea of cell-type colocalzations, we recommend to consider all cell types as potential senders and receivers. 

Here we will consider all cell types in the data:


```r
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]
```

In case you would have samples in your data that do not belong to one of the groups/conditions of interest, we recommend removing them and only keeping conditions of interst:


```r
conditions_keep = c("M", "S")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
          ]
```

# Running the MultiNicheNet core analysis

Now we will run the core of a MultiNicheNet analysis through the `multi_nichenet_analysis` wrapper function. Under the hood, the following steps are executed:

* 1. Cell-type filtering: determine which cell types are sufficiently present
* 2. Gene filtering: determine which genes are sufficiently expressed in each present cell type
* 3. Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type
* 4. Differential expression (DE) analysis: determine which genes are differentially expressed
* 5. Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes
* 6. Prioritization: rank cell-cell communication patterns through multi-criteria prioritization
* 7. Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
* 8. Prioritize communication patterns involving condition-specific cell types through an alternative prioritization scheme

After these steps, the output can be further explored as we will demonstrate in the "Downstream analysis of the MultiNicheNet output" section.

However, before we can start running the analysis, we need to define some parameters first.

## Defining parameters for the MultiNicheNet wrapper function

#### Parameters for step 1: Cell-type filtering:

Since MultiNicheNet will infer group differences at the sample level for each cell type (currently via Muscat - pseudobulking + EdgeR), we need to have sufficient cells per sample of a cell type, and this for all groups. In the following analysis we will set this minimum number of cells per cell type per sample at 10. Samples that have less than `min_cells` cells will be excluded from the analysis for that specific cell type.


```r
min_cells = 10
```

####  Parameters for step 2: Gene filtering

For each cell type, we will consider genes expressed if they are expressed in at least a `min_sample_prop` fraction of samples in the condition with the lowest number of samples. By default, we set `min_sample_prop = 0.50`. 


```r
min_sample_prop = 0.50
```

But how do we define which genes are expressed in a sample? For this we will consider genes as expressed if they have non-zero expression values in a `fraction_cutoff` fraction of cells of that cell type in that sample. By default, we set `fraction_cutoff = 0.05`, which means that genes should show non-zero expression values in at least 5% of cells in a sample. 


```r
fraction_cutoff = 0.05
```

####  Parameters for step 4: DE analysis

Output of the DE analysis step will be used to define the levels of differential expression of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver cell types (which will be used for the ligand activity analysis). We need to define whether we will use the regular or empirical p-values returned by this DE analysis. 
By default, we will use the regular p-values. After running the analysis, users can inspect the DE output p-value distributions to assess possible deviations to the underlying assumption. Only when you observe aberrant distributions we recommend redoing the analysis with the empirical p-values instead of the regular p-values. In our experience, this is only required for a limited number of datasets. 


```r
empirical_pval = FALSE
```

####  Parameters for step 5: Ligand activity prediction

One of the prioritization criteria is the predicted activity of ligands in receiver cell types. Similarly to base NicheNet (https://github.com/saeyslab/nichenetr), we use the DE output to create a "geneset of interest": here we assume that DE genes within a cell type may be DE because of differential cell-cell communication processes. To determine the genesets of interest based on DE output, we need to define which logFC and/or p-value thresholds we will use. 

By default, we will apply the p-value cutoff on the normal p-values, and not on the p-values corrected for multiple testing. This choice was made because most multi-sample single-cell transcriptomics datasets have just a few samples per group and we might have a lack of statistical power due to pseudobulking. But, if the smallest group >= 20 samples, we typically recommend using p_val_adj = TRUE. When the biological difference between the conditions is very large, we typically recommend increasing the logFC_threshold and/or using p_val_adj = TRUE.


```r
logFC_threshold = 0.50
p_val_threshold = 0.05
```


```r
p_val_adj = FALSE 
```

After the ligand activity prediction, we will also infer the predicted target genes of these ligands in each contrast. For this ligand-target inference procedure, we also need to select which top n of the predicted target genes will be considered (here: top 250 targets per ligand). This parameter will not affect the ligand activity predictions. It will only affect ligand-target visualizations and construction of the intercellular regulatory network during the downstream analysis. We recommend users to test other settings in case they would be interested in exploring fewer, but more confident target genes, or vice versa. 


```r
top_n_target = 250
```

The NicheNet ligand activity analysis can be run in parallel for each receiver cell type, by changing the number of cores as defined here. Using more cores will speed up the analysis at the cost of needing more memory. This is only recommended if you have many receiver cell types of interest. You can define here the maximum number of cores you want to be used. 


```r
n.cores = 8
```

####  Parameters for step 6: Prioritization

We will use the following criteria to prioritize ligand-receptor interactions:

* Upregulation of the ligand in a sender cell type and/or upregulation of the receptor in a receiver cell type - in the condition of interest.
* Cell-type specific expression of the ligand in the sender cell type and receptor in the receiver cell type in the condition of interest (to mitigate the influence of upregulated but still relatively weakly expressed ligands/receptors). 
* Sufficiently high expression levels of ligand and receptor in many samples of the same group.
* High NicheNet ligand activity, to further prioritize ligand-receptor pairs based on their predicted effect of the ligand-receptor interaction on the gene expression in the receiver cell type. 

We will combine these prioritization criteria in a single aggregated prioritization score. In the default setting, we will weigh each of these criteria equally (`scenario = "regular"`). This setting is strongly recommended. However, we also provide some additional setting to accomodate different biological scenarios. The setting `scenario = "lower_DE"` halves the weight for DE criteria and doubles the weight for ligand activity. This is recommended in case your hypothesis is that the differential CCC patterns in your data are less likely to be driven by DE (eg in cases of differential migration into a niche). The setting `scenario = "no_frac_LR_expr"` ignores the criterion "Sufficiently high expression levels of ligand and receptor in many samples of the same group". This may be interesting for users that have data with a limited number of samples and don’t want to penalize interactions if they are not sufficiently expressed in some samples. 

Here we will choose for the regular setting.


```r
scenario = "regular"
```

Finally, we still need to make one choice. For NicheNet ligand activity we can choose to prioritize ligands that only induce upregulation of target genes (`ligand_activity_down = FALSE`) or can lead potentially lead to both up- and downregulation (`ligand_activity_down = TRUE`). The benefit of `ligand_activity_down = FALSE` is ease of interpretability: prioritized ligand-receptor pairs will be upregulated in the condition of interest, just like their target genes.  `ligand_activity_down = TRUE` can be harder to interpret because target genes of some interactions may be upregulated in the other conditions compared to the condition of interest. This is harder to interpret, but may help to pick up interactions that can also repress gene expression. 

Here we will choose for setting `ligand_activity_down = FALSE` and focus specifically on upregulating ligands.


```r
ligand_activity_down = FALSE
```

## Running the MultiNicheNet wrapper function


```r
multinichenet_output = multi_nichenet_analysis(
  sce = sce, 
  celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
  batches = batches, covariates = covariates, 
  lr_network = lr_network, ligand_target_matrix = ligand_target_matrix, 
  contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl, 
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
  )
## [1] "Make diagnostic abundance plots + define expressed genes"
## [1] "condition-specific celltypes:"
## character(0)
## [1] "absent celltypes:"
## character(0)
## [1] "Samples are considered if they have more than 10 cells of the cell type of interest"
## [1] "Genes with non-zero counts in at least 5% of cells of a cell type of interest in a particular sample will be considered as expressed in that sample."
## [1] "Genes expressed in at least 2.5 samples will considered as expressed in the cell type: L_NK_CD56._CD16."
## [1] "Genes expressed in at least 2.5 samples will considered as expressed in the cell type: L_T_TIM3._CD38._HLADR."
## [1] "Genes expressed in at least 2.5 samples will considered as expressed in the cell type: M_Monocyte_CD16"
## [1] "6010 genes are considered as expressed in the cell type: L_NK_CD56._CD16."
## [1] "7589 genes are considered as expressed in the cell type: L_T_TIM3._CD38._HLADR."
## [1] "7798 genes are considered as expressed in the cell type: M_Monocyte_CD16"
## [1] "Calculate differential expression for all cell types"
## [1] "DE analysis is done:"
## [1] "included cell types are:"
## [1] "L_T_TIM3._CD38._HLADR." "L_NK_CD56._CD16."       "M_Monocyte_CD16"       
## [1] "retained cell types"
## [1] "L_T_TIM3._CD38._HLADR." "L_NK_CD56._CD16."       "M_Monocyte_CD16"       
## [1] "Calculate normalized average and pseudobulk expression"
## [1] "Calculate NicheNet ligand activities and ligand-target links"
## [1] "Combine all the information in prioritization tables"
## [1] "Calculate correlation between LR pairs and target genes"
## [1] "There are no condition specific cell types in the data. MultiNicheNet analysis is performed in the regular way for all cell types."
```

# Inspecting the MultiNicheNet output

## Inspecting output tables

Before incorporating the serum proteomics, we will quickly inspect the output of this regular MultiNicheNet analysis. Please check other vignettes for a more detailed breakdown of these different steps.

## Tables with the final prioritization scores (results per group and per sample)


```r
multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()
## # A tibble: 6 × 18
##   contrast group sender                 receiver               ligand  receptor lr_interaction id      scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##   <chr>    <chr> <chr>                  <chr>                  <chr>   <chr>    <chr>          <chr>               <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
## 1 M-S      M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. S100A8  CD69     S100A8_CD69    S100A8…             0.952                  0.910               0.963                  0.998               0.888             1.00               1.00                  0.857                0.943 M        
## 2 M-S      M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16        IFNG    IFNGR1   IFNG_IFNGR1    IFNG_I…             0.977                  0.885               0.848                  0.854               1.00              1.00               1.00                  0.857                0.940 M        
## 3 M-S      M     M_Monocyte_CD16        M_Monocyte_CD16        TIMP1   CD63     TIMP1_CD63     TIMP1_…             0.973                  0.993               0.985                  0.993               0.749             1.00               1.00                  0.857                0.930 M        
## 4 M-S      M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16        IFNG    IFNGR2   IFNG_IFNGR2    IFNG_I…             0.977                  0.885               0.726                  0.846               1.00              1.00               1.00                  0.857                0.930 M        
## 5 M-S      M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA LAG3     HLA.DRA_LAG3   HLA.DR…             0.791                  0.861               0.978                  0.981               0.766             1.00               1.00                  1                    0.929 M        
## 6 M-S      M     M_Monocyte_CD16        M_Monocyte_CD16        S100A9  CD68     S100A9_CD68    S100A9…             0.924                  0.906               0.790                  0.953               0.732             1.00               1.00                  1                    0.920 M        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
```

## Visualization of differential cell-cell interactions

### Summarizing ChordDiagram circos plots

In a first instance, we will look at the broad overview of prioritized interactions via condition-specific Chordiagram circos plots.

We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.


```r
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
  )
```


```r
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-88-1.png" width="100%" /><img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-88-2.png" width="100%" /><img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-88-3.png" width="100%" />

### Interpretable bubble plots

Whereas these ChordDiagrams show the most specific interactions per group, they don't give insights into the data behind these predictions. Therefore we will now look at visualizations that indicate the different prioritization criteria used in MultiNicheNet. 

In the next type of plots, we will 1) visualize the per-sample scaled product of normalized ligand and receptor pseudobulk expression, 2) visualize the scaled ligand activities, 3) cell-type specificity. 

We will now check the top 50 interactions specific for the MIS-C group


```r
group_oi = "M"
```


```r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  groups_oi = group_oi)
```


```r
plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50)
plot_oi
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-91-1.png" width="100%" />

# Using serum proteomics data to prioritize cell-cell communication patterns

Here we will demonstrate how to tailor the prioritization to additional data modalities that you may have. How to do this exactly will strongly depend on the type of data modality. Here we show an example on how you can incorporate "targeted" serum proteomics data.

As additional data modality, we will use serum proteomics data (generated through the OLINK platform) that was published by a different research group based on serum samples of a different cohort of MIS-C patients (Diorio et al. - https://www.nature.com/articles/s41467-021-27544-6). 
Ligands that are also upregulated in the serum at the protein level may be interesting candidates for further follow-up. Therefore, we can use this data to help in the prioritization. We will here add an additional score to each interaction related to how strongly the ligand is upregulated at the protein level in the serum of patients. 

Because we don't have OLINK information of each ligand, we will consider ligands without data as being not-differentially expressed at the protein level (p-val = 1, logFC = 0). By doing this, these ligands will still be kept in the total prioritization analysis such that they could still be ranked high if the previous criteria strongly point to their importance. 

Whereas we will only consider ligand upregulation for prioritization, we will still keep receptor information in the OLINK dataframe for later visualizations.

## Read in additional data modality - here: OLINK serum proteomics data 

Of this data, we will only keep information about ligands and receptors


```r
file_path = "summary_difference_diorioMIS-C-vs-HC.xlsx"
download.file("https://zenodo.org/records/10908003/files/summary_difference_diorioMIS-C-vs-HC.xlsx", 
              file_path, 
              mode = "wb")
olink_df = xlsx::read.xlsx(file_path, 1) %>% 
  as_tibble() %>% 
  dplyr::rename(gene = variables, logFC = mean_diff_d0, pval = adj_p_values) %>% 
  dplyr::select(gene, logFC, pval) %>% 
  dplyr::mutate(gene = gene %>% make.names()) %>% 
  mutate(olink_type = "Diorio")
olink_df %>% head()
## # A tibble: 6 × 4
##   gene         logFC     pval olink_type
##   <chr>        <dbl>    <dbl> <chr>     
## 1 PLA2G2A       7.01 2.16e-21 Diorio    
## 2 IL1RL1        4.59 3.77e-18 Diorio    
## 3 MPO           3.26 3.77e-18 Diorio    
## 4 DEFA1_DEFA1B  3.26 4.32e-18 Diorio    
## 5 FRZB         -2.24 4.32e-18 Diorio    
## 6 LILRB4        2.60 4.55e-18 Diorio
```

Let's inspect the nr of of ligands and receptors in the OLINK data frame:

```r
olink_df %>% filter(gene %in% lr_network$ligand) %>% nrow() # nr of ligands in OLINK data
## [1] 425
olink_df %>% filter(gene %in% lr_network$receptor) %>% nrow() # nr of receptors in OLINK data
## [1] 310
```

Now create a data frame with DE values for all ligands and receptors. If they were not assessed by the OLINK platform, we will consider them as being not-differentially expressed (p-val = 1, logFC = 0)


```r
olink_df_ligands_receptors = olink_df %>% 
  filter(gene %in% union(lr_network$ligand, lr_network$receptor))

olink_df_ligands_receptors_NAs = tibble(
  gene = union(
    lr_network$ligand, lr_network$receptor
    ) %>% setdiff(olink_df_ligands_receptors$gene), 
  logFC = 0, 
  pval = 1) 

olink_df_ligands_receptors = bind_rows(
  olink_df_ligands_receptors, 
  olink_df_ligands_receptors_NAs)
olink_df_ligands_receptors %>% head()
## # A tibble: 6 × 4
##   gene    logFC     pval olink_type
##   <chr>   <dbl>    <dbl> <chr>     
## 1 PLA2G2A  7.01 2.16e-21 Diorio    
## 2 IL1RL1   4.59 3.77e-18 Diorio    
## 3 LILRB4   2.60 4.55e-18 Diorio    
## 4 CXCL10   5.64 5.35e-18 Diorio    
## 5 CALCA    5.47 2.67e-16 Diorio    
## 6 HAVCR2   2.10 2.67e-16 Diorio
```

## Score sender/ligand-receiver/receptor interactions based on the alternative data modality

Important: the provided logFC values in the OLINK data represent the contrast MIS-C versus healthy siblings (M-vs-S). We will still need to add the reverse logFC for the contrast S-vs-M. 


```r
olink_df_ligand = multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  distinct(ligand) %>% 
  left_join(olink_df_ligands_receptors %>% rename(ligand = gene)) %>% 
  mutate(contrast = contrast_tbl$contrast[1])

olink_df_reverse_ligand = olink_df_ligand %>% 
  mutate(contrast = contrast_tbl$contrast[2], logFC = -1*logFC)
olink_df_ligand = olink_df_ligand %>% 
  bind_rows(olink_df_reverse_ligand)
olink_df_ligand %>% filter(ligand == "IFNG") # sanity check
## # A tibble: 2 × 5
##   ligand logFC      pval olink_type contrast
##   <chr>  <dbl>     <dbl> <chr>      <chr>   
## 1 IFNG    3.36 0.0000776 Diorio     M-S     
## 2 IFNG   -3.36 0.0000776 Diorio     S-M
olink_df_ligand %>% arrange(-logFC) # sanity check
## # A tibble: 388 × 5
##    ligand   logFC     pval olink_type contrast
##    <chr>    <dbl>    <dbl> <chr>      <chr>   
##  1 GZMB      4.09 5.64e- 7 Diorio     M-S     
##  2 CEACAM21  3.43 8.95e- 8 Diorio     M-S     
##  3 IFNG      3.36 7.76e- 5 Diorio     M-S     
##  4 CD14      3.34 2.47e- 6 Diorio     M-S     
##  5 SPON2     3.17 3.05e- 7 Diorio     M-S     
##  6 BST2      3.00 6.84e- 9 Diorio     M-S     
##  7 NID1      2.97 1.34e-12 Diorio     M-S     
##  8 CCL3      2.95 4.52e-11 Diorio     M-S     
##  9 TNFSF14   2.77 4.86e-14 Diorio     M-S     
## 10 S100A12   2.66 5.23e-11 Diorio     M-S     
## # ℹ 378 more rows
```

We will now score ligand-contrast combinations based on their OLINK logFC and p-values. We will do this in the exact same way as we did for differential expression at the RNA level.  


```r
olink_df_ligand = olink_df_ligand %>% 
  mutate(pval_adapted = -log10(pval)*sign(logFC)) 

olink_df_ligand_scaled = olink_df_ligand %>% 
  distinct(contrast, ligand, logFC, pval_adapted) %>% 
  dplyr::mutate(
    scaled_lfc_ligand_OLINK = 
      rank(logFC, ties.method = "average", na.last = FALSE)/
      max(rank(logFC, ties.method = "average", na.last = FALSE)), 
    scaled_p_val_ligand_adapted_OLINK = 
      rank(pval_adapted, ties.method = "average", na.last = FALSE)/
      max(rank(pval_adapted, ties.method = "average", na.last = FALSE))) 

olink_df_ligand_scaled = olink_df_ligand_scaled %>% 
  distinct(contrast, ligand, scaled_lfc_ligand_OLINK, scaled_p_val_ligand_adapted_OLINK) %>%
  dplyr::arrange(-scaled_lfc_ligand_OLINK)

olink_df_ligand_scaled %>% head()
## # A tibble: 6 × 4
##   contrast ligand   scaled_lfc_ligand_OLINK scaled_p_val_ligand_adapted_OLINK
##   <chr>    <chr>                      <dbl>                             <dbl>
## 1 M-S      GZMB                       1                                 0.892
## 2 M-S      CEACAM21                   0.997                             0.902
## 3 M-S      IFNG                       0.995                             0.856
## 4 M-S      CD14                       0.992                             0.879
## 5 M-S      SPON2                      0.990                             0.894
## 6 M-S      BST2                       0.987                             0.925
```

As you can see here, we scored each contrast-ligand combination by a between-0-and-1 scaled value for the logFC and p-value of differential expression in the OLINK proteomics data. The higher this value, the stronger the DE at protein level, and thus the stronger this interaction will be prioritized for this criterion.

## Calculate the new prioritization score by including the score of the criterion based on the additional data modality

Now we will calculate the aggregated prioritization score with inclusion of the scaled logFC and scaled p-value-based score from the OLINK data. 
To do this, we first add these new scores to the existing prioritization_table. 


```r
multinichenet_output$prioritization_tables$group_prioritization_tbl = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>% 
  inner_join(olink_df_ligand_scaled)
```

Then we will define the table of prioritization criteria and add the new criteria and their weights. It is recommended to give the new criteria the same weights as the original default. If you ran the analysis with the default weights (`scenario = "regular"` then these weights are 1)  


```r
regular_criteria_tbl = tibble(
  criterion = 
    c("scaled_lfc_ligand",
      "scaled_p_val_ligand_adapted",
      "scaled_lfc_receptor",
      "scaled_p_val_receptor_adapted", 
      "max_scaled_activity", 
      "scaled_pb_ligand", 
      "scaled_pb_receptor", 
      "fraction_expressing_ligand_receptor"
      ), 
  weight = NA, 
  regularization_factor = c(0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1)
  ) # do not change this  - 
# logFC and pval of ligand and receptor get 0.5 as regularization such that DE of the ligand counts as one score, 
# which is the sum of the scaled_score of the logFC and the p-value-based score. 
# This also applies to the receptor.

new_criteria_tbl = tibble(
  criterion = c("scaled_lfc_ligand_OLINK","scaled_p_val_ligand_adapted_OLINK"), 
  weight = c(1,1), 
  regularization_factor = c(0.5, 0.5)) # 0.5 as regularization factor
# such that DE of the ligand at the protein level counts as one score, 
# which is the sum of the scaled_score of the logFC and the p-value-based score. 
```

Now recalculate the aggregated prioritization score


```r
group_prioritization_tbl_withoutOLINK = multinichenet_output$prioritization_tables$group_prioritization_tbl 

prioritization_tables = add_extra_criterion(
  multinichenet_output$prioritization_tables, 
  new_criteria_tbl, 
  regular_criteria_tbl, 
  scenario = "regular"
  ) 
```


```r
prioritization_tables$group_prioritization_tbl %>% head(20)
## # A tibble: 20 × 20
##    contrast group sender          receiver ligand receptor lr_interaction id    scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group scaled_lfc_ligand_OL…⁴
##    <chr>    <chr> <chr>           <chr>    <chr>  <chr>    <chr>          <chr>             <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>                      <dbl>
##  1 M-S      M     L_T_TIM3._CD38… M_Monoc… IFNG   IFNGR1   IFNG_IFNGR1    IFNG…             0.977                  0.885               0.848                  0.854               1.00             1.00               1.00                   0.857                0.938 M                          0.995
##  2 M-S      M     M_Monocyte_CD16 M_Monoc… TIMP1  CD63     TIMP1_CD63     TIMP…             0.973                  0.993               0.985                  0.993               0.749            1.00               1.00                   0.857                0.936 M                          0.941
##  3 M-S      M     L_T_TIM3._CD38… M_Monoc… IFNG   IFNGR2   IFNG_IFNGR2    IFNG…             0.977                  0.885               0.726                  0.846               1.00             1.00               1.00                   0.857                0.929 M                          0.995
##  4 M-S      M     L_NK_CD56._CD1… L_T_TIM… CCL4   CCR5     CCL4_CCR5      CCL4…             0.876                  0.930               0.900                  0.769               0.923            1.00               1.00                   0.857                0.924 M                          0.956
##  5 M-S      M     L_NK_CD56._CD1… L_T_TIM… CCL3   CCR5     CCL3_CCR5      CCL3…             0.979                  0.989               0.900                  0.769               0.780            1.00               1.00                   0.857                0.918 M                          0.982
##  6 M-S      M     M_Monocyte_CD16 M_Monoc… LILRB2 IFNGR1   LILRB2_IFNGR1  LILR…             0.948                  0.987               0.848                  0.854               0.751            1.00               1.00                   0.857                0.911 M                          0.923
##  7 M-S      M     M_Monocyte_CD16 L_T_TIM… CXCL16 CXCR6    CXCL16_CXCR6   CXCL…             0.749                  0.789               0.973                  0.918               0.867            1.00               1.00                   0.857                0.909 M                          0.905
##  8 M-S      M     M_Monocyte_CD16 L_T_TIM… HLA.D… LAG3     HLA.DRA_LAG3   HLA.…             0.791                  0.861               0.978                  0.981               0.766            1.00               1.00                   1                    0.908 M                          0.789
##  9 S-M      S     L_T_TIM3._CD38… M_Monoc… CD28   CD4      CD28_CD4       CD28…             0.853                  0.852               0.815                  0.912               0.789            1.00               1.00                   1                    0.898 S                          0.781
## 10 M-S      M     L_NK_CD56._CD1… M_Monoc… BST2   LILRA5   BST2_LILRA5    BST2…             0.798                  0.894               0.959                  0.988               0.585            1.00               1.00                   0.857                0.888 M                          0.987
## 11 M-S      M     M_Monocyte_CD16 M_Monoc… C1QA   CR1      C1QA_CR1       C1QA…             0.998                  0.996               0.845                  0.808               0.821            1.00               1.00                   0.714                0.887 M                          0.851
## 12 M-S      M     L_T_TIM3._CD38… L_T_TIM… CCL4   CCR5     CCL4_CCR5      CCL4…             0.990                  0.984               0.900                  0.769               0.923            0.634              1.00                   0.857                0.883 M                          0.956
## 13 M-S      M     M_Monocyte_CD16 L_NK_CD… HLA.E  KLRC1    HLA.E_KLRC1    HLA.…             0.602                  0.714               0.963                  0.991               1.00             0.685              1.00                   1                    0.881 M                          0.838
## 14 M-S      M     M_Monocyte_CD16 L_T_TIM… IL15   IL2RG    IL15_IL2RG     IL15…             0.604                  0.586               0.847                  0.870               0.919            1.00               0.985                  0.857                0.881 M                          0.928
## 15 M-S      M     M_Monocyte_CD16 L_T_TIM… S100A8 CD69     S100A8_CD69    S100…             0.952                  0.910               0.963                  0.998               0.888            1.00               1.00                   0.857                0.880 M                          0.501
## 16 M-S      M     L_NK_CD56._CD1… M_Monoc… CCL3   CCR1     CCL3_CCR1      CCL3…             0.979                  0.989               0.961                  0.960               0.665            1.00               1.00                   0.571                0.878 M                          0.982
## 17 M-S      M     M_Monocyte_CD16 M_Monoc… LILRB2 SIGLEC9  LILRB2_SIGLEC9 LILR…             0.948                  0.987               0.773                  0.735               0.751            1.00               1.00                   0.714                0.877 M                          0.923
## 18 M-S      M     M_Monocyte_CD16 L_T_TIM… ICAM1  IL2RG    ICAM1_IL2RG    ICAM…             0.818                  0.792               0.847                  0.870               0.836            1.00               0.985                  0.714                0.876 M                          0.907
## 19 S-M      S     M_Monocyte_CD16 M_Monoc… ENG    BMPR2    ENG_BMPR2      ENG_…             0.925                  0.934               0.877                  0.844               0.743            1.00               1.00                   0.8                  0.875 S                          0.778
## 20 M-S      M     M_Monocyte_CD16 L_T_TIM… LGALS3 LAG3     LGALS3_LAG3    LGAL…             0.579                  0.606               0.978                  0.981               0.879            1.00               1.00                   0.857                0.874 M                          0.814
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor, ⁴​scaled_lfc_ligand_OLINK
## # ℹ 1 more variable: scaled_p_val_ligand_adapted_OLINK <dbl>
```

Here we see the updated prioritization table with updated scores after inclusion of the additional data modality.

We can also make a table where we compare the original prioritization scores with the updated ones.


```r
comparison_table = group_prioritization_tbl_withoutOLINK %>% select(group, id, prioritization_score) %>%
  inner_join(
    prioritization_tables$group_prioritization_tbl %>% 
      select(group, id, prioritization_score) %>% 
      rename(new_score = prioritization_score)) %>% 
  mutate(difference = new_score - prioritization_score) 
```

Let's now inspect the top interactions that got the highest increase in ranking because of the OLINK data


```r
comparison_table %>% 
  arrange(-difference) %>% 
  filter(new_score > 0.80)
## # A tibble: 168 × 5
##    group id                                                      prioritization_score new_score difference
##    <chr> <chr>                                                                  <dbl>     <dbl>      <dbl>
##  1 M     LILRB4_LAIR1_M_Monocyte_CD16_L_T_TIM3._CD38._HLADR.                    0.777     0.807     0.0296
##  2 M     TNFSF13_TNFRSF14_M_Monocyte_CD16_M_Monocyte_CD16                       0.771     0.800     0.0296
##  3 M     RETN_TLR4_M_Monocyte_CD16_M_Monocyte_CD16                              0.787     0.815     0.0282
##  4 M     CCL3_CCR5_L_T_TIM3._CD38._HLADR._L_T_TIM3._CD38._HLADR.                0.776     0.804     0.0271
##  5 M     CCL3_CCR1_L_NK_CD56._CD16._L_T_TIM3._CD38._HLADR.                      0.777     0.804     0.0271
##  6 M     TNFSF13_FAS_M_Monocyte_CD16_M_Monocyte_CD16                            0.795     0.821     0.0261
##  7 M     TIMD4_SIGLEC9_L_T_TIM3._CD38._HLADR._M_Monocyte_CD16                   0.789     0.815     0.0257
##  8 M     TIMD4_SIGLEC7_L_T_TIM3._CD38._HLADR._M_Monocyte_CD16                   0.796     0.821     0.0247
##  9 M     TIMP1_CD63_M_Monocyte_CD16_L_T_TIM3._CD38._HLADR.                      0.798     0.822     0.0245
## 10 M     BST2_LILRA5_M_Monocyte_CD16_M_Monocyte_CD16                            0.790     0.813     0.0238
## # ℹ 158 more rows
```

Let's now inspect the top interactions that got penalized the most because of the OLINK data. Note here: interactions can be penalized because they are less upregulated at the protein level, or because there was no protein expression information for them.


```r
comparison_table %>% 
  arrange(difference) %>% 
  filter(new_score > 0.80)
## # A tibble: 168 × 5
##    group id                                                 prioritization_score new_score difference
##    <chr> <chr>                                                             <dbl>     <dbl>      <dbl>
##  1 S     TGFBI_ITGA4_M_Monocyte_CD16_M_Monocyte_CD16                       0.919     0.809    -0.110 
##  2 S     TGFBI_ITGB1_M_Monocyte_CD16_L_T_TIM3._CD38._HLADR.                0.917     0.807    -0.109 
##  3 M     S100A8_CD69_M_Monocyte_CD16_L_T_TIM3._CD38._HLADR.                0.943     0.880    -0.0631
##  4 M     S100A9_CD68_M_Monocyte_CD16_M_Monocyte_CD16                       0.920     0.860    -0.0598
##  5 M     S100A8_CD68_M_Monocyte_CD16_M_Monocyte_CD16                       0.912     0.853    -0.0586
##  6 S     HLA.DQA1_CD4_M_Monocyte_CD16_M_Monocyte_CD16                      0.909     0.851    -0.0582
##  7 S     SEMA4A_PLXNB2_M_Monocyte_CD16_M_Monocyte_CD16                     0.905     0.848    -0.0577
##  8 M     HLA.F_KLRC1_M_Monocyte_CD16_L_NK_CD56._CD16.                      0.901     0.844    -0.0572
##  9 S     GNAS_ADRB2_M_Monocyte_CD16_L_NK_CD56._CD16.                       0.901     0.844    -0.0570
## 10 M     HLA.F_LILRB1_M_Monocyte_CD16_M_Monocyte_CD16                      0.900     0.843    -0.0570
## # ℹ 158 more rows
```

Let's now look at some visualizations of the updated prioritized condition-specific interactions.

## Visualization of the updated cell-cell communication prioritization

In a first instance, we will look at the broad overview of prioritized interactions via condition-specific Circos plots.

### Circos plot of top-prioritized links

We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.


```r
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  rank_per_group = FALSE)
```


```r
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(
  prioritized_tbl_oi$sender %>% unique(), 
  prioritized_tbl_oi$receiver %>% unique()) %>% 
  sort()

colors_sender = RColorBrewer::brewer.pal(
  n = length(senders_receivers), name = 'Spectral') %>% 
  magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(
  n = length(senders_receivers), name = 'Spectral') %>% 
  magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(
  prioritized_tbl_oi, 
  colors_sender, colors_receiver)
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-105-1.png" width="100%" /><img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-105-2.png" width="100%" /><img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-105-3.png" width="100%" />

### Visualization of scaled ligand-receptor pseudobulk products and ligand activity

Now we will visualize per sample the scaled product of ligand and receptor expression. Samples that were left out of the DE analysis are indicated with a smaller dot (this helps to indicate the samples that did not contribute to the calculation of the logFC, and thus not contributed to the final prioritization)

We will first visualize the interactions specific for the Sibling group that are part of the top50 differential interactions as shown in the above circos plots


```r
group_oi = "S"
```


```r
prioritized_tbl_oi_S_50 = prioritized_tbl_oi_all %>% filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_S_50)
plot_oi
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-107-1.png" width="100%" />

We will now check the MIS-C specific interactions


```r
group_oi = "M"
```


```r
prioritized_tbl_oi = prioritized_tbl_oi_all %>% 
  filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi)
plot_oi
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-109-1.png" width="100%" />

Whereas these are the classical plots with default multinichenetr code, we can also adapt this function to visualize the additional criteria that we used for prioritization. 

### Add additional data modality information to the interpretable bubble blot visualization

We will show example code to make this multi-panel plot step-by-step

First: ligand-receptor pseudobulk product expression panel


```r
# ligand-receptor pseudobulk product expression panel
sample_data = prioritization_tables$sample_prioritization_tbl %>% 
  dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
  dplyr::mutate(
    sender_receiver = paste(sender, receiver, sep = " --> "), 
    lr_interaction = paste(ligand, receptor, sep = " - ")) %>%
  dplyr::arrange(receiver) %>% 
  dplyr::group_by(receiver) %>%  
  dplyr::arrange(sender, .by_group = TRUE)

sample_data = sample_data %>% 
  dplyr::mutate(sender_receiver = factor(
    sender_receiver, 
    levels = sample_data$sender_receiver %>% unique()
    ))
  
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
  custom_scale_fill = scale_color_gradientn(
    colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
    values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  
    limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
p1
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-110-1.png" width="100%" />
Second: scaled ligand activity based on RNA data


```r
# scaled ligand activity based on RNA
group_data = prioritization_tables$group_prioritization_table_source  %>% 
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "), 
      lr_interaction = paste(ligand, receptor, sep = " - "))  %>% 
    dplyr::distinct(id, sender, receiver, sender_receiver, ligand, receptor, lr_interaction, group, activity_scaled, direction_regulation, prioritization_score) %>% 
  dplyr::filter(id %in% sample_data$id) %>% 
  dplyr::arrange(receiver) %>% 
  dplyr::group_by(receiver) %>% 
  dplyr::arrange(sender, .by_group = TRUE)

  group_data = group_data %>% dplyr::mutate(
    sender_receiver = factor(
      sender_receiver, 
      levels = group_data$sender_receiver %>% unique()
      ))
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
    ) + labs(fill = "Scaled Ligand\nActivity RNA")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(
    colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),
    values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  
    limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
p2
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-111-1.png" width="100%" />

Third: OLINK logFC MISC-vs-Sibling


```r
olink_df = olink_df %>% 
  mutate(contrast = contrast_tbl$contrast[1])

olink_df_reverse = olink_df %>% 
  mutate(contrast = contrast_tbl$contrast[2], logFC = -1*logFC)
olink_df = olink_df %>% 
  bind_rows(olink_df_reverse)
olink_df %>% filter(gene == "IFNG") # sanity check
## # A tibble: 2 × 5
##   gene  logFC      pval olink_type contrast
##   <chr> <dbl>     <dbl> <chr>      <chr>   
## 1 IFNG   3.36 0.0000776 Diorio     M-S     
## 2 IFNG  -3.36 0.0000776 Diorio     S-M
olink_df %>% arrange(-logFC) # sanity check
## # A tibble: 2,926 × 5
##    gene     logFC     pval olink_type contrast
##    <chr>    <dbl>    <dbl> <chr>      <chr>   
##  1 PLA2G2A   7.01 2.16e-21 Diorio     M-S     
##  2 NTproBNP  6.02 7.31e-14 Diorio     M-S     
##  3 CXCL11    5.79 2.03e-13 Diorio     M-S     
##  4 CCL7      5.72 2.16e-15 Diorio     M-S     
##  5 CXCL10    5.64 5.35e-18 Diorio     M-S     
##  6 MMP8      5.53 1.77e-16 Diorio     M-S     
##  7 CALCA     5.47 2.67e-16 Diorio     M-S     
##  8 VSIG4     5.34 7.42e-16 Diorio     M-S     
##  9 OSM       4.85 3.24e-12 Diorio     M-S     
## 10 CXCL9     4.64 8.18e-15 Diorio     M-S     
## # ℹ 2,916 more rows
```

```r
group_data_olink = 
  inner_join(
    
    group_data %>% 
      left_join(
        olink_df %>% filter(contrast == "M-S") %>% 
          select(gene, olink_type, logFC) %>% 
          dplyr::rename(ligand = gene, `ligand ` = logFC)
        ) %>% 
      left_join(
        olink_df %>% filter(contrast == "M-S") %>% 
          select(gene, olink_type, logFC) %>% 
          dplyr::rename(receptor = gene, `receptor ` = logFC)
        ) %>% 
      tidyr::gather(LR, logFC,`ligand `:`receptor `),
    
    group_data %>% 
      left_join(
        olink_df %>% filter(contrast == "M-S") %>% 
          select(gene, olink_type, pval) %>% 
          dplyr::rename(ligand = gene, `ligand ` = pval)
        ) %>% 
      left_join(
        olink_df %>% filter(contrast == "M-S") %>% 
          select(gene, olink_type, pval) %>% 
          dplyr::rename(receptor = gene, `receptor ` = pval)
        ) %>% 
      tidyr::gather(LR, pval,`ligand `:`receptor `) %>% mutate(neg_log_pval = -log10(pval))
    
  )
  
group_data_olink = group_data_olink %>% mutate(olink_type = "Diorio")
  
p_olink = group_data_olink %>% 
    ggplot(aes(LR , lr_interaction, color = logFC, size = neg_log_pval)) +
    geom_point() +
    facet_grid(sender_receiver ~ olink_type, scales = "free", space = "free") +
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
    ) + labs(color = "OLink LogFC M-vs-S") + labs(size = "OLink -log10(pval)")
  max_lfc = abs(group_data_olink$logFC) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_color_gradientn(
    colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
    values = c(0, 0.30, 0.40, 0.5, 0.60, 0.70, 1), 
    limits = c(-1*max_lfc, max_lfc))
  
p_olink = p_olink + custom_scale_fill + scale_size_binned_area(max_size = 4)
p_olink
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-113-1.png" width="100%" />

And now we will put everything together:

```r
p = patchwork::wrap_plots(
      p1,p2,p_olink,
      nrow = 1,guides = "collect",
      widths = c(sample_data$sample %>% unique() %>% length(), 2*(sample_data$group %>% unique() %>% length()), 2*length(unique(olink_df$olink_type)))
)
p
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-114-1.png" width="100%" />

Let's now put all this code to generate this visualization in a function we will use to visualize the interactions of which the ranking was most strongly affected by the incorporation of the additional data modality.


```r
make_sample_lr_prod_activity_OLINK_plots = function(prioritization_tables, prioritized_tbl_oi){
  # ligand-receptor pseudobulk product expression panel
  sample_data = prioritization_tables$sample_prioritization_tbl %>% 
    dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "), 
      lr_interaction = paste(ligand, receptor, sep = " - ")) %>%
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>%  
    dplyr::arrange(sender, .by_group = TRUE)
  
  sample_data = sample_data %>% 
    dplyr::mutate(sender_receiver = factor(
      sender_receiver, 
      levels = sample_data$sender_receiver %>% unique()
      ))
    
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
    custom_scale_fill = scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
      values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  
      limits = c(-1*max_lfc, max_lfc))
    
    p1 = p1 + custom_scale_fill
  p1
  # scaled ligand activity based on RNA
  group_data = prioritization_tables$group_prioritization_table_source  %>% 
      dplyr::mutate(
        sender_receiver = paste(sender, receiver, sep = " --> "), 
        lr_interaction = paste(ligand, receptor, sep = " - "))  %>% 
      dplyr::distinct(id, sender, receiver, sender_receiver, ligand, receptor, lr_interaction, group, activity_scaled, direction_regulation, prioritization_score) %>% 
    dplyr::filter(id %in% sample_data$id) %>% 
    dplyr::arrange(receiver) %>% 
    dplyr::group_by(receiver) %>% 
    dplyr::arrange(sender, .by_group = TRUE)
  
    group_data = group_data %>% dplyr::mutate(
      sender_receiver = factor(
        sender_receiver, 
        levels = group_data$sender_receiver %>% unique()
        ))
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
      ) + labs(fill = "Scaled Ligand\nActivity RNA")
    max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
    custom_scale_fill = scale_fill_gradientn(
      colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),
      values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  
      limits = c(-1*max_activity, max_activity))
    
    p2 = p2 + custom_scale_fill
  p2
  group_data_olink = 
    inner_join(
      
      group_data %>% 
        left_join(
          olink_df %>% filter(contrast == "M-S") %>% 
            select(gene, olink_type, logFC) %>% 
            dplyr::rename(ligand = gene, `ligand ` = logFC)
          ) %>% 
        left_join(
          olink_df %>% filter(contrast == "M-S") %>% 
            select(gene, olink_type, logFC) %>% 
            dplyr::rename(receptor = gene, `receptor ` = logFC)
          ) %>% 
        tidyr::gather(LR, logFC,`ligand `:`receptor `),
      
      group_data %>% 
        left_join(
          olink_df %>% filter(contrast == "M-S") %>% 
            select(gene, olink_type, pval) %>% 
            dplyr::rename(ligand = gene, `ligand ` = pval)
          ) %>% 
        left_join(
          olink_df %>% filter(contrast == "M-S") %>% 
            select(gene, olink_type, pval) %>% 
            dplyr::rename(receptor = gene, `receptor ` = pval)
          ) %>% 
        tidyr::gather(LR, pval,`ligand `:`receptor `) %>% mutate(neg_log_pval = -log10(pval))
      
    )
    
  group_data_olink = group_data_olink %>% mutate(olink_type = "Diorio")
    
  p_olink = group_data_olink %>% 
      ggplot(aes(LR , lr_interaction, color = logFC, size = neg_log_pval)) +
      geom_point() +
      facet_grid(sender_receiver ~ olink_type, scales = "free", space = "free") +
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
      ) + labs(color = "OLink LogFC M-vs-S") + labs(size = "OLink -log10(pval)")
    max_lfc = abs(group_data_olink$logFC) %>% max(na.rm = TRUE)
    custom_scale_fill = scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
      values = c(0, 0.30, 0.40, 0.5, 0.60, 0.70, 1), 
      limits = c(-1*max_lfc, max_lfc))
    
  p_olink = p_olink + custom_scale_fill + scale_size_binned_area(max_size = 4)
  p_olink
  p = patchwork::wrap_plots(
        p1,p2,p_olink,
        nrow = 1,guides = "collect",
        widths = c(sample_data$sample %>% unique() %>% length(), 2*(sample_data$group %>% unique() %>% length()), 2*length(unique(olink_df$olink_type)))
  )
  return(p)
}
```

Check whether we can reproduce the previous figure:

```r
make_sample_lr_prod_activity_OLINK_plots(
  prioritization_tables = prioritization_tables,
  prioritized_tbl_oi = prioritized_tbl_oi
  )
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-116-1.png" width="100%" />

For these top MIS-C specific interactions we can see: some interactions are supported by upregulation at the protein level of the ligand. However, other without complementary proteomics data did still end up in the top results overall because they scored very highly on the regular criteria.

Now visualize the interactions that most strongly benefited from the addition of OLINK data (top25 different interactions that are in the top1000 overall and have a new prioritization scorre > 0.80):


```r
ids_oi = comparison_table %>% 
    filter(new_score > 0.80) %>% top_n(25, difference) %>% pull(id)

prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  1000, 
  rank_per_group = TRUE, 
  groups_oi = "M"
  ) 

prioritized_tbl_oi = prioritized_tbl_oi_all %>% 
  filter(id %in% ids_oi) %>% filter(group == "M")
```


```r
plot_oi = make_sample_lr_prod_activity_OLINK_plots(
  prioritization_tables = prioritization_tables,
  prioritized_tbl_oi = prioritized_tbl_oi
  )
plot_oi
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-118-1.png" width="100%" />

As expected, these interactions are characterized by a strongly upregulated ligand at the protein ligand, whereas upregulation of the ligand-receptor pair at the RNA level and/or activity is not always very clear. 

Now visualize the interactions that were most strongly penalized by the addition of OLINK data - those were mainly interactions in the S-group:


```r
ids_oi = comparison_table %>% 
    filter(prioritization_score > 0.80) %>% top_n(25, -difference) %>% pull(id)

prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  1000, 
  rank_per_group = TRUE, groups_oi = "S") 
prioritized_tbl_oi = prioritized_tbl_oi_all %>% 
  filter(id %in% ids_oi) %>% filter(group == "S")
```


```r
plot_oi = make_sample_lr_prod_activity_OLINK_plots(
  prioritization_tables = prioritization_tables,
  prioritized_tbl_oi = prioritized_tbl_oi
  )
plot_oi
```

<img src="add_proteomics_MISC_files/figure-html/unnamed-chunk-120-1.png" width="100%" />

As expected for interactions penalized by the OLINK addition, we see an inconsistency between the RNA and protein level with respect to the direction of expression difference between M and S patients. Here, S-specific interactions at the RNA level have higher protein levels in the M-group instead of the S group. Because of this incongruency, these interactions are penalized in the final prioritization score.

# Tips for other data modalites

If you want to incorporate a different additional data modality, you can always open a github issue to discuss how to do this properly.

