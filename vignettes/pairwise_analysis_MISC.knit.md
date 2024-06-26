---
title: "MultiNicheNet analysis: MIS-C pairwise comparison - wrapper function"
author: "Robin Browaeys"
package: "multinichenetr 2.0.0"
output: 
  BiocStyle::html_document
output_dir: "/Users/robinb/Work/multinichenetr/vignettes"
vignette: >
  %\VignetteIndexEntry{MultiNicheNet analysis: MIS-C pairwise comparison - wrapper function}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8} 
date: 3 April 2024
link-citations: true
---

<style type="text/css">
.smaller {
  font-size: 10px
}
</style>

<!-- github markdown built using
rmarkdown::render("vignettes/pairwise_analysis_MISC.Rmd", clean = FALSE )
-->



In this vignette, you can learn how to perform a MultiNicheNet analysis to compare cell-cell communication between two conditions of interest (one-vs-one comparison). A MultiNicheNet analysis can be performed if you have multi-sample, multi-condition/group single-cell data. We strongly recommend having at least 4 samples in each of the groups/conditions you want to compare. With less samples, the benefits of performing a pseudobulk-based DE analysis are less clear. For those datasets, you can check and run our alternative workflow that makes use of cell-level sample-agnostic differential expression tools.

As input you need a SingleCellExperiment object containing at least the raw count matrix and metadata providing the following information for each cell: the **group**, **sample** and **cell type**.

As example expression data of interacting cells, we will here use scRNAseq data of immune cells in MIS-C patients and healthy siblings from this paper of Hoste et al.: [TIM3+ TRBV11-2 T cells and IFNγ signature in patrolling monocytes and CD16+ NK cells delineate MIS-C](https://rupress.org/jem/article/219/2/e20211381/212918/TIM3-TRBV11-2-T-cells-and-IFN-signature-in) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6362434.svg)](https://doi.org/10.5281/zenodo.6362434). 
MIS-C (multisystem inflammatory syndrome in children) is a novel rare immunodysregulation syndrome that can arise after SARS-CoV-2 infection in children. 

We will use MultiNicheNet to explore immune cell crosstalk enriched in MIS-C compared to healthy siblings. In this vignette we will demonstrate how to set the input of the analysis to perform this type of one-vs-one comparison. We will do this by using the MultiNicheNet wrapper function to perform all core steps of the analysis in one line of code. In general we would always recommend users to run the analysis step-by-step, but the wrapper function is applied here to reduce the length of the vignette. If you are new to MultiNicheNet and/or want to explore the different steps of MultiNicheNet one by one, we recommend reading and running this vignette: [basis_analysis_steps_MISC.knit.md](basis_analysis_steps_MISC.knit.md). 

In this vignette, we will first prepare the MultiNicheNet core analysis, then run the MultiNicheNet core analysis, and finally interpret the output.

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

# Downstream analysis of the MultiNicheNet output

## Inspecting the MultiNicheNet output

Before visualizing the major differential cell-cell communication patterns, we will have a look at the output of this regular MultiNicheNet analysis.

### Normalized pseudobulk expression for each cell type - sample combination


```r
multinichenet_output$celltype_info$pb_df %>% head()
## # A tibble: 6 × 4
##   gene    sample pb_sample celltype              
##   <chr>   <chr>      <dbl> <fct>                 
## 1 AAGAB   M1          5.35 L_T_TIM3._CD38._HLADR.
## 2 AAK1    M1          6.68 L_T_TIM3._CD38._HLADR.
## 3 ABCC1   M1          3.03 L_T_TIM3._CD38._HLADR.
## 4 ABHD13  M1          4.48 L_T_TIM3._CD38._HLADR.
## 5 ABHD17A M1          6.28 L_T_TIM3._CD38._HLADR.
## 6 ABI1    M1          6.14 L_T_TIM3._CD38._HLADR.
multinichenet_output$celltype_info$pb_df_group %>% head()
## # A tibble: 6 × 4
## # Groups:   group, celltype [1]
##   group celltype         gene   pb_group
##   <chr> <chr>            <chr>     <dbl>
## 1 M     L_NK_CD56._CD16. AAGAB     4.72 
## 2 M     L_NK_CD56._CD16. AAK1      6.60 
## 3 M     L_NK_CD56._CD16. ABCA1     0.896
## 4 M     L_NK_CD56._CD16. ABCC1     2.48 
## 5 M     L_NK_CD56._CD16. ABHD10    3.00 
## 6 M     L_NK_CD56._CD16. ABHD13    4.22
```

### DE information for each cell type - contrast combination


```r
multinichenet_output$celltype_de %>% head()
## # A tibble: 6 × 9
##   gene    cluster_id               logFC logCPM      F    p_val p_adj.loc  p_adj contrast
##   <chr>   <chr>                    <dbl>  <dbl>  <dbl>    <dbl>     <dbl>  <dbl> <chr>   
## 1 AAGAB   L_T_TIM3._CD38._HLADR.  0.679    5.02  5.1   0.036       0.524  0.524  M-S     
## 2 AAK1    L_T_TIM3._CD38._HLADR. -0.827    6.83 20     0.000269    0.0785 0.0785 M-S     
## 3 ABCC1   L_T_TIM3._CD38._HLADR.  0.271    3.77  0.255 0.62        0.976  0.976  M-S     
## 4 ABHD13  L_T_TIM3._CD38._HLADR. -0.173    4.31  0.217 0.646       0.98   0.98   M-S     
## 5 ABHD17A L_T_TIM3._CD38._HLADR.  0.131    6.2   0.226 0.64        0.978  0.978  M-S     
## 6 ABI1    L_T_TIM3._CD38._HLADR. -0.0616   6.15  0.107 0.748       1      1      M-S
```

### Output of the NicheNet ligand activity analysis, and the NicheNet ligand-target inference


```r
multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% head()
## # A tibble: 6 × 8
## # Groups:   receiver, contrast [1]
##   ligand activity contrast target ligand_target_weight receiver         direction_regulation
##   <chr>     <dbl> <chr>    <chr>                 <dbl> <chr>            <fct>               
## 1 A2M      0.0361 M-S      ACOT7               0.00632 L_NK_CD56._CD16. up                  
## 2 A2M      0.0361 M-S      AREG                0.00638 L_NK_CD56._CD16. up                  
## 3 A2M      0.0361 M-S      BCL2                0.0149  L_NK_CD56._CD16. up                  
## 4 A2M      0.0361 M-S      CD38                0.00621 L_NK_CD56._CD16. up                  
## 5 A2M      0.0361 M-S      CD55                0.00649 L_NK_CD56._CD16. up                  
## 6 A2M      0.0361 M-S      CD74                0.00634 L_NK_CD56._CD16. up                  
## # ℹ 1 more variable: activity_scaled <dbl>
```

### Tables with the final prioritization scores (results per group and per sample)


```r
multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()
## # A tibble: 6 × 18
##   contrast group sender       receiver ligand receptor lr_interaction id    scaled_lfc_ligand
##   <chr>    <chr> <chr>        <chr>    <chr>  <chr>    <chr>          <chr>             <dbl>
## 1 M-S      M     M_Monocyte_… L_T_TIM… S100A8 CD69     S100A8_CD69    S100…             0.952
## 2 M-S      M     L_T_TIM3._C… M_Monoc… IFNG   IFNGR1   IFNG_IFNGR1    IFNG…             0.977
## 3 M-S      M     M_Monocyte_… M_Monoc… TIMP1  CD63     TIMP1_CD63     TIMP…             0.973
## 4 M-S      M     L_T_TIM3._C… M_Monoc… IFNG   IFNGR2   IFNG_IFNGR2    IFNG…             0.977
## 5 M-S      M     M_Monocyte_… L_T_TIM… HLA.D… LAG3     HLA.DRA_LAG3   HLA.…             0.791
## 6 M-S      M     M_Monocyte_… M_Monoc… S100A9 CD68     S100A9_CD68    S100…             0.924
## # ℹ 9 more variables: scaled_p_val_ligand_adapted <dbl>, scaled_lfc_receptor <dbl>,
## #   scaled_p_val_receptor_adapted <dbl>, max_scaled_activity <dbl>, scaled_pb_ligand <dbl>,
## #   scaled_pb_receptor <dbl>, fraction_expressing_ligand_receptor <dbl>,
## #   prioritization_score <dbl>, top_group <chr>
```

Let's now generate visualize the major differential cell-cell communication patterns. 

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

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-26-1.png" width="100%" /><img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-26-2.png" width="100%" /><img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-26-3.png" width="100%" />

### Interpretable bubble plots

Whereas these ChordDiagrams show the most specific interactions per group, they don't give insights into the data behind these predictions. Therefore we will now look at visualizations that indicate the different prioritization criteria used in MultiNicheNet. 

In the next type of plots, we will 1) visualize the per-sample scaled product of normalized ligand and receptor pseudobulk expression, 2) visualize the scaled ligand activities, and 3) cell-type specificity. 

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

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-29-1.png" width="100%" />

Samples that were left out of the DE analysis are indicated with a smaller dot (this helps to indicate the samples that did not contribute to the calculation of the logFC, and thus not contributed to the final prioritization)

As a further help for further prioritization, we can assess the level of curation of these LR pairs as defined by the Intercellular Communication part of the Omnipath database


```r
prioritized_tbl_oi_M_50_omnipath = prioritized_tbl_oi_M_50 %>% 
  inner_join(lr_network_all)
```

Now we add this to the bubble plot visualization:

```r
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50_omnipath)
plot_oi
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-31-1.png" width="100%" />
As you can see, the CCL4L2-CCR5 interaction has no Omnipath DB scores. This is because this LR pair was not documented by the Omnipath LR database. Instead it was documented by the original NicheNet LR network as can be seen in the table.

Further note: Typically, there are way more than 50 differentially expressed and active ligand-receptor pairs per group across all sender-receiver combinations. Therefore it might be useful to zoom in on specific cell types as senders/receivers:

Eg M_Monocyte_CD16 as receiver:


```r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi, 
  receivers_oi = "M_Monocyte_CD16")
```


```r
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all))
plot_oi
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-33-1.png" width="100%" />

Eg M_Monocyte_CD16 as sender:


```r
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi, 
  senders_oi = "M_Monocyte_CD16")
```


```r
plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50 %>% inner_join(lr_network_all))
plot_oi
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-35-1.png" width="100%" />

You can make these plots also for the other condition(s), like we will illustrate now for the S group


```r
group_oi = "S"
```


```r
prioritized_tbl_oi_S_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_S_50 %>% inner_join(lr_network_all))
plot_oi
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-37-1.png" width="100%" />

## Visualization of differential ligand-target links

### Without filtering of target genes based on LR-target expression correlation

In another type of plot, we can visualize the ligand activities for a group-receiver combination, and show the predicted ligand-target links, and also the expression of the predicted target genes across samples.

For this, we now need to define a receiver cell type of interest. As example, we will take `L_T_TIM3._CD38._HLADR.` cells as receiver, and look at the top 10 senderLigand-receiverReceptor pairs with these cells as receiver.


```r
group_oi = "M"
receiver_oi = "M_Monocyte_CD16"
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  10, 
  groups_oi = group_oi, 
  receivers_oi = receiver_oi)
```


```r
combined_plot = make_ligand_activity_target_plot(
  group_oi, 
  receiver_oi, 
  prioritized_tbl_oi_M_10,
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix, 
  plot_legend = FALSE)
combined_plot
## $combined_plot
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-39-1.png" width="100%" />

```
## 
## $legends
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-39-2.png" width="100%" />


### With filtering of target genes based on LR-target expression correlation

In the previous plots, target genes were shown that are predicted as target gene of ligands based on prior knowledge. However, we can use the multi-sample nature of this data to filter target genes based on expression correlation between the upstream ligand-receptor pair and the downstream target gene. We will filter out correlated ligand-receptor --> target links that both show high expression correlation (spearman or pearson correlation > 0.50 in this example) and have some prior knowledge to support their link. Note that you can only make these visualization if you ran step 7 of the core MultiNicheNet analysis.


```r
group_oi = "M"
receiver_oi = "M_Monocyte_CD16"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
  inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% 
      distinct(ligand, target, direction_regulation, contrast)
    ) %>% 
  inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)

lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "up") %>% 
  filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "down") %>% 
  filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(
  lr_target_prior_cor_filtered_up, 
  lr_target_prior_cor_filtered_down)
```

Now we will visualize the top correlated target genes for the LR pairs that are also in the top 50 LR pairs discriminating the groups from each other:


```r
prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50, 
  groups_oi = group_oi, 
  receivers_oi = receiver_oi)
```


```r
lr_target_correlation_plot = make_lr_target_correlation_plot(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi,  
  lr_target_prior_cor_filtered , 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  receiver_oi,
  plot_legend = FALSE)
lr_target_correlation_plot$combined_plot
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-42-1.png" width="100%" />

You can also visualize the expression correlation in the following way for a selected LR pair and their targets:


```r
ligand_oi = "IFNG"
receptor_oi = "IFNGR2"
sender_oi = "L_T_TIM3._CD38._HLADR."
receiver_oi = "M_Monocyte_CD16"
lr_target_scatter_plot = make_lr_target_scatter_plot(
  multinichenet_output$prioritization_tables, 
  ligand_oi, receptor_oi, sender_oi, receiver_oi, 
  multinichenet_output$celltype_info, 
  multinichenet_output$grouping_tbl, 
  lr_target_prior_cor_filtered)
lr_target_scatter_plot
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-43-1.png" width="100%" />

## Intercellular regulatory network inference and visualization

In the plots before, we demonstrated that some DE genes have both expression correlation and prior knowledge support to be downstream of ligand-receptor pairs. Interestingly, some target genes can be ligands or receptors themselves. This illustrates that cells can send signals to other cells, who as a response to these signals produce signals themselves to feedback to the original sender cells, or who will effect other cell types. 

As last plot, we can generate a 'systems' view of these intercellular feedback and cascade processes than can be occuring between the different cell populations involved. In this plot, we will draw links between ligands of sender cell types their ligand/receptor-annotated target genes in receiver cell types. So links are ligand-target links (= gene regulatory links) and not ligand-receptor protein-protein interactions! We will infer this intercellular regulatory network here for the top100 interactions.(In practice, you can increase this number).


```r
prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  100, 
  rank_per_group = FALSE)

lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
        ) %>% 
      inner_join(contrast_tbl) %>% filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "up") %>% 
      filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "down") %>% 
      filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
      )
}) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 
```


```r
network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi)
network$links %>% head()
## # A tibble: 6 × 6
##   sender_ligand               receiver_target         direction_regulation group type  weight
##   <chr>                       <chr>                   <fct>                <chr> <chr>  <dbl>
## 1 M_Monocyte_CD16_HLA.A       L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
## 2 M_Monocyte_CD16_HLA.DRA     L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
## 3 M_Monocyte_CD16_LGALS1      L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
## 4 M_Monocyte_CD16_HLA.E       L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
## 5 L_T_TIM3._CD38._HLADR._CCL5 L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
## 6 M_Monocyte_CD16_HLA.DRB1    L_T_TIM3._CD38._HLADR.… up                   M     Liga…      1
network$nodes %>% head()
## # A tibble: 6 × 4
##   node                    celltype        gene    type_gene      
##   <chr>                   <chr>           <chr>   <chr>          
## 1 M_Monocyte_CD16_SIGLEC9 M_Monocyte_CD16 SIGLEC9 ligand/receptor
## 2 M_Monocyte_CD16_ENG     M_Monocyte_CD16 ENG     ligand/receptor
## 3 M_Monocyte_CD16_HLA.A   M_Monocyte_CD16 HLA.A   ligand         
## 4 M_Monocyte_CD16_HLA.DRA M_Monocyte_CD16 HLA.DRA ligand         
## 5 M_Monocyte_CD16_LGALS1  M_Monocyte_CD16 LGALS1  ligand         
## 6 M_Monocyte_CD16_HLA.E   M_Monocyte_CD16 HLA.E   ligand
```


```r
colors_sender["L_T_TIM3._CD38._HLADR."] = "pink" # the  original yellow background with white font is not very readable
network_graph = visualize_network(network, colors_sender)
network_graph$plot
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-46-1.png" width="100%" />

Interestingly, we can also use this network to further prioritize differential CCC interactions. Here we will assume that the most important LR interactions are the ones that are involved in this intercellular regulatory network. We can get these interactions as follows:

```r
network$prioritized_lr_interactions
## # A tibble: 60 × 5
##    group sender                 receiver               ligand   receptor
##    <chr> <chr>                  <chr>                  <chr>    <chr>   
##  1 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.A    CD8A    
##  2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA  LAG3    
##  3 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. LGALS1   CD69    
##  4 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.E    CD8A    
##  5 M     L_T_TIM3._CD38._HLADR. L_T_TIM3._CD38._HLADR. CCL5     CXCR3   
##  6 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRB1 LAG3    
##  7 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DPB1 LAG3    
##  8 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRB5 LAG3    
##  9 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. S100A8   CD69    
## 10 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. LGALS3   LAG3    
## # ℹ 50 more rows
```


```r
prioritized_tbl_oi_network = prioritized_tbl_oi %>% inner_join(
  network$prioritized_lr_interactions)
prioritized_tbl_oi_network
## # A tibble: 60 × 8
##    group sender       receiver ligand receptor id    prioritization_score prioritization_rank
##    <chr> <chr>        <chr>    <chr>  <chr>    <chr>                <dbl>               <dbl>
##  1 M     M_Monocyte_… L_T_TIM… S100A8 CD69     S100…                0.943                   1
##  2 M     L_T_TIM3._C… M_Monoc… IFNG   IFNGR1   IFNG…                0.940                   2
##  3 M     M_Monocyte_… M_Monoc… TIMP1  CD63     TIMP…                0.930                   3
##  4 M     L_T_TIM3._C… M_Monoc… IFNG   IFNGR2   IFNG…                0.930                   4
##  5 M     M_Monocyte_… L_T_TIM… HLA.D… LAG3     HLA.…                0.929                   5
##  6 M     M_Monocyte_… M_Monoc… S100A9 CD68     S100…                0.920                   6
##  7 S     M_Monocyte_… M_Monoc… TGFBI  ITGA4    TGFB…                0.919                   8
##  8 S     M_Monocyte_… L_T_TIM… TGFBI  ITGB1    TGFB…                0.917                  10
##  9 M     M_Monocyte_… M_Monoc… S100A8 CD68     S100…                0.912                  11
## 10 M     L_NK_CD56._… L_T_TIM… CCL3   CCR5     CCL3…                0.910                  13
## # ℹ 50 more rows
```

Visualize now the expression and activity of these interactions for the M group

```r
group_oi = "M"
```


```r
prioritized_tbl_oi_M = prioritized_tbl_oi_network %>% filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M %>% inner_join(lr_network_all)
  )
plot_oi
```

<img src="pairwise_analysis_MISC_files/figure-html/unnamed-chunk-50-1.png" width="100%" />

This was the end of this one-vs-one comparison vignettes. 
Check out the README page of the package to go to other vignettes that may build further upon this vignette (e.g., vignette showing additional downstream visualizations, how to incorporate additional data modalities in the prioritization etc). 
