---
title: "MultiNicheNet analysis: Workflow for condition-specific cell types"
author: "Robin Browaeys"
date: "2023-06-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{MultiNicheNet analysis: Workflow for condition-specific cell types}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style type="text/css">
.smaller {
  font-size: 10px
}
</style>

<!-- github markdown built using
rmarkdown::render("vignettes/condition_specific_celltype_MISC.Rmd", clean = FALSE )
-->


In this vignette, you can learn how to perform a MultiNicheNet analysis when you have one or more condition-specific cell types in your data. The differential expression analysis step in MultiNicheNet is the core issue for analyzing condition-specific cell types. This is because condition-specific cell types, by their very nature, do not have a counterpart in other conditions against which differential expression can be measured.

As the best solution, we recommend - if possible and biologically sensible to annotate your cells at one level higher in the cell type annotation hierarchy (eg sub cell type can be seen as a different cell state). Hereby, you will not have this problem of condition-specific cell types anymore. 

A a second solution, we propose the procedure as showcased in this vignette. This consists of 3 consecutive workflows: 

*1) General Workflow: This is the regular MultiNicheNet approach. In the general workflow, condition-specific cell types are excluded from the analysis. This approach focuses on achieving the best prioritization for other cell types that are present across different conditions, allowing for a standard DE analysis. The exclusion of condition-specific cell types simplifies the analysis but at the expense of potentially missing out on valuable insights that these unique cell types could provide. Therefore you can go further with the next two workflows: 
*2) Workflow A (Sender Cell Types): This workflow adapts the analysis to include condition-specific cell types as sender cell types. It circumvents the need for DE analysis of the ligand by only focusing the ligand-prioritization on ligand cell type specificity and ligand activity. 
*3) Workflow B (Receiver Cell Types): Workflow B allows for the inclusion of condition-specific cell types as receiver cell types by forgoing the DE analysis of the receiver, meaning that we cannot consider receptor differential expresion and ligand activity. Instead, receptor prioritization is based solely on the receptor cell type specificity. 

As example expression data of interacting cells, we will here use scRNAseq data of immune cells in MIS-C patients and healthy siblings from this paper of Hoste et al.: [TIM3+ TRBV11-2 T cells and IFNγ signature in patrolling monocytes and CD16+ NK cells delineate MIS-C](https://rupress.org/jem/article/219/2/e20211381/212918/TIM3-TRBV11-2-T-cells-and-IFN-signature-in) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6362434.svg)](https://doi.org/10.5281/zenodo.6362434)
. MIS-C (multisystem inflammatory syndrome in children) is a novel rare immunodysregulation syndrome that can arise after SARS-CoV-2 infection in children. We will use NicheNet to explore immune cell crosstalk enriched in MIS-C compared to healthy siblings and patients with adult COVID-19. 

# Preparation of the MultiNicheNet core analysis


```r
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
```

## Load NicheNet's ligand-receptor network and ligand-target matrix

MultiNicheNet builds upon the NicheNet framework and uses the same prior knowledge networks (ligand-receptor network and ligand-target matrix, currently v2 version).

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

## Read in SingleCellExperiment Objects 

In this vignette, we will load in a subset of the scRNAseq data of the MIS-C [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8010790.svg)](https://doi.org/10.5281/zenodo.8010790). For the sake of demonstration, this subset only contains 3 cell types. These celltypes are some of the cell types that were found to be most interesting related to MIS-C according to Hoste et al. 

If you start from a Seurat object, you can convert it easily to a SingleCellExperiment object via `sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")`.

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

In this case study, we want to study differences in cell-cell communication patterns between MIS-C patients (M), their healthy siblings (S) and adult patients with severe covid (A). The meta data columns that indicate this disease status(=group/condition of interest) is `MIS.C.AgeTier`. 

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

Here, we want to compare each patient group to the other groups, so the MIS-C (M) group vs healthy control siblings (S) and adult COVID19 patients (A) (= M vs S+A) and so on. We want to know which cell-cell communication patterns are specific for the M vs A+S group, the A vs M+S group and the S vs A+M group. 

To perform this comparison, we need to set the following contrasts:


```r
contrasts_oi = c("'M-(S+A)/2','S-(M+A)/2','A-(S+M)/2'")
```

__Very Important__ Note the format to indicate the contrasts! This formatting should be adhered to very strictly, and white spaces are not allowed! Check `?get_DE_info` for explanation about how to define this well. The most important points are that: 
*each contrast is surrounded by single quotation marks
*contrasts are separated by a comma without any white space 
*all contrasts together are surrounded by double quotation marks. 

If you compare against two groups, you should divide by 2 (as demonstrated here), if you compare against three groups, you should divide by 3 and so on.

For downstream visualizations and linking contrasts to their main condition, we also need to run the following:
This is necessary because we will also calculate cell-type+condition specificity of ligands and receptors. 


```r
contrast_tbl = tibble(contrast = 
                        c("M-(S+A)/2","S-(M+A)/2", "A-(S+M)/2"), 
                      group = c("M","S","A"))
```

If you want to compare only two groups (eg M vs S), you can use the following:
`contrasts_oi = c("'M-S','S-M'") `
`contrast_tbl = tibble(contrast = c("M-S","S-M"), group = c("M","S"))`

Other vignettes will demonstrate how to formalize different types of research questions.

### Define the sender and receiver cell types of interest.

If you want to focus the analysis on specific cell types (e.g. because you know which cell types reside in the same microenvironments based on spatial data), you can define this here. If you have sufficient computational resources and no specific idea of cell-type colocalzations, we recommend to consider all cell types as potential senders and receivers. Later on during analysis of the output it is still possible to zoom in on the cell types that interest you most, but your analysis is not biased to them.


```r
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
```

If the user wants it, it is possible to use only a subset of senders and receivers. Senders and receivers can be entirely different, but also overlapping, or the same. If you don't use all the cell types in your data, we recommend to continue with a subset of your data.


```r
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
```

In this vignette we will specifically demonstrate our workaround in case you are dealing with condition-specific cell types. These are cell types that are lacking in at least one condition, and are present in at least one other condition. 
For demonstration purposes, we will here artificially remove L_T_TIM3._CD38._HLADR. from the sibling and adult samples, making these cells specific for the MIS-C group. 


```r
sce_T = sce[, SummarizedExperiment::colData(sce)[,celltype_id] == "L_T_TIM3._CD38._HLADR."]
sce_T_Sibling =  sce_T[, SummarizedExperiment::colData(sce_T)[,group_id] %in% c("S","A")]
sce = sce[, setdiff(colnames(sce), colnames(sce_T_Sibling))]
```

# Running the MultiNicheNet core analysis

Now we will run the core of a MultiNicheNet analysis. This analysis consists of the following steps:

* 1. Cell-type filtering: determine which cell types are sufficiently present
* 2. Gene filtering: determine which genes are sufficiently expressed in each present cell type
* 3. Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type
* 4. Differential expression (DE) analysis: determine which genes are differentially expressed
* 5. Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes
* 6. Prioritization: rank cell-cell communication patterns through multi-criteria prioritization

Following these steps, one can optionally 
* 7. Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
* 8. Prioritize communication patterns involving condition-specific cell types through an alternative prioritization scheme

After these steps, the output can be further explored as we will demonstrate in the "Downstream analysis of the MultiNicheNet output" section. 

In this vignette, we will demonstrate these steps one-by-one, which offers the most flexibility to the user to assess intermediary results. 

## Cell-type filtering: determine which cell types are sufficiently present

In this step we will calculate and visualize cell type abundances. This will give an indication about which cell types will be retained in the analysis, and which cell types will be filtered out.  

Since MultiNicheNet will infer group differences at the sample level for each cell type (currently via Muscat - pseudobulking + EdgeR), we need to have sufficient cells per sample of a cell type, and this for all groups. In the following analysis we will set this minimum number of cells per cell type per sample at 10. Samples that have less than `min_cells` cells will be excluded from the analysis for that specific cell type.
 

```r
min_cells = 10
```

We recommend using `min_cells = 10`, except for datasets with several lowly abundant cell types of interest. For those datasets, we recommend using `min_cells = 5`.


```r
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
  )
```

First, we will check the cell type abundance diagnostic plots.

### Interpretation of cell type abundance information

The first plot visualizes the number of cells per celltype-sample combination, and indicates which combinations are removed during the DE analysis because there are less than `min_cells` in the celltype-sample combination. 


```r
abundance_info$abund_plot_sample
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-14-1.png)<!-- -->
The red dotted line indicates the required minimum of cells as defined above in `min_cells`. We can see here that L_T_TIM3._CD38._HLADR. are absent in the A and S groups, and only present in the M-group.

### Cell type filtering based on cell type abundance information

Running the following block of code can help you automatically determine which cell types are condition-specific and which cell types are absent. 


```r
sample_group_celltype_df = abundance_info$abundance_data %>% 
  filter(n > min_cells) %>% 
  ungroup() %>% 
  distinct(sample_id, group_id) %>% 
  cross_join(
    abundance_info$abundance_data %>% 
      ungroup() %>% 
      distinct(celltype_id)
    ) %>% 
  arrange(sample_id)

abundance_df = sample_group_celltype_df %>% left_join(
  abundance_info$abundance_data %>% ungroup()
  )

abundance_df$n[is.na(abundance_df$n)] = 0
abundance_df$keep[is.na(abundance_df$keep)] = FALSE
abundance_df_summarized = abundance_df %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)
  
print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## [1] "L_T_TIM3._CD38._HLADR."
  
print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)
## character(0)
```
Absent cell types will be filtered out, condition-specific cell types can be filtered out if you as a user do not want to run the alternative workflow for condition-specific cell types in the optional step 8 of the core MultiNicheNet analysis. However, in this analysis, we want to assess cell-cell communication involving condition-specific cell types, so we will keep condition-specific cell types by setting `analyse_condition_specific_celltypes = TRUE`


```r
analyse_condition_specific_celltypes = TRUE
```


```r
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]
```

## Gene filtering: determine which genes are sufficiently expressed in each present cell type

Before running the DE analysis, we will determine which genes are not sufficiently expressed and should be filtered out. 
We will perform gene filtering based on a similar procedure as used in `edgeR::filterByExpr`. However, we adapted this procedure to be more interpretable for single-cell datasets.  

For each cell type, we will consider genes expressed if they are expressed in at least a `min_sample_prop` fraction of samples in the condition with the lowest number of samples. By default, we set `min_sample_prop = 0.50`, which means that genes should be expressed in at least 2 samples if the group with lowest nr. of samples has 4 samples like this dataset. 


```r
min_sample_prop = 0.50
```

But how do we define which genes are expressed in a sample? For this we will consider genes as expressed if they have non-zero expression values in a `fraction_cutoff` fraction of cells of that cell type in that sample. By default, we set `fraction_cutoff = 0.05`, which means that genes should show non-zero expression values in at least 5% of cells in a sample. 


```r
fraction_cutoff = 0.05
```

We recommend using these default values unless there is specific interest in prioritizing (very) weakly expressed interactions. In that case, you could lower the value of `fraction_cutoff`. We explicitly recommend against using `fraction_cutoff > 0.10`.

Now we will calculate the information required for gene filtering with the following command:


```r
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
## [1] "Samples are considered if they have more than 10 cells of the cell type of interest"
## [1] "Genes with non-zero counts in at least 5% of cells of a cell type of interest in a particular sample will be considered as expressed in that sample."
## [1] "Genes expressed in at least 2 samples will considered as expressed in the cell type: L_NK_CD56._CD16."
## [1] "Genes expressed in at least 3.5 samples will considered as expressed in the cell type: L_T_TIM3._CD38._HLADR."
## [1] "Genes expressed in at least 2 samples will considered as expressed in the cell type: M_Monocyte_CD16"
## [1] "6621 genes are considered as expressed in the cell type: L_NK_CD56._CD16."
## [1] "7518 genes are considered as expressed in the cell type: L_T_TIM3._CD38._HLADR."
## [1] "8817 genes are considered as expressed in the cell type: M_Monocyte_CD16"
```

Now only keep genes that are expressed by at least one cell type:


```r
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]
```

## Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type

After filtering out absent cell types and genes, we will continue the analysis by calculating the different prioritization criteria that we will use to prioritize cell-cell communication patterns.

First, we will determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type. The function `process_abundance_expression_info` will link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types. This will later on allow us to define the cell-type specicificy criteria for ligands and receptors.


```r
abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
```

Normalized pseudobulk expression values per gene/celltype/sample can be inspected by:


```r
abundance_expression_info$celltype_info$pb_df %>% head()
## # A tibble: 6 × 4
##   gene  sample pb_sample celltype              
##   <chr> <chr>      <dbl> <fct>                 
## 1 A1BG  M1          4.70 L_T_TIM3._CD38._HLADR.
## 2 AAAS  M1          4.87 L_T_TIM3._CD38._HLADR.
## 3 AAGAB M1          5.33 L_T_TIM3._CD38._HLADR.
## 4 AAK1  M1          6.67 L_T_TIM3._CD38._HLADR.
## 5 AAMDC M1          3.46 L_T_TIM3._CD38._HLADR.
## 6 AAMP  M1          5.95 L_T_TIM3._CD38._HLADR.
```

An average of these sample-level expression values per condition/group can be inspected by:


```r
abundance_expression_info$celltype_info$pb_df_group %>% head()
## # A tibble: 6 × 4
## # Groups:   group, celltype [1]
##   group celltype         gene  pb_group
##   <chr> <chr>            <chr>    <dbl>
## 1 A     L_NK_CD56._CD16. A1BG      3.92
## 2 A     L_NK_CD56._CD16. AAAS      4.82
## 3 A     L_NK_CD56._CD16. AAGAB     4.83
## 4 A     L_NK_CD56._CD16. AAK1      7.04
## 5 A     L_NK_CD56._CD16. AAMDC     3.66
## 6 A     L_NK_CD56._CD16. AAMP      6.05
```

Inspecting these values for ligand-receptor interactions can be done by:


```r
abundance_expression_info$sender_receiver_info$pb_df %>% head()
## # A tibble: 6 × 8
##   sample sender                 receiver        ligand receptor pb_ligand pb_receptor ligand_receptor_pb_prod
##   <chr>  <chr>                  <chr>           <chr>  <chr>        <dbl>       <dbl>                   <dbl>
## 1 M4     M_Monocyte_CD16        M_Monocyte_CD16 B2M    LILRB1        14.1        11.4                    161.
## 2 M4     L_NK_CD56._CD16.       M_Monocyte_CD16 B2M    LILRB1        14.0        11.4                    160.
## 3 M4     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 B2M    LILRB1        13.8        11.4                    157.
## 4 M5     L_NK_CD56._CD16.       M_Monocyte_CD16 B2M    LILRB1        14.5        10.5                    153.
## 5 M5     M_Monocyte_CD16        M_Monocyte_CD16 B2M    LILRB1        14.4        10.5                    151.
## 6 M5     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 B2M    LILRB1        14.1        10.5                    149.
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()
## # A tibble: 6 × 8
## # Groups:   group, sender [5]
##   group sender                 receiver         ligand receptor pb_ligand_group pb_receptor_group ligand_receptor_pb_prod_group
##   <chr> <chr>                  <chr>            <chr>  <chr>              <dbl>             <dbl>                         <dbl>
## 1 M     L_NK_CD56._CD16.       M_Monocyte_CD16  B2M    LILRB1              14.1              9.99                          141.
## 2 M     M_Monocyte_CD16        M_Monocyte_CD16  B2M    LILRB1              14.1              9.99                          140.
## 3 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16  B2M    LILRB1              13.9              9.99                          139.
## 4 A     L_NK_CD56._CD16.       L_NK_CD56._CD16. B2M    KLRD1               14.3              9.72                          139.
## 5 S     L_NK_CD56._CD16.       L_NK_CD56._CD16. B2M    KLRD1               14.4              9.46                          136.
## 6 M     L_NK_CD56._CD16.       L_NK_CD56._CD16. B2M    KLRD1               14.1              9.49                          134.
```

## Differential expression (DE) analysis: determine which genes are differentially expressed

In this step, we will perform genome-wide differential expression analysis of receiver and sender cell types to define DE genes between the conditions of interest (as formalized by the `contrasts_oi`). Based on this analysis, we later can define the levels of differential expression of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver cell types (which will be used for the ligand activity analysis).

We will apply pseudobulking followed by EdgeR to perform multi-condition multi-sample differential expression (DE) analysis (also called 'differential state' analysis by the developers of Muscat). 


```r
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
## Error in perform_muscat_de_analysis(sce = sce_oi, sample_id = sample_id, : conditions written in contrasts should be in the condition-indicating column! This is not the case, which can lead to errors downstream.
```

Logically, DE analysis cannot be performed for the condition-specific cell types. This results here in an error for the "L_T_TIM3._CD38._HLADR." cell type. 

### Check DE results

Table with logFC and p-values for each gene-celltype-contrast:


```r
DE_info$celltype_de$de_output_tidy %>% head()
## # A tibble: 6 × 9
##   gene  cluster_id           logFC logCPM         F   p_val p_adj.loc p_adj contrast 
##   <chr> <chr>                <dbl>  <dbl>     <dbl>   <dbl>     <dbl> <dbl> <chr>    
## 1 A1BG  L_NK_CD56._CD16. -0.000822   4.04 0.0000126 0.997       0.999 0.999 M-(S+A)/2
## 2 AAAS  L_NK_CD56._CD16. -0.36       4.76 4.27      0.0529      0.457 0.457 M-(S+A)/2
## 3 AAGAB L_NK_CD56._CD16. -0.029      5    0.0254    0.875       0.999 0.999 M-(S+A)/2
## 4 AAK1  L_NK_CD56._CD16. -0.432      6.97 9.28      0.00612     0.217 0.217 M-(S+A)/2
## 5 AAMP  L_NK_CD56._CD16.  0.00309    6.13 0.000483  0.983       0.999 0.999 M-(S+A)/2
## 6 AAR2  L_NK_CD56._CD16.  0.251      4.78 2.58      0.139       0.61  0.61  M-(S+A)/2
```

Evaluate the distributions of p-values:


```r
DE_info$hist_pvals
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

These distributions look fine (uniform distribution, except peak at p-value <= 0.05), so we will continue using these regular p-values. In case these p-value distributions look irregular, you can estimate empirical p-values as we will demonstrate in another vignette.


```r
empirical_pval = FALSE
```


```r
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 
```

### Combine DE information for ligand-senders and receptors-receivers

To end this step, we will combine the DE information of senders and receivers by linking their ligands and receptors together based on the prior knowledge ligand-receptor network.


```r
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
```


```r
sender_receiver_de %>% head(20)
## # A tibble: 20 × 12
##    contrast  sender           receiver         ligand  receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor
##    <chr>     <chr>            <chr>            <chr>   <chr>         <dbl>        <dbl>                   <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
##  1 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  C3      VSIG4        2.96         5.76                      4.36      0.0964         0.439      0.00154           0.0923
##  2 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  IL10    IL10RB       5.61         0.715                     3.16      0.0052         0.142      0.0123            0.202 
##  3 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. IL10    IL10RA       5.61         0.0393                    2.82      0.0052         0.142      0.814             0.996 
##  4 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. IL10    IL10RB       5.61        -0.0727                    2.77      0.0052         0.142      0.78              0.987 
##  5 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  IL10    IL10RA       5.61        -0.316                     2.65      0.0052         0.142      0.233             0.621 
##  6 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  IL1B    IL1RAP       1.53         3.56                      2.54      0.0145         0.216      0.00294           0.114 
##  7 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  THBS1   CD47         4.58         0.446                     2.51      0.00784        0.166      0.0682            0.387 
##  8 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  THBS1   ITGA6        4.58         0.382                     2.48      0.00784        0.166      0.527             0.83  
##  9 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. THBS1   ITGA6        4.58         0.352                     2.47      0.00784        0.166      0.0844            0.524 
## 10 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. THBS1   ITGA4        4.58         0.283                     2.43      0.00784        0.166      0.0206            0.341 
## 11 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  HMGB1   CD163        0.154        4.54                      2.35      0.287          0.673      0.00403           0.127 
## 12 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. THBS1   CD47         4.58         0.0384                    2.31      0.00784        0.166      0.704             0.965 
## 13 M-(S+A)/2 L_NK_CD56._CD16. M_Monocyte_CD16  HMGB1   CD163        0.0745       4.54                      2.31      0.445          0.857      0.00403           0.127 
## 14 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  THBS1   ITGA4        4.58        -0.228                     2.18      0.00784        0.166      0.179             0.556 
## 15 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  C3      C3AR1        2.96         1.38                      2.17      0.0964         0.439      0.000331          0.0553
## 16 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  THBS1   LRP1         4.58        -0.28                      2.15      0.00784        0.166      0.129             0.492 
## 17 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. THBS1   ITGB1        4.58        -0.328                     2.13      0.00784        0.166      0.0774            0.508 
## 18 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  TNFSF12 CD163       -0.346        4.54                      2.10      0.515          0.822      0.00403           0.127 
## 19 M-(S+A)/2 M_Monocyte_CD16  M_Monocyte_CD16  THBS1   ITGB1        4.58        -0.412                     2.08      0.00784        0.166      0.0187            0.236 
## 20 M-(S+A)/2 M_Monocyte_CD16  L_NK_CD56._CD16. C3      ITGAX        2.96         0.751                     1.86      0.0964         0.439      0.0000237         0.0358
```

## Ligand activity prediction: use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes

In this step, we will predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results. We do this to prioritize interactions based on their predicted effect on a receiver cell type. We will assume that the most important group-specific interactions are those that lead to group-specific gene expression changes in a receiver cell type.

Similarly to base NicheNet (https://github.com/saeyslab/nichenetr), we use the DE output to create a "geneset of interest": here we assume that DE genes within a cell type may be DE because of differential cell-cell communication processes. In the ligand activity prediction, we will assess the enrichment of target genes of ligands within this geneset of interest. In case high-probabiliy target genes of a ligand are enriched in this set compared to the background of expressed genes, we predict that this ligand may have a high activity. 

Because the ligand activity analysis is an enrichment procedure, it is important that this geneset of interest should contain a sufficient but not too large number of genes. The ratio geneset_oi/background should ideally be between 1/200 and 1/10 (or close to these ratios).

To determine the genesets of interest based on DE output, we need to define some logFC and/or p-value thresholds per cell type/contrast combination. In general, we recommend inspecting the nr. of DE genes for all cell types based on the default thresholds and adapting accordingly. By default, we will apply the p-value cutoff on the normal p-values, and not on the p-values corrected for multiple testing. This choice was made because most multi-sample single-cell transcriptomics datasets have just a few samples per group and we might have a lack of statistical power due to pseudobulking. But, if the smallest group >= 20 samples, we typically recommend using p_val_adj = TRUE. When the biological difference between the conditions is very large, we typically recommend increasing the logFC_threshold and/or using p_val_adj = TRUE.

### Assess geneset_oi-vs-background ratios for different DE output tresholds prior to the NicheNet ligand activity analysis 

We will first inspect the geneset_oi-vs-background ratios for the default tresholds:


```r
logFC_threshold = 0.50
p_val_threshold = 0.05
```


```r
p_val_adj = FALSE 
```


```r
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment
## # A tibble: 6 × 12
##   cluster_id       n_background n_geneset_up n_geneset_down prop_geneset_up prop_geneset_down in_range_up in_range_down contrast  logFC_threshold p_val_threshold adjusted
##   <chr>                   <int>        <int>          <int>           <dbl>             <dbl> <lgl>       <lgl>         <chr>               <dbl>           <dbl> <lgl>   
## 1 L_NK_CD56._CD16.         6621          162             82          0.0245            0.0124 TRUE        TRUE          M-(S+A)/2             0.5            0.05 FALSE   
## 2 M_Monocyte_CD16          8817          647            438          0.0734            0.0497 TRUE        TRUE          M-(S+A)/2             0.5            0.05 FALSE   
## 3 L_NK_CD56._CD16.         6621          150            219          0.0227            0.0331 TRUE        TRUE          S-(M+A)/2             0.5            0.05 FALSE   
## 4 M_Monocyte_CD16          8817          368            254          0.0417            0.0288 TRUE        TRUE          S-(M+A)/2             0.5            0.05 FALSE   
## 5 L_NK_CD56._CD16.         6621          118            110          0.0178            0.0166 TRUE        TRUE          A-(S+M)/2             0.5            0.05 FALSE   
## 6 M_Monocyte_CD16          8817          262            464          0.0297            0.0526 TRUE        TRUE          A-(S+M)/2             0.5            0.05 FALSE
```
We can see here that for all cell type / contrast combinations, all geneset/background ratio's are within the recommended range (`in_range_up` and `in_range_down` columns). When these geneset/background ratio's would not be within the recommended ranges, we should interpret ligand activity results for these cell types with more caution, or use different thresholds (for these or all cell types). 

For the sake of demonstration, we will also calculate these ratio's in case we would use the adjusted p-value as threshold.


```r
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
    ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval
## # A tibble: 6 × 12
##   cluster_id       n_background n_geneset_up n_geneset_down prop_geneset_up prop_geneset_down in_range_up in_range_down contrast  logFC_threshold p_val_threshold adjusted
##   <chr>                   <int>        <int>          <int>           <dbl>             <dbl> <lgl>       <lgl>         <chr>               <dbl>           <dbl> <lgl>   
## 1 L_NK_CD56._CD16.         6621            7              0        0.00106           0        FALSE       FALSE         M-(S+A)/2             0.5            0.05 TRUE    
## 2 M_Monocyte_CD16          8817           25             11        0.00284           0.00125  FALSE       FALSE         M-(S+A)/2             0.5            0.05 TRUE    
## 3 L_NK_CD56._CD16.         6621           28             50        0.00423           0.00755  FALSE       TRUE          S-(M+A)/2             0.5            0.05 TRUE    
## 4 M_Monocyte_CD16          8817           10             15        0.00113           0.00170  FALSE       FALSE         S-(M+A)/2             0.5            0.05 TRUE    
## 5 L_NK_CD56._CD16.         6621           36             19        0.00544           0.00287  TRUE        FALSE         A-(S+M)/2             0.5            0.05 TRUE    
## 6 M_Monocyte_CD16          8817            4              3        0.000454          0.000340 FALSE       FALSE         A-(S+M)/2             0.5            0.05 TRUE
```
We can see here that for most cell type / contrast combinations, the geneset/background ratio's are not within the recommended range. Therefore, we will proceed with the default tresholds for the ligand activity analysis

### Perform the ligand activity analysis and ligand-target inference

After the ligand activity prediction, we will also infer the predicted target genes of these ligands in each contrast. For this ligand-target inference procedure, we also need to select which top n of the predicted target genes will be considered (here: top 250 targets per ligand). This parameter will not affect the ligand activity predictions. It will only affect ligand-target visualizations and construction of the intercellular regulatory network during the downstream analysis. We recommend users to test other settings in case they would be interested in exploring fewer, but more confident target genes, or vice versa. 


```r
top_n_target = 250
```

The NicheNet ligand activity analysis can be run in parallel for each receiver cell type, by changing the number of cores as defined here. Using more cores will speed up the analysis at the cost of needing more memory. This is only recommended if you have many receiver cell types of interest. 


```r
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 
```

Running the ligand activity prediction will take some time (the more cell types and contrasts, the more time)


```r
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))
```

You can check the output of the ligand activity and ligand-target inference here:


```r
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
## # A tibble: 20 × 8
## # Groups:   receiver, contrast [1]
##    ligand activity contrast  target   ligand_target_weight receiver         direction_regulation activity_scaled
##    <chr>     <dbl> <chr>     <chr>                   <dbl> <chr>            <fct>                          <dbl>
##  1 A2M      0.0282 M-(S+A)/2 ACOT7                 0.00632 L_NK_CD56._CD16. up                             0.794
##  2 A2M      0.0282 M-(S+A)/2 AREG                  0.00638 L_NK_CD56._CD16. up                             0.794
##  3 A2M      0.0282 M-(S+A)/2 CD55                  0.00649 L_NK_CD56._CD16. up                             0.794
##  4 A2M      0.0282 M-(S+A)/2 CD74                  0.00634 L_NK_CD56._CD16. up                             0.794
##  5 A2M      0.0282 M-(S+A)/2 FKBP5                 0.00723 L_NK_CD56._CD16. up                             0.794
##  6 A2M      0.0282 M-(S+A)/2 FOS                   0.0146  L_NK_CD56._CD16. up                             0.794
##  7 A2M      0.0282 M-(S+A)/2 GADD45A               0.0110  L_NK_CD56._CD16. up                             0.794
##  8 A2M      0.0282 M-(S+A)/2 H2AC6                 0.00747 L_NK_CD56._CD16. up                             0.794
##  9 A2M      0.0282 M-(S+A)/2 H2BC12                0.00692 L_NK_CD56._CD16. up                             0.794
## 10 A2M      0.0282 M-(S+A)/2 ISG20                 0.00738 L_NK_CD56._CD16. up                             0.794
## 11 A2M      0.0282 M-(S+A)/2 LMNB1                 0.00699 L_NK_CD56._CD16. up                             0.794
## 12 A2M      0.0282 M-(S+A)/2 MYC                   0.0199  L_NK_CD56._CD16. up                             0.794
## 13 A2M      0.0282 M-(S+A)/2 NFKB1                 0.00895 L_NK_CD56._CD16. up                             0.794
## 14 A2M      0.0282 M-(S+A)/2 NFKBIZ                0.00661 L_NK_CD56._CD16. up                             0.794
## 15 A2M      0.0282 M-(S+A)/2 PPP1R15A              0.00776 L_NK_CD56._CD16. up                             0.794
## 16 A2M      0.0282 M-(S+A)/2 SLC1A5                0.00710 L_NK_CD56._CD16. up                             0.794
## 17 A2M      0.0282 M-(S+A)/2 SMAD3                 0.00776 L_NK_CD56._CD16. up                             0.794
## 18 A2M      0.0282 M-(S+A)/2 SOCS3                 0.00968 L_NK_CD56._CD16. up                             0.794
## 19 A2M      0.0282 M-(S+A)/2 TFRC                  0.00712 L_NK_CD56._CD16. up                             0.794
## 20 A2M      0.0282 M-(S+A)/2 WARS1                 0.00743 L_NK_CD56._CD16. up                             0.794
```

## Prioritization: rank cell-cell communication patterns through multi-criteria prioritization

In the previous steps, we calculated expression, differential expression and NicheNet ligand activity. In the final step, we will now combine all calculated information to rank all sender-ligand---receiver-receptor pairs according to group/condition specificity. We will use the following criteria to prioritize ligand-receptor interactions:

* Upregulation of the ligand in a sender cell type and/or upregulation of the receptor in a receiver cell type - in the condition of interest.
* Cell-type specific expression of the ligand in the sender cell type and receptor in the receiver cell type in the condition of interest (to mitigate the influence of upregulated but still relatively weakly expressed ligands/receptors). 
* Sufficiently high expression levels of ligand and receptor in many samples of the same group.
* High NicheNet ligand activity, to further prioritize ligand-receptor pairs based on their predicted effect of the ligand-receptor interaction on the gene expression in the receiver cell type. 

We will combine these prioritization criteria in a single aggregated prioritization score. In the default setting, we will weigh each of these criteria equally (`scenario = "regular"`). This setting is strongly recommended. However, we also provide some additional setting to accomodate different biological scenarios. The setting `scenario = "lower_DE"` halves the weight for DE criteria and doubles the weight for ligand activity. This is recommended in case your hypothesis is that the differential CCC patterns in your data are less likely to be driven by DE (eg in cases of differential migration into a niche). The setting `scenario = "no_frac_LR_expr"` ignores the criterion "Sufficiently high expression levels of ligand and receptor in many samples of the same group". This may be interesting for users that have data with a limited number of samples and don’t want to penalize interactions if they are not sufficiently expressed in some samples. 

Finally, we still need to make one choice. For NicheNet ligand activity we can choose to prioritize ligands that only induce upregulation of target genes (`ligand_activity_down = FALSE`) or can lead potentially lead to both up- and downregulation (`ligand_activity_down = TRUE`). The benefit of `ligand_activity_down = FALSE` is ease of interpretability: prioritized ligand-receptor pairs will be upregulated in the condition of interest, just like their target genes.  `ligand_activity_down = TRUE` can be harder to interpret because target genes of some interactions may be upregulated in the other conditions compared to the condition of interest. This is harder to interpret, but may help to pick up interactions that can also repress gene expression. 

Here we will choose for setting `ligand_activity_down = FALSE` and focus specifically on upregulating ligands.


```r
ligand_activity_down = FALSE
```


```r
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular", # all prioritization criteria will be weighted equally
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down
  ))
```

Check the output tables

First: group-based summary table


```r
prioritization_tables$group_prioritization_tbl %>% head(20)
## # A tibble: 20 × 18
##    contrast  group sender           receiver         ligand receptor lr_interaction id                 scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>            <chr>            <chr>  <chr>    <chr>          <chr>                          <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16  L_NK_CD56._CD16. HLA.E  KLRC1    HLA.E_KLRC1    HLA.E_KLRC1_M_Mon…             0.812                  0.983               0.951                  0.989               0.978            0.866               1.00                  1                    0.952 M        
##  2 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16  VCAN   TLR2     VCAN_TLR2      VCAN_TLR2_M_Monoc…             0.938                  0.805               0.947                  0.937               0.714            1.00                1.00                  1                    0.922 A        
##  3 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  S100A9 CD68     S100A9_CD68    S100A9_CD68_M_Mon…             0.886                  0.873               0.863                  0.987               0.704            1.00                1.00                  1                    0.918 M        
##  4 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16  VCAN   ITGB1    VCAN_ITGB1     VCAN_ITGB1_M_Mono…             0.938                  0.805               0.824                  0.896               0.714            1.00                1.00                  1                    0.908 A        
##  5 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  HLA.F  LILRB1   HLA.F_LILRB1   HLA.F_LILRB1_M_Mo…             0.885                  0.986               0.983                  0.995               0.642            1.00                1.00                  0.857                0.904 M        
##  6 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  TNF    LTBR     TNF_LTBR       TNF_LTBR_M_Monocy…             0.984                  0.948               0.904                  0.921               0.962            1.00                1.00                  0.571                0.902 M        
##  7 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16  S100A8 CD36     S100A8_CD36    S100A8_CD36_M_Mon…             0.827                  0.692               0.966                  0.872               0.716            1.00                1.00                  1                    0.900 A        
##  8 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  HLA.F  LILRB2   HLA.F_LILRB2   HLA.F_LILRB2_M_Mo…             0.885                  0.986               0.940                  0.980               0.642            1.00                1.00                  0.857                0.899 M        
##  9 M-(S+A)/2 M     M_Monocyte_CD16  L_NK_CD56._CD16. HLA.C  KIR2DL1  HLA.C_KIR2DL1  HLA.C_KIR2DL1_M_M…             0.893                  0.971               0.796                  0.898               0.761            0.833               1.00                  1                    0.896 M        
## 10 S-(M+A)/2 S     L_NK_CD56._CD16. M_Monocyte_CD16  TGFB1  ENG      TGFB1_ENG      TGFB1_ENG_L_NK_CD…             0.669                  0.866               0.972                  0.967               0.632            1.00                1.00                  1                    0.895 S        
## 11 M-(S+A)/2 M     L_NK_CD56._CD16. M_Monocyte_CD16  IFNG   IFNGR1   IFNG_IFNGR1    IFNG_IFNGR1_L_NK_…             0.826                  0.719               0.872                  0.891               1.00             1.00                1.00                  0.714                0.895 M        
## 12 S-(M+A)/2 S     L_NK_CD56._CD16. L_NK_CD56._CD16. TGFB1  TGFBR3   TGFB1_TGFBR3   TGFB1_TGFBR3_L_NK…             0.669                  0.866               0.746                  0.957               0.725            1.00                1.00                  1                    0.891 S        
## 13 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16  NAMPT  ITGB1    NAMPT_ITGB1    NAMPT_ITGB1_M_Mon…             0.857                  0.856               0.824                  0.896               0.616            1.00                1.00                  1                    0.889 A        
## 14 M-(S+A)/2 M     M_Monocyte_CD16  L_NK_CD56._CD16. HLA.F  KLRC1    HLA.F_KLRC1    HLA.F_KLRC1_M_Mon…             0.885                  0.986               0.951                  0.989               0.567            1.00                1.00                  0.857                0.889 M        
## 15 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  HEBP1  FPR2     HEBP1_FPR2     HEBP1_FPR2_M_Mono…             0.874                  0.950               0.988                  0.977               0.718            1.00                1.00                  0.714                0.888 M        
## 16 S-(M+A)/2 S     M_Monocyte_CD16  L_NK_CD56._CD16. GNAS   ADRB2    GNAS_ADRB2     GNAS_ADRB2_M_Mono…             0.823                  0.931               0.898                  0.990               0.503            1.00                1.00                  1                    0.888 S        
## 17 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16  HLA.C  LILRB1   HLA.C_LILRB1   HLA.C_LILRB1_M_Mo…             0.893                  0.971               0.983                  0.995               0.567            0.833               1.00                  1                    0.887 M        
## 18 A-(S+M)/2 A     L_NK_CD56._CD16. L_NK_CD56._CD16. HLA.B  KIR2DL3  HLA.B_KIR2DL3  HLA.B_KIR2DL3_L_N…             0.775                  0.900               0.899                  0.913               0.569            1.00                1.00                  1                    0.886 A        
## 19 M-(S+A)/2 M     M_Monocyte_CD16  L_NK_CD56._CD16. SELPLG SELL     SELPLG_SELL    SELPLG_SELL_M_Mon…             0.836                  0.942               0.879                  0.979               0.482            1.00                1.00                  1                    0.884 M        
## 20 S-(M+A)/2 S     M_Monocyte_CD16  M_Monocyte_CD16  TGFBI  ITGA4    TGFBI_ITGA4    TGFBI_ITGA4_M_Mon…             0.892                  0.943               0.832                  0.955               0.468            1.00                1.00                  1                    0.880 S        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
```
This table gives the final prioritization score of each interaction, and the values of the individual prioritization criteria.

With this step, all required steps are finished. Now, we can optionally still run the following steps
* Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
* Prioritize communication patterns involving condition-specific cell types through an alternative prioritization scheme

Here we will first focus on the expression correlation step:

## Calculate the across-samples expression correlation between ligand-receptor pairs and target genes

In multi-sample datasets, we have the opportunity to look whether expression of ligand-receptor across all samples is correlated with the expression of their by NicheNet predicted target genes. This is what we will do with the following line of code:


```r
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
  )
```

## Save all the output of MultiNicheNet 

To avoid needing to redo the analysis later, we will here to save an output object that contains all information to perform all downstream analyses. 

```r
path = "./"

multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
multinichenet_output = make_lite_output(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))

}
```

# Interpreting the MultiNicheNet analysis output

## Visualization of differential cell-cell interactions

### Summarizing ChordDiagram circos plots

In a first instance, we will look at the broad overview of prioritized interactions via condition-specific Chordiagram circos plots. The aim of this visualizatin is to provide a summary of the top prioritized senderLigand-receiverReceptor interactions per condition (between all cell types or between cell type pairs of interest). 

We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.


```r
prioritized_tbl_oi_all = get_top_n_lr_pairs(prioritization_tables, 50, rank_per_group = FALSE)
```


```r
prioritized_tbl_oi = prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0 
```


```r
senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-48-1.png)<!-- -->![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-48-2.png)<!-- -->![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-48-3.png)<!-- -->![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-48-4.png)<!-- -->

The above steps and prioritization are the same as for a regular MultiNicheNet analysis, and only focus on cell types present in all conditions. Now we will showcase two workflows that enable studying condition-specific cell types (**"Prioritize communication patterns involving condition-specific cell types through an alternative prioritization scheme"**).

The first one (workflow A) looks at condition-specific cell types as sender cell types. The second one (workflow B) looks at condition-specific cell types as receiver cell types.

# Condition-specificity workflow mode A: investigating the role of condition-specific cell types as SENDER cell types

## A-1: Add interactions potentially involving condition-specific cell types

Here we will demonstrate how to look whether the condition-specific cell type(s) may contribute to important cell-cell communication patterns in your system of interest. 
In MultiNicheNet, we prioritize ligand-receptor interactions based on: 1) ligand expression (condition/cell-type specific); 2) receptor expression (condition/cell-type specific); 3) ligand activity (downstream signature in receiver)

For condition-specific cell type(s) we cannot calculate DE, but we can still calculate cell-type specificity and their ligand activity on other cell types. 

We will thus assume that following ligands from condition-specific cell types may be important:
* ligands are cell-type specific for the condition-specific cell type(s); meaning that this cell type is potentially the predominant source of this ligand
* these ligands have strong downstream target gene enrichment (high ligand activity) 

Based on these criteria, we will add ligand-receptor pairs to the prioritization_table involving the condition-specific cell type as sender. 

Now we will prioritize the interactions with the function: `prioritize_condition_specific_sender`. This function will not use the logFC & p-value for DE of the ligand anymore. For the ligand, only celltype specific expression in the condition of interest, and ligand activity will be considered.


```r
prioritization_tables_with_condition_specific_celltype_sender = prioritize_condition_specific_sender(
  abundance_info = abundance_info,
  abundance_expression_info = abundance_expression_info, 
  condition_specific_celltypes = condition_specific_celltypes, 
  grouping_tbl = grouping_tbl, 
  fraction_cutoff = fraction_cutoff, 
  contrast_tbl = contrast_tbl, 
  sender_receiver_de = sender_receiver_de, 
  lr_network = lr_network, 
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  scenario = "regular",
  ligand_activity_down = FALSE
)
```


```r
prioritization_tables_with_condition_specific_celltype_sender$group_prioritization_tbl %>% head(20)
## # A tibble: 20 × 18
##    contrast  group sender                 receiver         ligand receptor lr_interaction id           scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>                  <chr>            <chr>  <chr>    <chr>          <chr>                    <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16        L_NK_CD56._CD16. HLA.E  KLRC1    HLA.E_KLRC1    HLA.E_KLRC1…           0.871                  0.988                 0.951                  0.989               0.978            0.866               1.00                  1                    0.963 M        
##  2 M-(S+A)/2 M     L_NK_CD56._CD16.       L_NK_CD56._CD16. HLA.E  KLRC1    HLA.E_KLRC1    HLA.E_KLRC1…           0.585                  0.464                 0.951                  0.989               0.978            0.794               1.00                  1                    0.948 M        
##  3 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16  IFNG   IFNGR1   IFNG_IFNGR1    IFNG_IFNGR1…           0.129                  0.129                 0.872                  0.891               1.00             1.00                1.00                  0.857                0.948 M        
##  4 A-(S+M)/2 A     L_NK_CD56._CD16.       L_NK_CD56._CD16. IFITM1 CD81     IFITM1_CD81    IFITM1_CD81…           0.694                  0.735                 0.686                  0.741               1.00             0.979               1.00                  1                    0.939 A        
##  5 A-(S+M)/2 A     L_NK_CD56._CD16.       M_Monocyte_CD16  TGFB1  ITGB1    TGFB1_ITGB1    TGFB1_ITGB1…           0.608                  0.550                 0.824                  0.896               1.00             0.817               1.00                  1                    0.936 A        
##  6 A-(S+M)/2 A     M_Monocyte_CD16        M_Monocyte_CD16  VCAN   TLR2     VCAN_TLR2      VCAN_TLR2_M…           0.958                  0.866                 0.947                  0.937               0.714            1.00                1.00                  1                    0.932 A        
##  7 A-(S+M)/2 A     M_Monocyte_CD16        M_Monocyte_CD16  S100A8 CD36     S100A8_CD36    S100A8_CD36…           0.881                  0.789                 0.966                  0.872               0.716            1.00                1.00                  1                    0.927 A        
##  8 M-(S+A)/2 M     M_Monocyte_CD16        M_Monocyte_CD16  S100A9 CD68     S100A9_CD68    S100A9_CD68…           0.921                  0.913                 0.863                  0.987               0.704            1.00                1.00                  1                    0.926 M        
##  9 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16  IFNG   IFNGR2   IFNG_IFNGR2    IFNG_IFNGR2…           0.129                  0.129                 0.701                  0.828               1.00             1.00                1.00                  0.857                0.925 M        
## 10 S-(M+A)/2 S     L_NK_CD56._CD16.       M_Monocyte_CD16  TGFB1  ENG      TGFB1_ENG      TGFB1_ENG_L…           0.772                  0.908                 0.972                  0.967               0.632            1.00                1.00                  1                    0.921 S        
## 11 S-(M+A)/2 S     L_NK_CD56._CD16.       L_NK_CD56._CD16. TGFB1  TGFBR3   TGFB1_TGFBR3   TGFB1_TGFBR…           0.772                  0.908                 0.746                  0.957               0.725            1.00                1.00                  1                    0.916 S        
## 12 A-(S+M)/2 A     M_Monocyte_CD16        M_Monocyte_CD16  VCAN   ITGB1    VCAN_ITGB1     VCAN_ITGB1_…           0.958                  0.866                 0.824                  0.896               0.714            1.00                1.00                  1                    0.915 A        
## 13 M-(S+A)/2 M     L_NK_CD56._CD16.       M_Monocyte_CD16  HLA.F  LILRB1   HLA.F_LILRB1   HLA.F_LILRB…           0.646                  0.650                 0.983                  0.995               0.642            0.934               1.00                  1                    0.913 M        
## 14 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16  SIRPG  CD47     SIRPG_CD47     SIRPG_CD47_…           0.00136                0.00136               0.834                  0.890               0.644            1.00                1.00                  1                    0.902 M        
## 15 S-(M+A)/2 S     M_Monocyte_CD16        M_Monocyte_CD16  LGALS9 CD44     LGALS9_CD44    LGALS9_CD44…           0.790                  0.825                 0.849                  0.962               0.584            1.00                1.00                  1                    0.898 S        
## 16 M-(S+A)/2 M     M_Monocyte_CD16        M_Monocyte_CD16  HLA.F  LILRB1   HLA.F_LILRB1   HLA.F_LILRB…           0.921                  0.991                 0.983                  0.995               0.642            1.00                1.00                  0.857                0.898 M        
## 17 M-(S+A)/2 M     L_NK_CD56._CD16.       M_Monocyte_CD16  IFNG   IFNGR1   IFNG_IFNGR1    IFNG_IFNGR1…           0.880                  0.807                 0.872                  0.891               1.00             0.890               1.00                  0.714                0.898 M        
## 18 A-(S+M)/2 A     M_Monocyte_CD16        M_Monocyte_CD16  TIMP1  CD63     TIMP1_CD63     TIMP1_CD63_…           0.839                  0.820                 0.964                  0.972               0.575            0.940               1.00                  1                    0.897 A        
## 19 S-(M+A)/2 S     M_Monocyte_CD16        L_NK_CD56._CD16. TNF    TNFRSF1A TNF_TNFRSF1A   TNF_TNFRSF1…           0.802                  0.699                 0.766                  0.824               0.701            0.985               1.00                  1                    0.896 S        
## 20 S-(M+A)/2 S     L_NK_CD56._CD16.       M_Monocyte_CD16  TGFB1  TGFBR2   TGFB1_TGFBR2   TGFB1_TGFBR…           0.772                  0.908                 0.792                  0.895               0.632            1.00                1.00                  1                    0.895 S        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
prioritization_tables_with_condition_specific_celltype_sender$group_prioritization_tbl %>% filter(sender %in% condition_specific_celltypes) %>% head(10)
## # A tibble: 10 × 18
##    contrast  group sender                 receiver        ligand   receptor lr_interaction  id         scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>                  <chr>           <chr>    <chr>    <chr>           <chr>                  <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 IFNG     IFNGR1   IFNG_IFNGR1     IFNG_IFNG…           0.129                  0.129                 0.872                  0.891               1.00              1.00              1.00                   0.857                0.948 M        
##  2 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 IFNG     IFNGR2   IFNG_IFNGR2     IFNG_IFNG…           0.129                  0.129                 0.701                  0.828               1.00              1.00              1.00                   0.857                0.925 M        
##  3 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 SIRPG    CD47     SIRPG_CD47      SIRPG_CD4…           0.00136                0.00136               0.834                  0.890               0.644             1.00              1.00                   1                    0.902 M        
##  4 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CD82     CD37     CD82_CD37       CD82_CD37…           0.0257                 0.0257                0.716                  0.900               0.653             1.00              1.00                   1                    0.893 M        
##  5 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CEACAM21 LILRA6   CEACAM21_LILRA6 CEACAM21_…           0.185                  0.185                 0.897                  0.798               0.705             1.00              1.00                   0.857                0.882 M        
##  6 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CD40LG   CD40     CD40LG_CD40     CD40LG_CD…           0.0115                 0.0115                0.911                  0.866               0.612             1.00              1.00                   0.857                0.872 M        
##  7 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CD82     ITGB2    CD82_ITGB2      CD82_ITGB…           0.0257                 0.0257                0.707                  0.827               0.653             1.00              0.925                  1                    0.869 M        
##  8 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CD2      CD53     CD2_CD53        CD2_CD53_…           0.110                  0.110                 0.920                  0.999               0.592             1.00              0.780                  1                    0.866 M        
##  9 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 SELL     SELPLG   SELL_SELPLG     SELL_SELP…           0.167                  0.167                 0.858                  0.948               0.404             1.00              1.00                   1                    0.862 M        
## 10 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16 CD40LG   ITGB2    CD40LG_ITGB2    CD40LG_IT…           0.0115                 0.0115                0.707                  0.827               0.612             1.00              0.925                  1                    0.861 M        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
```

We see that some top-predicted interactions are coming from the condition-specific cell type(s) (L_T_TIM3._CD38._HLADR.). This means that these ligands are the most strongly expressed by that cell type, and that they have quite high ligand activity.

# Condition-specificity workflow mode B: investigating the role of condition-specific cell types as RECEIVER cell types

## B-1: Add interactions potentially involving condition-specific cell types

For condition-specific receiver cell type(s) we cannot calculate DE, meaning that we can only focus on cell-type specificity of receptors as prioritization criteria since ligand activity calculation requires DE analysis. 

We will thus assume that following ligands from condition-specific cell types may be important:
* receptors are cell-type specific for the condition-specific cell type(s); meaning that this cell type is potentially the predominant "sink" of this ligand

Based on these criteria, we will add ligand-receptor pairs to the prioritization_table involving the condition-specific cell type as receiver 

Now we will prioritize the interactions with the function: `prioritize_condition_specific_receiver`. This function will not use the logFC & p-value for DE of the receptor anymore, nor ligand activity. For the receptor, only celltype specific expression in the condition of interest will be considered.


```r
prioritization_tables_with_condition_specific_celltype_receiver = prioritize_condition_specific_receiver(
  abundance_info = abundance_info,
  abundance_expression_info = abundance_expression_info, 
  condition_specific_celltypes = condition_specific_celltypes, 
  grouping_tbl = grouping_tbl, 
  fraction_cutoff = fraction_cutoff, 
  contrast_tbl = contrast_tbl, 
  sender_receiver_de = sender_receiver_de, 
  lr_network = lr_network, 
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  scenario = "regular",
  ligand_activity_down = FALSE
)
```


```r
prioritization_tables_with_condition_specific_celltype_receiver$group_prioritization_tbl %>% head(20)
## # A tibble: 20 × 18
##    contrast  group sender           receiver               ligand   receptor lr_interaction id         scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>            <chr>                  <chr>    <chr>    <chr>          <chr>                  <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16        HLA.DRA  CD37     HLA.DRA_CD37   HLA.DRA_C…             0.883                  0.978              0.807                  0.932                    0             1.00              1.00                       1                0.983 M        
##  2 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRA  CD82     HLA.DRA_CD82   HLA.DRA_C…             0.883                  0.978              0.0343                 0.0343                   0             1.00              1.00                       1                0.983 M        
##  3 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRA  LAG3     HLA.DRA_LAG3   HLA.DRA_L…             0.883                  0.978              0.0952                 0.0952                   0             1.00              1.00                       1                0.983 M        
##  4 M-(S+A)/2 M     L_NK_CD56._CD16. M_Monocyte_CD16        SELL     SELPLG   SELL_SELPLG    SELL_SELP…             0.866                  0.984              0.903                  0.965                    0             1.00              1.00                       1                0.982 M        
##  5 S-(M+A)/2 S     M_Monocyte_CD16  M_Monocyte_CD16        TGFBI    ITGA4    TGFBI_ITGA4    TGFBI_ITG…             0.892                  0.943              0.886                  0.969                    0             1.00              1.00                       1                0.980 S        
##  6 S-(M+A)/2 S     M_Monocyte_CD16  M_Monocyte_CD16        HLA.DQA1 CD4      HLA.DQA1_CD4   HLA.DQA1_…             0.938                  0.949              0.712                  0.742                    0             1.00              0.971                      1                0.979 S        
##  7 M-(S+A)/2 M     M_Monocyte_CD16  L_NK_CD56._CD16.       HLA.DRA  CD81     HLA.DRA_CD81   HLA.DRA_C…             0.883                  0.978              0.818                  0.856                    0             1.00              0.967                      1                0.975 M        
##  8 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. SELPLG   SELL     SELPLG_SELL    SELPLG_SE…             0.836                  0.942              0.139                  0.139                    0             1.00              1.00                       1                0.973 M        
##  9 M-(S+A)/2 M     M_Monocyte_CD16  M_Monocyte_CD16        S100A9   CD68     S100A9_CD68    S100A9_CD…             0.886                  0.873              0.907                  0.991                    0             1.00              1.00                       1                0.970 M        
## 10 S-(M+A)/2 S     M_Monocyte_CD16  L_NK_CD56._CD16.       GNAS     ADRB2    GNAS_ADRB2     GNAS_ADRB…             0.823                  0.931              0.931                  0.993                    0             1.00              1.00                       1                0.970 S        
## 11 S-(M+A)/2 S     M_Monocyte_CD16  L_NK_CD56._CD16.       GNAS     PTGDR    GNAS_PTGDR     GNAS_PTGD…             0.823                  0.931              0.870                  0.954                    0             1.00              1.00                       1                0.970 S        
## 12 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16        VCAN     TLR2     VCAN_TLR2      VCAN_TLR2…             0.938                  0.805              0.964                  0.957                    0             1.00              1.00                       1                0.968 A        
## 13 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16        VCAN     ITGB1    VCAN_ITGB1     VCAN_ITGB…             0.938                  0.805              0.880                  0.929                    0             1.00              1.00                       1                0.968 A        
## 14 A-(S+M)/2 A     L_NK_CD56._CD16. L_NK_CD56._CD16.       GZMB     IGF2R    GZMB_IGF2R     GZMB_IGF2…             0.832                  0.902              0.810                  0.889                    0             1.00              1.00                       1                0.967 A        
## 15 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16        NAMPT    ITGB1    NAMPT_ITGB1    NAMPT_ITG…             0.857                  0.856              0.880                  0.929                    0             1.00              1.00                       1                0.965 A        
## 16 S-(M+A)/2 S     M_Monocyte_CD16  L_NK_CD56._CD16.       PTPRC    CD247    PTPRC_CD247    PTPRC_CD2…             0.784                  0.944              0.744                  0.848                    0             1.00              0.992                      1                0.964 S        
## 17 A-(S+M)/2 A     L_NK_CD56._CD16. M_Monocyte_CD16        SPON2    ITGB1    SPON2_ITGB1    SPON2_ITG…             0.814                  0.892              0.880                  0.929                    0             1.00              1.00                       1                0.964 A        
## 18 A-(S+M)/2 A     M_Monocyte_CD16  M_Monocyte_CD16        LGALS3   ANXA2    LGALS3_ANXA2   LGALS3_AN…             0.915                  0.981              0.550                  0.493                    0             1.00              0.903                      1                0.963 A        
## 19 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. CD47     SIRPG    CD47_SIRPG     CD47_SIRP…             0.816                  0.876              0.0497                 0.0497                   0             1.00              1.00                       1                0.962 M        
## 20 A-(S+M)/2 A     L_NK_CD56._CD16. L_NK_CD56._CD16.       HLA.B    KIR3DL1  HLA.B_KIR3DL1  HLA.B_KIR…             0.775                  0.900              0.934                  0.880                    0             1.00              1.00                       1                0.960 A        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
prioritization_tables_with_condition_specific_celltype_receiver$group_prioritization_tbl %>% filter(receiver %in% condition_specific_celltypes) %>% head(10)
## # A tibble: 10 × 18
##    contrast  group sender           receiver               ligand   receptor lr_interaction id         scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>            <chr>                  <chr>    <chr>    <chr>          <chr>                  <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRA  CD82     HLA.DRA_CD82   HLA.DRA_C…             0.883                  0.978              0.0343                 0.0343                   0            1.00               1.00                   1                    0.983 M        
##  2 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRA  LAG3     HLA.DRA_LAG3   HLA.DRA_L…             0.883                  0.978              0.0952                 0.0952                   0            1.00               1.00                   1                    0.983 M        
##  3 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. SELPLG   SELL     SELPLG_SELL    SELPLG_SE…             0.836                  0.942              0.139                  0.139                    0            1.00               1.00                   1                    0.973 M        
##  4 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. CD47     SIRPG    CD47_SIRPG     CD47_SIRP…             0.816                  0.876              0.0497                 0.0497                   0            1.00               1.00                   1                    0.962 M        
##  5 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRB1 CD82     HLA.DRB1_CD82  HLA.DRB1_…             0.925                  0.996              0.0343                 0.0343                   0            1.00               1.00                   0.857                0.955 M        
##  6 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.DRB1 LAG3     HLA.DRB1_LAG3  HLA.DRB1_…             0.925                  0.996              0.0952                 0.0952                   0            1.00               1.00                   0.857                0.955 M        
##  7 M-(S+A)/2 M     L_NK_CD56._CD16. L_T_TIM3._CD38._HLADR. CCL3     CCR5     CCL3_CCR5      CCL3_CCR5…             0.933                  0.968              0.0532                 0.0532                   0            1.00               1.00                   0.857                0.952 M        
##  8 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.F    CD8A     HLA.F_CD8A     HLA.F_CD8…             0.885                  0.986              0.0399                 0.0399                   0            1.00               1.00                   0.857                0.949 M        
##  9 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. TNFSF13  FAS      TNFSF13_FAS    TNFSF13_F…             0.898                  0.947              0.0798                 0.0798                   0            1.00               0.999                  0.857                0.945 M        
## 10 M-(S+A)/2 M     M_Monocyte_CD16  L_T_TIM3._CD38._HLADR. HLA.E    CD8A     HLA.E_CD8A     HLA.E_CD8…             0.812                  0.983              0.0399                 0.0399                   0            0.866              1.00                   1                    0.941 M        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
```
We see that some top-predicted interactions are coming from the condition-specific cell type. This means that these receptors are the most strongly expressed by that cell type.

# Combine the condition-specific workflow with the regular one

Finally we can combine the prioritization scores of the regular workflow with the scores from the workflows focusing on the condition-specific cell types. 

This involves a strategic approach where:

Regular scores are retained for interactions involving cell types that are present across all conditions.
Workflow A scores are utilized for interactions that involve condition-specific cell types functioning as senders.
Workflow B scores are applied for interactions with condition-specific cell types in the receiver role.

Note that these scores are only semi-comparable because the scores are calculated from different combinations of prioritization criteria. We want to emphasize that you can only use these scores for visualizing all interactions together, not to make statements concerning which cell types are most likely to interact given the number of top-prioritized interactions between cell types. 


```r
combined_prioritization_tables = list(
  group_prioritization_tbl = bind_rows(
    prioritization_tables_with_condition_specific_celltype_receiver$group_prioritization_tbl %>% filter(receiver %in% condition_specific_celltypes),
    prioritization_tables_with_condition_specific_celltype_sender$group_prioritization_tbl %>% filter(sender %in% condition_specific_celltypes)
  ) %>% 
    bind_rows(
      prioritization_tables$group_prioritization_tbl
  ) %>% arrange(-prioritization_score) %>% distinct()

)
combined_prioritization_tables %>% head()
## $group_prioritization_tbl
## # A tibble: 4,940 × 18
##    contrast  group sender                 receiver               ligand  receptor lr_interaction id    scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>                  <chr>                  <chr>   <chr>    <chr>          <chr>             <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA CD82     HLA.DRA_CD82   HLA.…             0.883                  0.978              0.0343                 0.0343               0                1.00                1.00                  1                    0.983 M        
##  2 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA LAG3     HLA.DRA_LAG3   HLA.…             0.883                  0.978              0.0952                 0.0952               0                1.00                1.00                  1                    0.983 M        
##  3 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. SELPLG  SELL     SELPLG_SELL    SELP…             0.836                  0.942              0.139                  0.139                0                1.00                1.00                  1                    0.973 M        
##  4 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. CD47    SIRPG    CD47_SIRPG     CD47…             0.816                  0.876              0.0497                 0.0497               0                1.00                1.00                  1                    0.962 M        
##  5 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DR… CD82     HLA.DRB1_CD82  HLA.…             0.925                  0.996              0.0343                 0.0343               0                1.00                1.00                  0.857                0.955 M        
##  6 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DR… LAG3     HLA.DRB1_LAG3  HLA.…             0.925                  0.996              0.0952                 0.0952               0                1.00                1.00                  0.857                0.955 M        
##  7 M-(S+A)/2 M     L_NK_CD56._CD16.       L_T_TIM3._CD38._HLADR. CCL3    CCR5     CCL3_CCR5      CCL3…             0.933                  0.968              0.0532                 0.0532               0                1.00                1.00                  0.857                0.952 M        
##  8 M-(S+A)/2 M     M_Monocyte_CD16        L_NK_CD56._CD16.       HLA.E   KLRC1    HLA.E_KLRC1    HLA.…             0.812                  0.983              0.951                  0.989                0.978            0.866               1.00                  1                    0.952 M        
##  9 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.F   CD8A     HLA.F_CD8A     HLA.…             0.885                  0.986              0.0399                 0.0399               0                1.00                1.00                  0.857                0.949 M        
## 10 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16        IFNG    IFNGR1   IFNG_IFNGR1    IFNG…             0.129                  0.129              0.872                  0.891                1.00             1.00                1.00                  0.857                0.948 M        
## # ℹ 4,930 more rows
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
combined_prioritization_tables$group_prioritization_tbl %>% filter(sender %in% condition_specific_celltypes  | receiver %in% condition_specific_celltypes) %>% head(10)
## # A tibble: 10 × 18
##    contrast  group sender                 receiver               ligand  receptor lr_interaction id    scaled_lfc_ligand scaled_p_val_ligand_…¹ scaled_lfc_receptor scaled_p_val_recepto…² max_scaled_activity scaled_pb_ligand scaled_pb_receptor fraction_expressing_…³ prioritization_score top_group
##    <chr>     <chr> <chr>                  <chr>                  <chr>   <chr>    <chr>          <chr>             <dbl>                  <dbl>               <dbl>                  <dbl>               <dbl>            <dbl>              <dbl>                  <dbl>                <dbl> <chr>    
##  1 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA CD82     HLA.DRA_CD82   HLA.…             0.883                  0.978              0.0343                 0.0343                0                1.00              1.00                   1                    0.983 M        
##  2 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DRA LAG3     HLA.DRA_LAG3   HLA.…             0.883                  0.978              0.0952                 0.0952                0                1.00              1.00                   1                    0.983 M        
##  3 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. SELPLG  SELL     SELPLG_SELL    SELP…             0.836                  0.942              0.139                  0.139                 0                1.00              1.00                   1                    0.973 M        
##  4 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. CD47    SIRPG    CD47_SIRPG     CD47…             0.816                  0.876              0.0497                 0.0497                0                1.00              1.00                   1                    0.962 M        
##  5 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DR… CD82     HLA.DRB1_CD82  HLA.…             0.925                  0.996              0.0343                 0.0343                0                1.00              1.00                   0.857                0.955 M        
##  6 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.DR… LAG3     HLA.DRB1_LAG3  HLA.…             0.925                  0.996              0.0952                 0.0952                0                1.00              1.00                   0.857                0.955 M        
##  7 M-(S+A)/2 M     L_NK_CD56._CD16.       L_T_TIM3._CD38._HLADR. CCL3    CCR5     CCL3_CCR5      CCL3…             0.933                  0.968              0.0532                 0.0532                0                1.00              1.00                   0.857                0.952 M        
##  8 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. HLA.F   CD8A     HLA.F_CD8A     HLA.…             0.885                  0.986              0.0399                 0.0399                0                1.00              1.00                   0.857                0.949 M        
##  9 M-(S+A)/2 M     L_T_TIM3._CD38._HLADR. M_Monocyte_CD16        IFNG    IFNGR1   IFNG_IFNGR1    IFNG…             0.129                  0.129              0.872                  0.891                 1.00             1.00              1.00                   0.857                0.948 M        
## 10 M-(S+A)/2 M     M_Monocyte_CD16        L_T_TIM3._CD38._HLADR. TNFSF13 FAS      TNFSF13_FAS    TNFS…             0.898                  0.947              0.0798                 0.0798                0                1.00              0.999                  0.857                0.945 M        
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​fraction_expressing_ligand_receptor
```

# Add the new prioritization tables to the saved output of MultiNicheNet 


```r
path = "./"

multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    prioritization_tables_with_condition_specific_celltype_sender = prioritization_tables_with_condition_specific_celltype_sender, 
    prioritization_tables_with_condition_specific_celltype_receiver = prioritization_tables_with_condition_specific_celltype_receiver, 
    combined_prioritization_tables = combined_prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 

multinichenet_output = make_lite_output_condition_specific(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_with_condition_specific_celltypes.rds"))

}
```

# Visualizations of the output of the condition-specific cell types workflows

## Visualization of differential cell-cell interactions -  results of the cell-cell communication analysis for the condition-specific cell type as SENDER

We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.


```r
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables_with_condition_specific_celltype_sender, 50, rank_per_group = FALSE)
```

### Interpretable bubble plots

Now we will visualize per sample the scaled product of ligand and receptor expression. Samples that were left out of the DE analysis are indicated with a smaller dot (this helps to indicate the samples that did not contribute to the calculation of the logFC, and thus not contributed to the final prioritization)

We will first visualize the interactions specific for the Sibling group that are part of the top50 differential interactions as shown in the above circos plots


```r
group_oi = "S"
```


```r
prioritized_tbl_oi_S_50 = prioritized_tbl_oi_all %>% filter(group == group_oi) %>% inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(multinichenet_output$prioritization_tables_with_condition_specific_celltype_sender, prioritized_tbl_oi_S_50)
plot_oi
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-57-1.png)<!-- -->

We will now check the MIS-C specific interactions


```r
group_oi = "M"
```


```r
prioritized_tbl_oi_M_50 = prioritized_tbl_oi_all %>% filter(group == group_oi) %>% inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(multinichenet_output$prioritization_tables_with_condition_specific_celltype_sender, prioritized_tbl_oi_M_50)
plot_oi
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

This plot deviates from the regular plots: in case the sender cell type is the condition-specific cell type "L_T_TIM3._CD38._HLADR." we see gray dots in the A and S condition, indicating it is not expressed there. However, we stil see the scaled cell type specicifity and ligand activity values.

Note: interpreting scaled expression values does not make much sense if only one condition where cell type is present.

## Visualization of differential cell-cell interactions -  results of the cell-cell communication analysis for the condition-specific cell type as SENDER

We will look here at the top 50 predictions across all contrasts, senders, and receivers of interest.


```r
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables_with_condition_specific_celltype_receiver, 50, rank_per_group = FALSE)
```

### Interpretable bubble plots

Now we will visualize per sample the scaled product of ligand and receptor expression. Samples that were left out of the DE analysis are indicated with a smaller dot (this helps to indicate the samples that did not contribute to the calculation of the logFC, and thus not contributed to the final prioritization)

We will first visualize the interactions specific for the Sibling group that are part of the top50 differential interactions as shown in the above circos plots


```r
group_oi = "S"
```


```r
prioritized_tbl_oi_S_50 = prioritized_tbl_oi_all %>% filter(group == group_oi) %>% inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(multinichenet_output$prioritization_tables_with_condition_specific_celltype_receiver, prioritized_tbl_oi_S_50)
plot_oi
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-62-1.png)<!-- -->

We will now check the MIS-C specific interactions


```r
group_oi = "M"
```


```r
prioritized_tbl_oi_M_50 = prioritized_tbl_oi_all %>% filter(group == group_oi) %>% inner_join(lr_network_all)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(multinichenet_output$prioritization_tables_with_condition_specific_celltype_receiver, prioritized_tbl_oi_M_50)
plot_oi
```

![](condition_specific_celltype_MISC_files/figure-html/unnamed-chunk-64-1.png)<!-- -->
This plot deviates from the regular plots: in case the receiver cell type is the condition-specific cell type "L_T_TIM3._CD38._HLADR." we see gray dots in the A and S condition for ligand-receptor expression, and also for the ligand activity prediction. However, we stil see the scaled cell type specicifity values.



```r
print("done")
## [1] "done"
```

# Note

Be cautious and critical while interpreting these results. These workflows mainly serve to still extract potentially interesting biological signal from your data in case you really have condition-specific cell types.
