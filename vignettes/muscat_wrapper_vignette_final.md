Muscat wrapper: HNSCC application
================
Robin Browaeys
2021-07-08

<!-- github markdown built using 
rmarkdown::render("vignettes/muscat_wrapper_vignette_final.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a muscat differential
state (DS) analysis. A DS analysis can be performed if you have
multi-sample, multi-group single-cell data. For each cell type of
interest, muscat will compare the sample-wise expression of all genes
between groups of interest. Therefore, the absolute minimum of meta data
you need to have, are following columns indicating for each cell: the
**group**, **sample** and **cell type**.

As example expression data, we will use data from Puram et al. to
explore intercellular communication in the tumor microenvironment in
head and neck squamous cell carcinoma (HNSCC) \[See
@puram\_single-cell\_2017\]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4675430.svg)](https://doi.org/10.5281/zenodo.4675430).
The groups we have here are tumors scoring high for a partial
epithelial-mesenschymal transition (p-EMT) program vs low-scoring
tumors.

The different steps of the MultiNicheNet analysis are the following:

-   0.  Preparation of the analysis: load packages, read in the
        single-cell expression data, and define the main settings of the
        muscat analysis

-   1.  Check cell type abundance for the cell types of interest

-   2.  Perform genome-wide differential expression (DS) analysis

-   3.  Downstream analysis of the DS output, including visualization

In this vignette, we will demonstrate all these steps in detail.

# Step 0: Preparation of the analysis: load packages, read in the single-cell expression data, and define the main settings of the muscat analysis

The current implementation of the muscat wrapper starts from a
SingleCellExperiment object, therefore we will need to load the
SingleCellExperiment library. If you start from a Seurat object, you can
convert it easily to a SingleCellExperiment via
`sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")`.

``` r
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
```

In this case study, we want to study differences in expression between
pEMT-high and pEMT-low tumors. The meta data columns that indicate the
pEMT status of tumors are ‘pEMT’ and ‘pEMT\_fine’, cell type is
indicated in the ‘celltype’ column, and the sample is indicated by the
‘tumor’ column.

**User adaptation required**

``` r
sce = readRDS(url("https://zenodo.org/record/5196144/files/sce_hnscc.rds"))
scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "celltype")
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "tumor")
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "pEMT")
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "pEMT_fine")
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

Now we will define in which metadata columns we can find the **group**,
**sample** and **cell type** IDs

For the group\_id, we now choose for the ‘pEMT’ column instead of
‘pEMT\_fine’, which we will select in a subsequent analysis.

**User adaptation required**

``` r
sample_id = "tumor"
group_id = "pEMT"
celltype_id = "celltype"
```

Now we will go to the first real step of the muscat analysis

# Step 1: Check cell type abundance

## Step 1a: Number of cells per celltype-sample combination

We will now check the number of cells per cell type condition
combination, and the number of patients per condition. This is important
because muscat performs pseudobulking to infer group differences at the
sample level for each cell type. This means that we will group the
information of all cells of a cell type in a sample together to get 1
sample-celltype estimate. The more cells we have, the more accurate this
aggregated expression measure will be.

**User adaptation required**

``` r
table(SummarizedExperiment::colData(sce)$celltype, SummarizedExperiment::colData(sce)$tumor) # cell types vs samples
##                
##                 HN16 HN17 HN18 HN20 HN22 HN25 HN26 HN28 HN5 HN6
##   CAF             47   37   36    3    9   73   19  157  37  82
##   Endothelial     43   17   18    1    1    1    0   14  11  52
##   Malignant       82  353  263  331  123  153   61   49  70 157
##   Myeloid         15    2    7    0    1    8    1    1  58   6
##   myofibroblast   84    6   14   10   45   88   45  140   5   6
##   T.cell         300   61  207    0    0   93    3    0  28   0
table(SummarizedExperiment::colData(sce)$celltype, SummarizedExperiment::colData(sce)$pEMT) # cell types vs conditions
##                
##                 High  Low
##   CAF            396  104
##   Endothelial    105   53
##   Malignant     1093  549
##   Myeloid         92    7
##   myofibroblast  382   61
##   T.cell         689    3
table(SummarizedExperiment::colData(sce)$tumor, SummarizedExperiment::colData(sce)$pEMT) # samples vs conditions
##       
##        High Low
##   HN16  571   0
##   HN17  476   0
##   HN18  545   0
##   HN20    0 345
##   HN22  179   0
##   HN25  416   0
##   HN26    0 129
##   HN28  361   0
##   HN5   209   0
##   HN6     0 303
```

As you can see in the upper table, some Celltype-Sample combinations
have 0 cells. It is possible that during DE analysis, some cell types
will be removed from the analysis if there is not enough information to
do a DE analysis. (More info later)

We can define the minimum number of cells we require per celltype-sample
combination. It is recommened to have at least 10 (and preferably more)
cells in each sample-celltype combination. Therefore we will set the
`min_cells` parameter in the analysis to 10. Celltype-sample
combinations with less cells will not be considered during the muscat DS
analysis.

**User adaptation possible**

``` r
min_cells = 10
```

To visually see which celltype-sample combinations won’t be considered,
you can run the following code:

``` r
metadata_abundance = SummarizedExperiment::colData(sce) %>% tibble::as_tibble() %>% .[,c(sample_id, group_id, celltype_id)]
colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
  
abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% dplyr::distinct(sample_id , group_id ), by = "sample_id")
abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
  
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
abund_plot
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
Celltype-sample combinations that won’t be considered are indicated in
red (because they have less cells than the `min_cells` threshold
indicated by the red dashed line)

If too many celltype-sample combinations don’t pass this threshold, we
recommend to define your cell types in a more general way it this would
still be possible and make sense biologically (–&gt; use one level
higher of the cell type ontology hierarchy; eg TH17 CD4T cells –&gt;
CD4T cells \| but not myeloid + T.cell together).

We can see here that quite many sample-celltype combinations are left
out. For Endothelial, Myeloid, and T cells, we don’t even have two or
more samples that have enough cells of those cell types. When we don’t
have two or more samples per group left, we cannot do a group comparison
(we need at least 2 replicates per group for a statistical analysis).
Therefore, those cell types will be removed before the DE analysis.

As stated before when seeing this, we would recommend to use a
higher-level cell type annotation if possible. But the annotation here
is already high-level, and grouping Endothelial cells, T cells and
Myeloid cells eg would not make sense biologically. That we won’t be
able to include these cell types in our analysis is a limitation of the
MultiNicheNet approach compared to classic cell-level-based approaches
(like Seurat::FindMarkers). On the contrary, those cell-level-based
approaches don’t reveal the lack of cells in many samples, and might
lead to biased results.

## Step 1b: Differential cell type abundance between the groups of interest

In another visualization, we can compare the cell type abundances
between the groups of interest. This can be interesting because too
strong abundance differences might have an effect on the DS analysis.
Downstream results of these cell types should then be considered with
some caution.

``` r
abund_plot_group = abundance_data %>% ggplot(aes(group_id, n, group = group_id, color = group_id)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(aes(alpha = keep), width = 0.15, height = 0.05) + scale_alpha_manual(values = c(1,0.30)) + facet_wrap( ~ celltype_id, scales = "free") + theme_bw() + 
    scale_color_discrete("tomato","steelblue2") + geom_hline(yintercept = min_cells, color = "red", linetype  = "longdash") + ggtitle("Cell type abundances per group") + ylab("# cells per sample-celltype combination") + xlab("Group")
abund_plot_group
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
Differential abundance looks quite OK for the cell types kept for the DE
analysis (i.e. CAF, Malignant and myofibroblast)

### Conclusion of this step:

**Important**: Based on the cell type abundance diagnostics, we
recommend users to change their analysis settings if required, before
proceeding with the rest of the analysis.

# Step 2: Perform genome-wide differential expression analysis

Now we will go over to the multi-group, multi-sample differential
expression (DE) analysis (also called ‘differential state’ analysis by
the developers of Muscat).

### Define the contrasts and covariates of interest for the DE analysis.

Here, we want to compare the p-EMT-high vs the p-EMT-low group and find
cell-cell communication events that are higher in high than low pEMT. We
don’t have other covariates to correct for in this dataset. If you would
have covariates you can correct for, we strongly recommend doing this,
since this is one of the main unique possibilities of the MultiNicheNet
approach.

#### about covariates:

Note that it is only possible to add a covariate if the different
covariate categories are present in all your groups of interest as
defined in the contrasts. Eg adding the covariate ‘sex’ is possible if
both group 1 and 2 contain male and female samples. It would not be
possible if group 2 does not contain male samples for example.

If you have paired data, meaning that you have samples in group 1 (eg
steady-state) and group 2 (eg treatment) coming from the same patient,
we strongly recommend exploiting this benefit in your experimental
design by using the patient id as covariate.

#### about contrasts and how to set them:

Note the format to indicate the contrasts! (This formatting should be
adhered to very strictly, and white spaces are not allowed)

**User adaptation required**

``` r
covariates = NA
contrasts_oi = c("'High-Low','Low-High'")
contrast_tbl = tibble(contrast = 
                        c("High-Low","Low-High"), 
                      group = c("High","Low")) # not necessary to be the same as the groups in your data, but recommended for later visualization
```

### Perform the DE analysis for each cell type.

``` r
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)
## [1] "excluded cell types are:"
## [1] "Endothelial" "Myeloid"     "T.cell"     
## [1] "These celltypes are not considered in the analysis. After removing samples that contain less cells than the required minimal, some groups don't have 2 or more samples anymore. As a result the analysis cannot be run. To solve this: decrease the number of min_cells or change your group_id and pool all samples that belong to groups that are not of interest! "
```

### Check DE results

The `DE_info` object contains both the default output as given by the
`muscat::pbDS()` function, and a cleaner output table (through
`muscat::resDS`).

Table with logFC and p-values for each gene-celltype-contrast:

``` r
DE_info$celltype_de$de_output_tidy
## # A tibble: 67,972 x 9
##    gene         cluster_id   logFC logCPM       F  p_val p_adj.loc p_adj contrast
##    <chr>        <chr>        <dbl>  <dbl>   <dbl>  <dbl>     <dbl> <dbl> <chr>   
##  1 RPS11        CAF        -0.0218   9.3  0.00697 0.935          1 0.993 High-Low
##  2 ELMO2        CAF         0.639    6.52 2.27    0.156          1 0.704 High-Low
##  3 CREB3L1      CAF         0.264    5.5  0.0617  0.808          1 0.974 High-Low
##  4 PNMA1        CAF        -1.01     6.51 4.01    0.0673         1 0.549 High-Low
##  5 MMP2         CAF         0.115    9.2  0.213   0.653          1 0.945 High-Low
##  6 TMEM216      CAF        -1.28     6.02 5.94    0.0304         1 0.424 High-Low
##  7 TRAF3IP2-AS1 CAF        -0.438    5.54 0.654   0.434          1 0.889 High-Low
##  8 ZHX3         CAF        -0.113    5.5  0.0391  0.846          1 0.979 High-Low
##  9 ERCC5        CAF        -0.066    5.85 0.0139  0.908          1 0.989 High-Low
## 10 APBB2        CAF        -0.211    6.46 0.232   0.638          1 0.942 High-Low
## # ... with 67,962 more rows
```

We can also show the distribution of the p-values:

``` r
DE_info$hist_pvals
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
(Note: this p-value histograms are the same for High-Low and Low-High
because we only have two groups and compare them to each other - a DE
gene in one comparison will then also be DE in the other comparison,
with just a reversed sign of the logFC)

In order to trust the p-values, the p-value distributions should be
uniform distributions, with a peak allowed between 0 and 0.05 if there
would be a clear biological effect in the data. This clear effect
(=clear DE) seems to be present here in the Malignant cell type
populations, although the histogram is not very uniformly distributed
for p-values between 0.05 and 0.25. This might point to issues in the DE
model definition (eg we did not add all important covariates,
substructure present,…)

Because there might be some issues, and we anticipate this could be
present in other datasets, we will now use the empiricall null
procedure. This is a procedure that will define empirical p-values based
on the observed distribution of the test statistic (here: logFC) and not
based on the theoretical distribution. We only recommend this if the
p-value distributions point to possible issues. (\[ADD REFERENCE\])
NOTE: IF USING THIS: COMMUNICATE WITH JEROEN GILLIS

### Empirical Null procedure

**User adaptation recommended**

``` r
empirical_pval = TRUE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
} 
```

Table with logFC and p-values for each gene-celltype-contrast:

``` r
DE_info_emp$de_output_tidy_emp
## # A tibble: 67,972 x 11
##    gene         cluster_id   logFC logCPM       F  p_val p_adj.loc p_adj contrast  p_emp p_adj_emp
##    <chr>        <chr>        <dbl>  <dbl>   <dbl>  <dbl>     <dbl> <dbl> <chr>     <dbl>     <dbl>
##  1 RPS11        CAF        -0.0218   9.3  0.00697 0.935          1 0.993 High-Low 0.863      0.995
##  2 ELMO2        CAF         0.639    6.52 2.27    0.156          1 0.704 High-Low 0.128      0.994
##  3 CREB3L1      CAF         0.264    5.5  0.0617  0.808          1 0.974 High-Low 0.846      0.995
##  4 PNMA1        CAF        -1.01     6.51 4.01    0.0673         1 0.549 High-Low 0.0317     0.994
##  5 MMP2         CAF         0.115    9.2  0.213   0.653          1 0.945 High-Low 0.668      0.994
##  6 TMEM216      CAF        -1.28     6.02 5.94    0.0304         1 0.424 High-Low 0.0115     0.994
##  7 TRAF3IP2-AS1 CAF        -0.438    5.54 0.654   0.434          1 0.889 High-Low 0.335      0.994
##  8 ZHX3         CAF        -0.113    5.5  0.0391  0.846          1 0.979 High-Low 0.765      0.994
##  9 ERCC5        CAF        -0.066    5.85 0.0139  0.908          1 0.989 High-Low 0.833      0.994
## 10 APBB2        CAF        -0.211    6.46 0.232   0.638          1 0.942 High-Low 0.541      0.994
## # ... with 67,962 more rows
```

The following plot shows the distribution of these corrected, empirical
p-values:

``` r
DE_info_emp$hist_pvals_emp
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
The following plots show how well the correction worked. The green
fitted curve should fit well with the histogram. If not, this might
point to some issues in the DE model definition.

``` r
DE_info_emp$z_distr_plots_emp_pval
## $`CAF.High-Low`
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

    ## 
    ## $`CAF.Low-High`

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

    ## 
    ## $`Malignant.High-Low`

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

    ## 
    ## $`Malignant.Low-High`

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

    ## 
    ## $`myofibroblast.High-Low`

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

    ## 
    ## $`myofibroblast.Low-High`

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->

In general, these plots looks fine, except for the Malignant cells. As
discussed in the previous plots: there might be an issue here. One
possible explanation might be that there is additional substructure in
the data. This could make sense because the pEMT high group consists
both of pEMT very-high and pEMT high samples according to the finer
subdivision/metadata column `pEMT_fine`. Another reason for the possible
substructure: malignant cells of each tumor are very different because
of genetic aberrations.

If we compare the empirical p-values to the original ones before for the
malignant cells, we see that through the empirical null procedure, we
lost some DE genes in malignant cells. But it is likely that the DE
genes that er kept are the most bona fide ones.

As additional check, we will look for the concordance between p-values
rankings of the original and empirical DE analysis (via ranking-line and
upset plots):

``` r
comparison_plots = DE_info$celltype_de$de_output_tidy$cluster_id %>% unique() %>% lapply(function(celltype_oi, adjusted = FALSE){
  if(adjusted == TRUE){
      de_genes_normal = DE_info$celltype_de$de_output_tidy %>% filter(cluster_id == celltype_oi) %>% filter(p_adj.glb <= 0.05) %>% pull(gene) %>% unique()
      de_genes_emp = DE_info_emp$de_output_tidy_emp %>% filter(cluster_id == celltype_oi) %>% filter(p_adj_emp <= 0.05) %>% pull(gene) %>% unique()

  } else {
      de_genes_normal = DE_info$celltype_de$de_output_tidy %>% filter(cluster_id == celltype_oi) %>% filter(p_val <= 0.05) %>% pull(gene) %>% unique()
      de_genes_emp = DE_info_emp$de_output_tidy_emp %>% filter(cluster_id == celltype_oi) %>% filter(p_emp <= 0.05) %>% pull(gene) %>% unique()

  }

  upset_df = tibble(gene = union(de_genes_normal, de_genes_emp), normal = as.double(gene %in% de_genes_normal), empirical = as.double(gene %in% de_genes_emp)) %>% data.frame() %>% magrittr::set_rownames(.$gene) %>% dplyr::select(-gene)
  colnames(upset_df) = paste(colnames(upset_df), celltype_oi, sep = "-")
  p_upset = UpSetR::upset(upset_df, sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on") 
  
  p_ranking = DE_info_emp$de_output_tidy_emp %>% filter(gene %in% union(de_genes_normal, de_genes_emp) & cluster_id == celltype_oi) %>% group_by(cluster_id, contrast) %>% mutate(normal = rank(p_val), empirical = rank(p_emp)) %>% filter(normal != empirical) %>% mutate(empirical_lower = empirical < normal) %>% tidyr::gather(rank_type, rank, normal:empirical) %>% dplyr::select(gene, rank_type, rank, empirical_lower)  %>% 
    ggplot(aes(rank_type, rank, group = gene, color = empirical_lower)) + geom_line(aes(group = gene)) + facet_grid(cluster_id ~ contrast) + theme_bw()
  
  return(list(p_upset, p_ranking))
  
}, adjusted = FALSE) 
comparison_plots
## [[1]]
## [[1]][[1]]
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## 
    ## [[1]][[2]]

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

    ## 
    ## 
    ## [[2]]
    ## [[2]][[1]]

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

    ## 
    ## [[2]][[2]]

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

    ## 
    ## 
    ## [[3]]
    ## [[3]][[1]]

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-5.png)<!-- -->

    ## 
    ## [[3]][[2]]

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-16-6.png)<!-- -->

### Conclusion of the diagnostic plots concerning the DE analysis

P-value histograms of both the normal and empirical p-values indicate
there might be some problems in the DE model definition, certainly for
malignant cells. If possible it might be a good idea to include more
covariates in the model, or use the `pEMT_fine` group definition instead
(which we will do later in this vignette).

# Step 3: Downstream analysis and visualization

``` r
receiver_oi = "Malignant"
group_oi = "High"

targets_oi = DE_info$celltype_de$de_output_tidy %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% filter(cluster_id == receiver_oi) %>% filter(p_adj <= 0.05) %>% arrange(p_adj) %>% pull(gene) %>% unique()
```

First, make a violin plot

``` r
target_oi = targets_oi[1]

make_target_violin_plot(sce_receiver = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id_receiver = celltype_id)
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Then a feature plot

``` r
target_oi = targets_oi[1]

make_target_feature_plot(sce_receiver = sce, target_oi = target_oi, group_oi = group_oi, group_id = group_id, celltype_id_receiver = celltype_id, receivers_oi = c("Malignant","myofibroblast","CAF")) 
```

![](muscat_wrapper_vignette_final_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Then, a dotplot THE FOLLOWING REQUIRES THE PSEUDOBULKED COUNTS
INFORMATION ETC

``` r
# p_dotplot = make_sample_target_plots(receiver_info = multinichenet_output$celltype_info, targets_oi, receiver_oi, output$grouping_tbl)
#p_dotplot + ggtitle(paste0("DE genes in ",group_oi, " in celltype ",receiver_oi))
```
