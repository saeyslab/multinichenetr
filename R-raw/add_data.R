library(tidyverse)

# in the past, we used seurat objects as input
# now this was changed to single-cell experiment
# seurat_obj = readRDS("C:/Users/rbrowaey/work/Research/NicheNet/current_projects/CRC_NicheNet/data/seurat_obj_lite_hnscc.rds")
# usethis::use_data(seurat_obj,overwrite = T, compress = "bzip2")
# convert seurat to SCE object
# sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce = readRDS("C:/Users/rbrowaey/work/Research/NicheNet/sce_hnscc.rds")
set.seed(1919)
sample_id = "tumor"
extra_metadata = SummarizedExperiment::colData(sce)  %>% tibble::as_tibble()  %>% dplyr::select(all_of(sample_id)) %>% dplyr::distinct() %>% mutate(batch = sample(c("A","B"),nrow(.), replace = TRUE))
new_metadata = SummarizedExperiment::colData(sce)  %>% tibble::as_tibble()  %>% inner_join(extra_metadata)
new_metadata = new_metadata %>% data.frame()
rownames(new_metadata) = new_metadata$cell

sce = SingleCellExperiment::SingleCellExperiment(list(counts=SingleCellExperiment::counts(sce), logcounts=SingleCellExperiment::logcounts(sce)),
                                                 reducedDims = SingleCellExperiment::reducedDims(sce), 
                                                 colData=new_metadata,
                                                 rowData=SingleCellExperiment::rowData(sce),
                                                 metadata=sce@metadata
)
usethis::use_data(sce,overwrite = T, compress = "bzip2")


# usethis::use_package("Seurat") # not used anymore
usethis::use_package("dplyr")
usethis::use_package("ggplot2")
usethis::use_package("circlize")
usethis::use_package("patchwork")
usethis::use_package("SingleCellExperiment") # because of the data input dependency!
usethis::use_package("ggbeeswarm")
usethis::use_package("Scater")
usethis::use_package("Muscat")
usethis::use_package("tibble")
usethis::use_package("tidyr")
usethis::use_package("purrr")
usethis::use_package("ComplexHeatmap")
usethis::use_package("parallel")

usethis::use_package("circlize")
usethis::use_package("stringr")
usethis::use_package("generics")
usethis::use_package("ComplexHeatmap")
usethis::use_package("grid")
usethis::use_package("Nebulosa")


usethis::use_package("muscat")
usethis::use_package("limma")
usethis::use_package("SummarizedExperiment")
usethis::use_package("S4Vectors")
usethis::use_package("magrittr")
usethis::use_package("nichenetr")
# usethis::use_package("scater")

usethis::use_package("RColorBrewer")
usethis::use_package("ggpubr")
usethis::use_package("edgeR")
usethis::use_package("sva")



usethis::use_package("knitr", type = "Suggests")
usethis::use_package("testthat", type = "Suggests")
usethis::use_package("covr", type = "Suggests")
usethis::use_package("rmarkdown", type = "Suggests")
usethis::use_package("tidyverse", type = "Suggests")
usethis::use_package("locfdr")

usethis::use_package("UpSetR", type = "Suggests")
usethis::use_package("scran")
