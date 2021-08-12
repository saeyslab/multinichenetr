library(Seurat)
library(tidyverse)

# add visium
seurat_obj = readRDS("C:/Users/rbrowaey/work/Research/NicheNet/current_projects/CRC_NicheNet/data/seurat_obj_lite_hnscc.rds")
usethis::use_data(seurat_obj,overwrite = T, compress = "bzip2")


usethis::use_package("Seurat")
usethis::use_package("dplyr")
usethis::use_package("ggplot2")
usethis::use_package("circlize")
usethis::use_package("patchwork")

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
usethis::use_package("scater")

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
