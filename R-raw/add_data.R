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

usethis::use_package("circlize")
usethis::use_package("stringr")
usethis::use_package("generics")
usethis::use_package("ComplexHeatmap")
usethis::use_package("grid")
usethis::use_package("Nebulosa")



