library(nichenetr)
library(tidyr)
library(dplyr)
library(Seurat)
library(testthat)

test_check("spatialnichenetr", filter = "visium_scrnaseq_ccc")
