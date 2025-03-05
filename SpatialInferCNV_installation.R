#R  SpatialInferCNV Installation and Breast Cancer example
library(tidyverse)
library(infercnv)
library(data.table)
library(SpatialInferCNV)
install.packages("BiocManager")
install.packages("SpatialInferCNV")
BiocManager::install("tidyverse")
BiocManager::install("SpatialInferCNV")
Yes

install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
install.packages("tidyverse")
install.packages("Seurat")
install.packages("phylogram")
install.packages("ape")
install.packages("hdf5r")

install.packages("devtools")
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)

download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5", "./V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5", mode = "wb")

Breast_ENSBMLID_Counts <- ImportCountData("Breast10X", "./V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
install.packages("remotes")
remotes::install_github("aerickso/SpatialInferCNV")

remotes::install_version("tidyverse", version = "1.3.1")
install.packages("tidyverse", version = "1.3.1")
getwd()
list.files("/Users/elena_s_kim")
setwd("/Users/elena_s_kim/SpatialInferCNV")

Breast_ENSBMLID_Counts <- ImportCountData("Breast10X", "./V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
install.packages("remotes")
remotes::install_version("tidyverse", version = "1.3.1")
install.packages("tidyverse", version = "1.3.1")
No
getwd()
#/Users/elena_s_kim/SpatialInferCNV/10xBreast_UserguidHistologyAnnotations.csv
#Userguide_10xBreast_Histology <- ImportHistologicalAnnotations("Breast10X", "./UserGuideFiles/10xBreast_UserguideHistologyAnnotations.csv")
Userguide_10xBreast_Histology <- ImportHistologicalAnnotations("Breast10X", "/Users/elena_s_kim/SpatialInferCNV/10xBreast_UserguidHistologyAnnotations.csv")
head(Userguide_10xBreast_Histology)
