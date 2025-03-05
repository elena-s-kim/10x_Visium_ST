#SpatialnferCNV, followed https://aerickso.github.io/SpatialInferCNV/
#https://github.com/aerickso/SpatialInferCNV

install.packages("devtools")
library(devtools)
install.packages("usethis",force = TRUE)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)
library(tidyverse)

download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5", "./V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5", mode = "wb")
Breast_ENSBMLID_Counts <- ImportCountData("Breast10X", "./V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")

head(Breast_ENSBMLID_Counts)
Userguide_10xBreast_Histology <- ImportHistologicalAnnotations("Breast10X", "/Users/elena_s_kim/Documents/SpatialInferCNV/10xBreast_UserguideHistologyAnnotations.csv")
head(Userguide_10xBreast_Histology)

Userguide_10xBreast_Joined_Counts <- MergingCountAndAnnotationData("Breast10X",Userguide_10xBreast_Histology, Breast_ENSBMLID_Counts)

Userguide_10xBreast_Joined_Counts <- Userguide_10xBreast_Joined_Counts %>% column_to_rownames(var = "Genes")

rm(Breast_ENSBMLID_Counts)

head(Userguide_10xBreast_Joined_Counts)

dim(Userguide_10xBreast_Histology)

Userguide_10xBreast_Joined_Counts <- Userguide_10xBreast_Joined_Counts %>% select(-Breast10X_AATAACGTCGCGCCCA.1)

FinalAnnotationsForExport <- FinalAnnotations(Userguide_10xBreast_Histology, Userguide_10xBreast_Joined_Counts)

dim(FinalAnnotationsForExport)

write.table(Userguide_10xBreast_Joined_Counts, "Userguide_10xBreast_Joined_Counts.tsv", sep = "\t")

write.table(FinalAnnotationsForExport, "Userguide_10xBreastFinalAnnotationsForExport.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

getwd()

#download.file("https://raw.githubusercontent.com/aerickso/SpatialInferCNV/main/FigureScripts/siCNV_GeneOrderFile.tsv?token=GHSAT0AAAAAABRUPXQSKNJYB6FRBX2GLUGMYR77RCQ", "./siCNV_GeneOrderFile.tsv", mode = "wb")
#downloaded from https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/siCNV_GeneOrderFile.tsv

Userguide_10xBreastCancer_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Userguide_10xBreast_Joined_Counts.tsv", 
                                                                   gene_order_file="/Users/elena_s_kim/Documents/SpatialInferCNV/siCNV_GeneOrderFile.tsv",#./siCNV_GeneOrderFile.tsv",
                                                                   annotations_file="./Userguide_10xBreastFinalAnnotationsForExport.tsv",
                                                                   delim="\t",
                                                                   ref_group_names=NULL, 
                                                                   chr_exclude = c("chrM")) 

Userguide_10xBreastCancer_infCNV = infercnv::run(Userguide_10xBreastCancer_infCNV, 
                                                 cutoff=0.1, #(see infercnv::run documentation)
                                                 out_dir="./InferCNVrun_outputs", 
                                                 cluster_by_groups=FALSE, #unsupervised analysis
                                                 HMM = FALSE, 
                                                 denoise=TRUE) #denoising applies noise reduction for the plot 


library("devtools")
devtools::install_github("broadinstitute/infercnvApp")

infercnvApp::infercnvApp()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

library(infercnv)
library(Seurat)
getwd()
setwd("/Users/elena_s_kim/Documents/SpatialInferCNV")
output_dir_full = "output_dir_full"

#devtools::install_github("broadinstitute/infercnvNGCHM") #intrvtv vslshn
#devtools::install_github("bmbroom/tsvio")
#devtools::install_github("bmbroom/NGCHMR", ref="stable")

#Contitue https://aerickso.github.io/SpatialInferCNV/
library(ape)# R version 4.1.3 
library(phylogram)# R version 4.1.3 

BreastCancer10x_for_clustering <- read.dendrogram(file = "/Users/elena_s_kim/Documents/SpatialInferCNV/InferCNVrun_outputs/infercnv.observations_dendrogram.txt")#./UserGuideFiles/infercnv.21_denoised.observations_dendrogram.txt")
#Convert to a phylo object using as.phylo()
BreastCancer10x_for_clustering_phylo <- as.phylo(BreastCancer10x_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(BreastCancer10x_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("BreastCancer10x_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(BreastCancer10x_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:BreastCancer10x_for_clustering_phylo$Nnode,node=1:BreastCancer10x_for_clustering_phylo$Nnode+Ntip(BreastCancer10x_for_clustering_phylo))
dev.off()
