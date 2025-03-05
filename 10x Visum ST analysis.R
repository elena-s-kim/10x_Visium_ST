#R script for the 10x Visum ST analysis, ***Work In Progress***
#Refs: Satijalab, LLM, UPItt CRC workshops (contains RMD parts)

library(Seurat)
library(SeuratData)
library(ggplot2)
library(spatstat.explore)
library(tidyverse)
library(patchwork)
library(dplyr)
library(openxlsx)
library(tools) 

# Define the paths: for the ST there are 2 of them per sample
matrix_dir <- "/ix/<...>/SpaceRanger_OUT/PBC-PR_6840-6A/outs/filtered_feature_bc_matrix/"
image_dir <- "/ix/<...>/SpaceRanger_OUT/PBC-PR_6840-6A/outs/spatial/"

# Step 1: Load the expression matrix
counts <- Read10X(data.dir = matrix_dir)

# Step 2: Create the Seurat object
seurat_obj6840 <- CreateSeuratObject(counts = counts, assay = "Spatial")

# Step 3: Load the spatial image
image <- Read10X_Image(image.dir = image_dir)

# Step 4: Add the spatial image to the Seurat object
DefaultAssay(image) <- "Spatial"
seurat_obj6840[["image"]] <- image

# Step 5: Visualize spatial features to confirm the image is integrated correctly
SpatialFeaturePlot(seurat_obj6840, features = "nCount_Spatial")

# Optional: Add metadata
#seurat_obj$sample <- "SampleName"

# Print object summary
print(seurat_obj6840)

vln.plot <- VlnPlot(seurat_obj6840, features = 'nCount_Spatial', pt.size = 0.1) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(seurat_obj6840, features = 'nCount_Spatial') + theme(legend.position = "right")
print(wrap_plots(vln.plot, count.plot))


seurat_obj6840 <- NormalizeData(seurat_obj6849, assay = "Spatial", normalization.method = "LogNormalize")

DefaultAssay(seurat_obj6837) <- "Spatial"
seurat_obj6840 <- SCTransform(seurat_obj6840, assay = "Spatial", verbose = FALSE)#double check here since SCT includes normalization (omit NormalizeData)

#TROUBLESHOOTING#

print(seurat_obj6840)
head(colnames(seurat_obj6840))

# Print a summary of the Seurat object
print(seurat_obj6840)

# Check the number of cells (columns) and features (rows)
num_cells <- ncol(seurat_obj6840)
num_features <- nrow(seurat_obj6840)

cat("Number of Cells:", num_cells, "\n")
cat("Number of Features (Genes):", num_features, "\n")


# Violin plot of UMI counts per spot
VlnPlot(seurat_obj6840, features = "nCount_Spatial", pt.size = 0.1) + 
  ggtitle("Distribution of UMI Counts per Cell/Spot") + NoLegend()

# Histogram of counts
metadata <- seurat_obj6840@meta.data
hist(metadata$nCount_Spatial, breaks = 50, main = "Distribution of nCount_Spatial",
     xlab = "UMI Counts per Spot", col = "lightblue")

######

# perform clustering workflow
seurat_obj6840 <- FindVariableFeatures(seurat_obj6840)
objseurat_obj6840ect <- ScaleData(seurat_obj6840)
seurat_obj6840 <- RunPCA(seurat_obj6840, assay="sketch", reduction.name = "pca.sketch")
seurat_obj6840 <- FindNeighbors(seurat_obj6840, assay="sketch", reduction = "pca.sketch", dims = 1:50)
seurat_obj6840 <- FindClusters(seurat_obj6840, cluster.name="seurat_cluster.sketched", resolution = 3)
seurat_obj6840 <- RunUMAP(seurat_obj6840, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
```

Now we can project the cluster labels, and dimensional reductions (PCA and UMAP) that we learned from the 50,000 sketched cells - to the entire dataset, using the `ProjectData` function. 

In the resulting object, for all cells:

* cluster labels will be stored in `object$seurat_cluster.projected`
* Projected PCA embeddings will be stored in `object[["pca.008um"]]`
* Projected UMAP embeddings will be stored in `object[["umap.sketch"]]`

```{r project}
object <- ProjectData(
  object = object,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
```

We can visualize the clustering results for the sketched cells, as well as the projected clustering results for the full dataset:

```{r project.plots, fig.height = 8}
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2
```

Of course, we can now also visualize the unsupervised clusters based on their spatial location. Note that running `SpatialDimPlot(object, interactive = TRUE)`, also enables interactive visualization and exploration.

```{r dim.plot}
SpatialDimPlot(object, label=T, repel=T, label.size = 4)
```

When there are many different clusters (some of which are spatially restricted and others are mixed), plotting the spatial location of all clusters can be challenging to interpret. We find it helpful to plot the spatial location of different clusters individually. For example, we highlight the spatial localization of a few clusters below, which happen to correspond to different cortical layers:

```{r cluster.plot, fig.height = 7}
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents=c(1,20,26))
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T) + NoLegend()
p

```

We can also find and visualize the top gene expression markers for each cluster:

```{r heatmap, fig.height = 12}
# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[['Spatial']]), downsample=1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = 'Spatial', only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
# p

# write file to CSV
write.csv(top5, file="top5.marker.gene.csv")

```


## Session information

```{r SessionInfo, echo=TRUE}
sessionInfo()

```

###### #####
####
# clear the working environment
rm(list=ls())

#options(future.globals.maxSize= 1e12)

# set the working directory
setwd("/ihome/hosmanbeyoglu/ysk13/CRC_Wrkshp_24/10X_Visium/analysis_10Xvisium")

# set the samples and conditions
SAMPLENAME = c("PBC_PR_6840_1023") #c("GSM6963120",  "GSM6963116") # sample name, unique name
LIBRARY = c("benign") #c("healthy",  "APAP") # condition, corresponding to SAMPLENAME

data=list() # data list

### read in the data
for(sInd in 1:length(SAMPLENAME)){
  print(c(sInd, SAMPLENAME[sInd]))
  
  ## read in data with hdf5 file 
  # data_dir <- paste0("../SpaceRanger_OUT/",SAMPLENAME[sInd],"/outs/")
  # list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
  # img = Read10X_Image(paste0("../SpaceRanger_OUT/", SAMPLENAME[sInd], "/outs/spatial/"), image.name = "tissue_hires_image.png")
  # data[[sInd]] <- Load10X_Spatial(data.dir = data_dir, image = img)
  
  
  ### read in data, based on gene-cell expression matrix file
  # count matrix file
  matrix_dir = paste0("/ix/hosmanbeyoglu/DATA/hosmanbeyoglu_wam30/Osmanbeyoglu_Visium_Lung/SpaceRanger_OUT/PBC-PR_6840-6A/outs/filtered_feature_bc_matrix/")
  counts = Seurat::Read10X(data.dir = matrix_dir)  
  
  data[[sInd]] = Seurat::CreateSeuratObject(
    counts = counts , 
    assay = 'Spatial')
  
  # image file
  imgpath = paste0("/ix/hosmanbeyoglu/DATA/hosmanbeyoglu_wam30/Osmanbeyoglu_Visium_Lung/SpaceRanger_OUT/PBC-PR_6840-6A/outs/spatial")
  img = Seurat::Read10X_Image(image.dir = imgpath)  
  
  Seurat::DefaultAssay(object = img) <- 'Spatial'  
  
  img = img[colnames(x = data)]  
  data[[sInd]][['image']] = img  
  
  # add in the metadata
  data[[sInd]]$slice = sInd
  data[[sInd]]$region = SAMPLENAME[sInd]
  data[[sInd]]$sample = SAMPLENAME[sInd]
  data[[sInd]]$library = LIBRARY[sInd]
  
  print(data[[sInd]])
  
}
names(data)=SAMPLENAME



#QC/ NORM

for(sInd in 1:length(data)){
  print(c(sInd, SAMPLENAME[sInd]))
  

  data[[sInd]] <- NormalizeData(data[[sInd]], assay = "Spatial", normalization.method = "LogNormalize")#I added EK
  
  ## SCT transformation / normalization within each slide
  data[[sInd]] <- SCTransform(data[[sInd]], assay = "Spatial", verbose = FALSE)
  
}


## data visualization
#pdf(paste0("Vln.",sInd,".",LIBRARY[sInd],".",SAMPLENAME[sInd],".pdf"), width=10, height=5)
plot1 <- VlnPlot(data[[sInd]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data[[sInd]], features = "nCount_Spatial") + theme(legend.position = "right")
print(wrap_plots(plot1, plot2))
#dev.off()


##########
#sample *** 6837 ***
SpatialFeaturePlot(data[[sInd]], features = c("BAP1", "CD86"))

seurat_obj6837 <- NormalizeData(seurat_obj6837, assay = "Spatial", normalization.method = "LogNormalize")



DefaultAssay(seurat_obj6837) <- "Spatial"
seurat_obj6837 <- SCTransform(seurat_obj6837, assay = "Spatial", verbose = FALSE)


#data_dir <- paste0("/ix/hosmanbeyoglu/DATA/hosmanbeyoglu_wam30/Osmanbeyoglu_Visium_05_2024/SpaceRanger_OUT/PBC_PR_6836/outs/")
#list.files(data_dir) # Should show filtered_feature_bc_matrix.h5

# Check alignment
print(all(rownames(seurat_obj6837@meta.data) == colnames(seurat_obj6837)))

# Fix alignment if necessary
rownames(seurat_obj6837@meta.data) <- colnames(seurat_obj6837)
######
DefaultAssay(seurat_obj6837) <- "Spatial"

# Normalize only the gene expression data
seurat_obj6837 <- NormalizeData(seurat_obj6837, assay = "Spatial", normalization.method = "LogNormalize", layer = "counts.Gene Expression")

# Perform SCTransform on the gene expression data
seurat_obj6837 <- SCTransform(seurat_obj6837, assay = "Spatial", verbose = FALSE, layer = "counts.Gene Expression")





SpatialFeaturePlot(seurat_obj6837, features = c("BAP1", "LAG3") + theme(legend.text = element_text(size = 0), legend.title = element_text(size = 20), legend.key.size = unit(1, "cm")))

SpatialFeaturePlot(seurat_obj6837, features = c("BAP1", "LAG3")) +
  theme(legend.text = element_text(size = 10), # Adjust size as needed
        legend.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"))


SpatialFeaturePlot(seurat_obj6837, features = c("CD4", "CD8A"))
SpatialFeaturePlot(seurat_obj6837, features = c("KRT19", "KRT8"))
SpatialFeaturePlot(seurat_obj6837, features = c("CD68", "CD30", "CD20"))#no CD20 found




SpatialFeaturePlot(seurat_obj6837, features = c("BAP1", "LAG3")) +
  theme(legend.text = element_text(size = 10), # Adjust size as needed
        legend.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"))
#### Modifying the image ####
SpatialFeaturePlot(seurat_obj6837, features = c("BAP1", "LAG3"), 
                   pt.size.factor = 3, # Increase dot size
                   alpha = c(0.2, 0.8)) + # Make background pale
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        strip.text = element_text(size = 15, face = "bold")) # Set facet label font size
######

# Store the plot
spatial_plot <- SpatialFeaturePlot(seurat_obj6837, features = c("BAP1", "LAG3")) +
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 20),
        legend.key.size = unit(1, "cm"))

# Print the plot
print(spatial_plot)



dir.create("../output/images", recursive = TRUE, showWarnings = FALSE)

print(
  SpatialFeaturePlot(seurat_obj, features = c("BAP1", "LAG3")) +
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 20),
          legend.key.size = unit(1, "cm"))
)


seurat_obj6837 <- RunPCA(seurat_obj6837, assay = "SCT", verbose = FALSE)
seurat_obj6837 <- FindNeighbors(seurat_obj6837, reduction = "pca", dims = 1:30)
print(seurat_obj6837@graphs)
seurat_obj6837 <- FindClusters(seurat_obj6837, graph.name = "SCT_snn", verbose = FALSE)
seurat_obj6837 <- RunUMAP(seurat_obj6837, reduction = "pca", dims = 1:30)

p1 <- DimPlot(seurat_obj6837, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat_obj6837, label = TRUE, label.size = 3)
p1 + p2
combined_plot <- p1 + p2

# Display the combined plot
print(combined_plot)
# Save to a PDF
ggsave("combined_plot.pdf", plot = combined_plot, width = 12, height = 6)
# Save to a PNG
ggsave("combined_plot.png", plot = combined_plot, width = 12, height = 6, dpi = 300)
getwd()

SpatialDimPlot<-SpatialDimPlot(seurat_obj6837, cells.highlight = CellsByIdentities(object = seurat_obj6837, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)), facet.highlight = TRUE, ncol = 3)
#ggsave("SpatialDimPlot.pdf", plot = SpatialDimPlot, width = 12, height = 6)
print(SpatialDimPlot)


SpatialDimPlot <- SpatialDimPlot(
  seurat_obj6837,
  cells.highlight = CellsByIdentities(object = seurat_obj6837, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
  facet.highlight = TRUE,
  ncol = 3,
  pt.size.factor = 5
) +
  scale_alpha(range = c(0.2, 0.9))  # Adjust to make the background dimmer
print(SpatialDimPlot)
#ggsave("SpatialDimPlot.pdf", plot = SpatialDimPlot, width = 12, height = 6)

data("seurat_obj6837")
VlnPlot(object = seurat_obj6837, features = 'BAP1', split.by = clusters)
print(VlnPlot)


plot1 <- VlnPlot(seurat_obj6837, features = "idents", pt.size = 0.1) + NoLegend()



++
  # SpatialDimPlot
  SpatialDimPlot <- SpatialDimPlot(
    seurat_obj6837,
    cells.highlight = CellsByIdentities(object = seurat_obj6837, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
    facet.highlight = TRUE,
    ncol = 3,
    pt.size.factor = 5
  ) +
  scale_alpha(range = c(0.2, 0.9))  # Adjust background dimness
print(SpatialDimPlot)

# Save plot
ggsave("SpatialDimPlot.pdf", plot = SpatialDimPlot, width = 12, height = 6)

# Violin Plot for BAP1 Expression (by cluster identities)
plot2 <- VlnPlot(object = seurat_obj6837, features = 'BAP1', split.by = "ident")
print(plot2)

# Violin Plot for cluster identities (metadata: Idents)
plot3 <- VlnPlot(seurat_obj6837, features = "BAP1", pt.size = 0.1) + NoLegend()
print(plot3)

++
  cluster5_markers <- FindMarkers(
    seurat_obj6837, 
    ident.1 = 5, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )

# Top 10 Upregulated and Downregulated Genes
cluster5_markers <- FindMarkers(
  seurat_obj6837, 
  ident.1 = 5, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)
cluster5_up <- head(cluster5_markers[order(-cluster5_markers$avg_log2FC), ], 10)
cluster5_down <- head(cluster5_markers[order(cluster5_markers$avg_log2FC), ], 10)
# Cluster 5
print("Cluster 5 - Top 10 Upregulated Genes")
print(cluster5_up)

print("Cluster 5 - Top 10 Downregulated Genes")
print(cluster5_down)

FeaturePlot(seurat_obj6837, features = c("CD209", "TRAC", "MS4A7", "FOLR2", "KLRK1"))

plot3 <- VlnPlot(seurat_obj6837, features = c("CD209", "TRAC", "MS4A7", "FOLR2", "KLRK1"), pt.size = 0.1) + NoLegend()
print(plot3)

install.packages("patchwork")

library(patchwork)  # Ensure the patchwork package is loaded

# Generate violin plots for multiple features
plot3 <- VlnPlot(
  seurat_obj6837, 
  features = c("CD209", "TRAC", "MS4A7", "FOLR2", "KLRK1"), 
  pt.size = 0.1
) & NoLegend()  # Apply NoLegend to all subplots

# Print the combined plot
print(plot3)

