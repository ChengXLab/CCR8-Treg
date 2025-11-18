setwd("D:/R/monocle")
library(monocle3)
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)

Treg.intergrated <- readRDS("D:/R/monocle/Treg.intergrated.rds")
DimPlot(Treg.intergrated, reduction = "umap")

#RUN monocle3
cds2 <- as.cell_data_set(x = Treg.intergrated)

#mapping cluster
cds2@clusters$UMAP$clusters <- Idents(Treg.intergrated)[rownames(colData(cds2))]
cds2@clusters$UMAP$orig.ident <- Idents(Treg.intergrated)[rownames(colData(cds2))]
cds2@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(cds2)))), levels = 1)
names(cds2@clusters$UMAP$partitions) <- rownames(colData(cds2))
cds2 <- estimate_size_factors(cds2)

## Add gene names into CDS
cds2@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds2)
rownames(cds2@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@reducedDims$UMAP) <- NULL
colnames(cds2@int_colData@listData$reducedDims@listData$UMAP) <- NULL

cds2 <- learn_graph(cds2, use_partition = F)
cds2 <- order_cells(cds2)

plot_cells(cds2,color_cells_by = "seurat_clusters",
           cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=.1,
           trajectory_graph_segment_size = 1.5)+ 
  theme(legend.position = "none")+
  scale_color_manual(values = c("#71ACD8","#C14F58", "#E9A66E", "#F4DA90","#60C2A6","#DC7656","#DCE49F"))& NoAxes()

ggsave("monocle3.pdf", width = 3.5 , height = 3.5, dpi=600)

