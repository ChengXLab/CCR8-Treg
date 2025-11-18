BiocManager::available()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

BiocManager::install(c('lme4', 'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3',force = TRUE)

remotes::install_github('satijalab/seurat-wrappers')

setwd("D:/4.ccr8/R/0902/monocle")
library(monocle3)
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)

Treg.intergrated <- readRDS("D:/4.ccr8/0719/Treg/Treg.intergrated.rds")
DimPlot(Treg.intergrated, reduction = "umap")
#Treg.intergrated <- Heart
#Treg.intergrated = Treg.intergrated[,Treg.intergrated@meta.data$orig.ident %in% c("LN")]

#Treg.intergrated <- scRNA

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

plot_cells(cds2, color_cells_by = "pseudotime",
           #show_trajectory_graph=F,
           cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0.5,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1.3)+
  theme(legend.position = "none") & NoAxes()

ggsave("monocle3_pse.pdf", width = 3.5 , height = 3.5, dpi=600)

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

plot_cells(cds2,color_cells_by = "orig.ident",
           cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=.1,
           trajectory_graph_segment_size = 1.5)+ 
  theme(legend.position = "none")+
  scale_color_manual(values = c( "#C14F58","#ffad73","#26b3ff","#DC7656","#DCE49F","#B77B13"))& NoAxes()

ggsave("monocle3_orig.ident.pdf", width = 3.5 , height = 3.5, dpi=600)


ciliated_cds_pr_test_res <- graph_test(cds2, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- ciliated_cds_pr_test_res %>% filter(q_value < 0.05) %>% arrange(-morans_I)

plot_cells(cds2, genes=c("Ccr8","Klrg1","Sell","Gzmb"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

plot_genes_in_pseudotime(cds2[c("Ccr8","Klrg1","Sell","Gzmb"),],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5)


data <- GetAssayData(Treg.intergrated,assay = "RNA",slot = "counts")
cell_metadata <- Treg.intergrated@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
saveRDS(cds, "cds_raw.rds")
cds <- readRDS("D:/4.ccr8/0719/monocle/cds_raw.rds")
#预处理，相当于Seurat流程normalize过程
#前面数据已经过normalize处理，故选择norm_method = c("none")
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("none"))
#去除批次效应,多个样本时可以通过此方式去除批次效应
cds <- align_cds(cds, alignment_group = "orig.ident")
## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds,cores=5)
## 聚类分群
cds <- cluster_cells(cds)
## 拟时序
cds <- learn_graph(cds)

##选择特定细胞作为起点，这里选择0群细胞为起点
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}

cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'seurat_clusters',my_select = "0"))

##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Treg.intergrated, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
##不同细胞类型拟时序数值，拟时序值越高表示细胞分化程度越高
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=T) + plot_cells(cds,
                                                 color_cells_by = "seurat_clusters",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=1)+ 
  scale_color_manual(values = c("#71ACD8", "#F4DA90", "#E9A66E", "#DC7656","#C14F58"))

ggsave("monocle3_1.png", width = 12 , height = 4, dpi=600)

