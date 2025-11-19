library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(harmony)

setwd("D:/ccr8/R/")

rm(list = ls())

ln <- read.csv(file = "D:/ccr8/lnscRNA.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")
Heart <- read.csv(file = "D:/4.ccr8/HscRNA.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")

# Create Seurat Object
ln <- CreateSeuratObject(counts = ln, project = "ln", min.cells = 3, min.features = 200)
Heart <- CreateSeuratObject(counts = Heart, project = "Heart", min.cells = 3, min.features = 200)

#MetaData
ln$sample <- "LN"
Heart$sample <- "Heart"

#Merge
scRNA <- merge(ln,y = c(Heart),project = "scRNA")
saveRDS(scRNA, "scRNA_all.rds") 
scRNA
table(scRNA$orig.ident)

##
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^Rp[sl]")
scRNA$Complexity <- log10(scRNA$nFeature_RNA / scRNA$nCount_RNA)

sce <- scRNA
sce <- as.SingleCellExperiment(sce)
sce <- scDblFinder(sce)
scRNA$Doublet_score <- sce$scDblFinder.score
scRNA$scDblFinder_class <- sce$scDblFinder.class

theme.set2 <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y  = element_text(size = 10),
    strip.text   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.position = "none")
group = "orig.ident"
plot.features = c("nFeature_RNA", "nCount_RNA","percent.mt", 
                   "Complexity","Doublet_score")
plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA,
                       group.by = group,
                       pt.size = 0,
                       features = plot.features[i],
                       cols = c("#71ACD8","#F4DA90" )) + 
    theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, ncol = 5)
dir.create("QC", showWarnings = FALSE)
ggsave("QC/vlnplot_before_qc_all_metrics_H.pdf", 
       plot = violin, width = 8, height = 2.5)

## QC---------------------
scRNA <- subset(
  scRNA,
  subset = 
    nFeature_RNA > 500 &
    nFeature_RNA < 5000 &
    nCount_RNA < 30000 &
    percent.mt < 5 &
    percent.rb < 50 &
    scDblFinder_class == "singlet")

plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA,
                       group.by = group,
                       pt.size = 0,
                       features = plot.features[i],
                       cols = c("#71ACD8","#F4DA90")) + 
    theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, ncol = 5)
ggsave("QC/vlnplot_after_qc_all_metrics_H.pdf",
       plot = violin, width = 8, height = 2.5)


### SCT
scRNA.intergrated <- SCTransform(scRNA)

### PCA
scRNA.intergrated <- RunPCA(scRNA.intergrated, npcs=50, verbose=FALSE)

### RunHarmony
system.time({scRNA.intergrated <- RunHarmony(scRNA.intergrated, group.by.vars = "orig.ident",
                                            assay.use="SCT",max.iter.harmony = 20,lambda=1)})

ElbowPlot(scRNA.intergrated, ndims = 50)
pc.num=1:10

scRNA.intergrated <- FindNeighbors(scRNA.intergrated, reduction = "harmony", dims = pc.num)
scRNA.intergrated <- FindClusters(scRNA.intergrated, resolution = 0.6)
scRNA.intergrated <- RunUMAP(scRNA.intergrated, reduction = "harmony", dims = pc.num)

colours_cluster = c("#71ACD8","#F4DA90", "#E9A66E","#C14F58", "#60C2A6","#DC7656","#DCE49F")
DimPlot(scRNA.intergrated, reduction = "umap", label = F, pt.size = 0.8,
        cols =colours_cluster) & NoAxes()& NoLegend()
ggsave("umap.pdf",width = 3, height = 3, dpi=600)

saveRDS(scRNA.intergrated, "scRNA.intergrated.rds") 







