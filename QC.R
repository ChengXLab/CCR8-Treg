library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(harmony)

setwd("D:/ccr8/R/")

rm(list = ls())

ln <- read.csv(file = "D:/ccr8/lntreg.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")
Heart <- read.csv(file = "D:/4.ccr8/Htreg.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")

# Create Seurat Object
ln <- CreateSeuratObject(counts = ln, project = "ln", min.cells = 3, min.features = 200)
Heart <- CreateSeuratObject(counts = Heart, project = "Heart", min.cells = 3, min.features = 200)

#MetaData
ln$sample <- "LN"
Heart$sample <- "Heart"

#Merge
Treg <- merge(ln,y = c(Heart),project = "Treg")
saveRDS(Treg, "Treg_all.rds") 
Treg
table(Treg$orig.ident)

##质控
Treg[["percent.mt"]] <- PercentageFeatureSet(Treg, pattern = "^mt-")
Treg[["percent.rb"]] <- PercentageFeatureSet(Treg, pattern = "^Rp[sl]")
Treg[["complexity"]] <- Treg$nFeature_RNA / Treg$nCount_RNA

dir.create("QC")
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb")
group = "orig.ident"

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)


##质控后
Treg <- subset(Treg, subset = nCount_RNA < 10000 & nFeature_RNA > 500 
                & nFeature_RNA < 7500 & percent.mt < 10 & percent.rb < 40 )

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8) 


### SCT
Treg.intergrated <- SCTransform(Treg)

### PCA
Treg.intergrated <- RunPCA(Treg.intergrated, npcs=50, verbose=FALSE)

### RunHarmony
system.time({Treg.intergrated <- RunHarmony(Treg.intergrated, group.by.vars = "orig.ident",
                                            assay.use="SCT",max.iter.harmony = 20,lambda=1)})

p1 <- DimPlot(object = Treg.intergrated,reduction = "harmony",pt.size =.1,group.by = "orig.ident")
p2 <- VlnPlot(object = Treg.intergrated,features = "harmony_1",pt.size =.1,group.by = "orig.ident")
plot_grid(p1,p2)
ggsave("harmony_after.png",width = 8, height = 3, dpi=600)

ElbowPlot(Treg.intergrated, ndims = 50)
pc.num=1:20

Treg.intergrated <- FindNeighbors(Treg.intergrated, reduction = "harmony", dims = pc.num)
Treg.intergrated <- FindClusters(Treg.intergrated, resolution = 0.6)
Treg.intergrated <- RunUMAP(Treg.intergrated, reduction = "harmony", dims = pc.num)
#Treg.intergrated <- RunTSNE(Treg.intergrated, reduction = "harmony", dims = pc.num)

FeaturePlot(Treg.intergrated,reduction = "umap", features = c('Foxp3','Ccr8'), 
            ncol = 2, pt.size = 1.3)& NoAxes()
ggsave("Ccr8.png",width = 7, height = 3.5, dpi=600)

colours_cluster = c("#71ACD8","#F4DA90", "#E9A66E","#C14F58", "#60C2A6","#DC7656","#DCE49F")
DimPlot(Treg.intergrated, reduction = "umap", label = F, pt.size = 0.8,
        cols =colours_cluster) & NoAxes()& NoLegend()
ggsave("umap.pdf",width = 3, height = 3, dpi=600)

DimPlot(Treg.intergrated, split.by = "orig.ident", reduction='umap',label= F,ncol = 2, pt.size = 0.8,
        cols =colours_cluster) & NoAxes()& NoLegend()
ggsave("umap_split.pdf",width = 6, height = 3.2, dpi=600)

colours_sample = c("#71ACD8", "#F4DA90")
DimPlot(Treg.intergrated, group.by = "orig.ident", reduction = "umap", label = F, pt.size = 0.8,
        cols =colours_sample) & NoAxes()
ggsave("umap_sample.png",width = 4, height = 3, dpi=600)

saveRDS(Treg.intergrated, "Treg.intergrated.rds") 





