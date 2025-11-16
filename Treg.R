library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(harmony)
setwd("D:/4.ccr8/R/0902/")

rm(list = ls())

ln <- read.csv(file = "D:/4.ccr8/lntreg.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")
Heart <- read.csv(file = "D:/4.ccr8/Htreg.csv",row.names = 1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE, sep = ",")

# Create Seurat Object
ln <- CreateSeuratObject(counts = ln, project = "ln", min.cells = 3, min.features = 200)
Heart <- CreateSeuratObject(counts = Heart, project = "Heart", min.cells = 3, min.features = 200)

#添加样本信息的MetaData
ln$sample <- "LN"
Heart$sample <- "Heart"


#####
Treg <- merge(ln,y = c(Heart),project = "Treg")
saveRDS(Treg, "Treg_all.rds") 
Treg
table(Treg$orig.ident)

##质控
scRNA <- Treg
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-") 
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^Rp[sl]")

theme.set2 = theme(axis.title.x=element_blank())
#
plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb")
group = "orig.ident"
#
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

##质控后
scRNA <- subset(scRNA, subset = nCount_RNA < 10000 & nFeature_RNA > 500 
                & nFeature_RNA < 7500 & percent.mt < 10 & percent.rb < 40 )

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8) 

###
Treg.intergrated <- NormalizeData(scRNA,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)

p1 <- DimPlot(object = Treg.intergrated,reduction = "pca",pt.size =.1,group.by = "orig.ident")
p2 <- VlnPlot(object = Treg.intergrated,features = "PC_1",pt.size =.1,group.by = "orig.ident")
plot_grid(p1,p2)
ggsave("harmony_before.png",width = 8, height = 3, dpi=600)

###
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb"))
Treg.intergrated <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)

### SCT标准化数据
Treg.intergrated <- SCTransform(Treg.intergrated)

### PCA
Treg.intergrated <- RunPCA(Treg.intergrated, npcs=50, verbose=FALSE)

#saveRDS(Treg.intergrated, "Treg.intergrated_Harmony_before.rds") 
###RunHarmony
system.time({Treg.intergrated <- RunHarmony(Treg.intergrated, group.by.vars = "orig.ident",
                                            assay.use="SCT",max.iter.harmony = 20,lambda=0.000000001)})

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

#colours_cluster = c("#71ACD8", "#F4DA90", "#E9A66E","#60C2A6","#C14F58","#DC7656","#DCE49F")
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

table(Treg.intergrated$seurat_clusters)
# Cellratio
table(Treg.intergrated$orig.ident)#查看各组细胞数
prop.table(table(Idents(Treg.intergrated)))
table(Idents(Treg.intergrated), Treg.intergrated$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(Treg.intergrated), Treg.intergrated$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = colours_cluster)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  theme(text = element_text(size = 15))
ggsave("Cellratio.pdf",width = 7, height = 2)


# 获取Ccr8基因的表达数据
ccr8_expression <- FetchData(Treg.intergrated, vars = "Ccr8")
# 将群体信息和Ccr8表达信息结合
cluster_info <- Treg.intergrated@meta.data$seurat_clusters
ccr8_data <- data.frame(cluster_info, ccr8_expression)
# 计算每个群体中Ccr8阳性细胞的比例
ccr8_positive_ratio <- function(cluster) {
  subset_data <- ccr8_data[ccr8_data$cluster_info == cluster, ]
  sum(subset_data$Ccr8 > 0) / nrow(subset_data)
}
# 计算每个群体的Ccr8表达比例
clusters <- unique(Treg.intergrated@meta.data$seurat_clusters)
ccr8_ratios <- sapply(clusters, ccr8_positive_ratio)

ccr8_ratio_df <- read.csv("Ccr8.csv",sep = ',',check.names = FALSE)
ccr8_ratio_df

ggplot(ccr8_ratio_df) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",
           width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = colours_cluster)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  theme(text = element_text(size = 15))
ggsave("ccr8_ratio_df.png",width = 4, height = 5)


# Differential gene expression analysis
library(dplyr)
# 用wilcox检验进行差异分析，耗时较长，这一步是将每一群与除这一群之外的所有细胞进行差异分析
diff.wilcox = FindAllMarkers(Treg.intergrated, test.use = "wilcox")
# 设定差异基因为p_val<0.05
#all.markers2 = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
all.markers <- subset(diff.wilcox, p_val<0.05)
# 差异基因中排名前15的基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 保存分析结果表格
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "diff_genes_wilcox_top10.csv", row.names = F)

# 用TOP10的基因做热图
top10 <- read.csv("diff_genes_wilcox_top10.csv",sep = ',',check.names = FALSE)
top10_genes <- CaseMatch(search = as.vector(top10$gene), match = rownames(Treg.intergrated)) 

top_genes <- all.markers%>%group_by(cluster)%>%top_n(n=300,wt=avg_log2FC)
DoHeatmap(Treg.intergrated,features = top_genes$gene)+NoLegend()
ggsave("heatmap2.pdf", width=6, height=12)

Heatmap <- DoHeatmap(Treg.intergrated, features = top10_genes, group.by = "seurat_clusters",
                     group.bar = T, size = 4,
                     group.colors = c("#FF7869", "#A9C181", "#B77B13", "#84A9FF"))+
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))
ggsave("heatmap.pdf",plot=Heatmap, width=6, height=12)





##Vlnplot
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

markers <- c('Ccr7',"Sell","Tcf7",'Igfbp4',
             'Gzmb','Ctla2a','Itgae', 'Lmna',
             'Ccr8','Klrg1','Areg','Il1rl1')
markers <- c("S1pr1",'Klf2')

vln.dat=FetchData(Treg.intergrated,c(markers,"orig.ident","seurat_clusters"))
vln.dat$Cell <- rownames(vln.dat)
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","seurat_clusters"), 
                               measure.vars = markers,
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(seurat_clusters,gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值
vln.dat.melt <- na.omit(vln.dat.melt)
head(vln.dat.melt,10)

vln.dat.melt$gene <- factor(vln.dat.melt$gene, levels = unique(vln.dat.melt$gene))

ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  # 分面
  facet_grid(rows = vars(seurat_clusters), scales = "free", switch = "y")+
  # 自定义填充色
  scale_fill_manual(values = my36colors) + 
  # 使用 cowplot 主题
  theme_cowplot(font_size = 20) +
  theme_bw() + theme(panel.grid=element_blank())+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=1.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  ggtitle("VlnPlot")

ggsave("VlnPlot_markers10.png", width = 8, height = 3, dpi = 600)



###Vlnplot2
#宽转长
colours_cluster = c("#71ACD8","#F4DA90", "#E9A66E", "#60C2A6","#DC7656","#DCE49F")
vln.dat=FetchData(Treg.intergrated,c(markers,"orig.ident","seurat_clusters"))
vln.dat$Cell <- rownames(vln.dat)
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("orig.ident","seurat_clusters"), 
                               measure.vars = markers,
                               variable.name = "gene", 
                               value.name = "Expr") %>%
  group_by(seurat_clusters,gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值

#write.csv(vln.dat.melt, "vln.dat.melt.csv", row.names = F)
#vln.dat.melt <- read.csv("vln.dat.melt.csv",sep = ',',check.names = FALSE)
#colours_sample = c("#E95C59","#00BFC4","#7CAE00","#F8766D")
ggplot(vln.dat.melt, aes(factor(seurat_clusters), Expr, fill = seurat_clusters)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(
    limits = c(0, 2),          # y轴范围依旧可以是 0 到 2
    breaks = c(0, 1,2),        # 只显示0和1.5
    expand = c(0, 0),
    position = "right"
  )+
  facet_grid(rows = vars(gene), scales = "free", switch = "y")+
  scale_fill_manual(values = colours_cluster) + 
  theme_cowplot(font_size = 20) +
  theme_bw() + theme(panel.grid=element_blank())+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),     # 标题居中
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),  # 全局分面文字加粗
        strip.text.y.left = element_text(angle = 0, face = "bold.italic", size = 12),  # 左侧标签：加粗 + 斜体
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  labs(y = 'Expression')+
  ggtitle("VlnPlot")

ggsave("VlnPlot_markers.pdf", width = 4.5, height = 7.5, dpi = 1200)


#Featureplot
rm(list=ls())
setwd("D:/4.ccr8/Treg")
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
# install.packages("tidydr")
# devtools::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")
# remotes::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")
#install.packages("scCustomize")
#devtools::install_github("xmc811/Scillus", ref = "development")
#library(Scillus)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

Treg.intergrated <- readRDS("D:/4.ccr8/R/0902/Treg.intergrated.rds")
markers <- c("Klf2",'Ccr7','Igfbp4',"Sell",
             'Ccr2',"Cxcr3",'Icos','Ctla2a',
             'Gzmb','Nkg7','Ifng','Pdcd1',
             'Klrg1','Ccr8','Il1rl1','Areg')

markers <- c('Ccr8','Il1r2')

i=1
plots=list()
for (i in 1:length(markers)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = Treg.intergrated, 
                                  colors_use = viridis_magma_dark_high, 
                                  features = markers[i])+NoLegend()+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
p<-wrap_plots(plots, ncol = 2);p
ggsave(p,file="Featureplot2.pdf",width = 6,height = 3.5)


for (i in 1:length(markers)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = Treg.intergrated, 
                                  colors_use = colorRampPalette(c("#3288BD", "white", "#D53E4F" ))(50), 
                                  features = markers[i])+NoLegend()+NoAxes()+
  theme(panel.border = element_rect(fill = NA,color = "black",
  size=1.5,linetype = "solid"))
}
p<-wrap_plots(plots, ncol = 4);p
ggsave(p,file="featureplot3.pdf",width = 8,height = 6)




#给图例定义范围和颜色
col_fun = colorRamp2(c(0,2,4), c("#3288BD", "white", "#D53E4F"))
lgd3 = Legend(at = c(0,2,4), col_fun = col_fun, 
              nrow = 1,#border = "black",direction = "horizontal",
              legend_height = unit(4, "cm")
)
#可以调整图例的方向，高度，还有加不加border。
#用draw把图例画在合适的位置。
draw(lgd3,x = unit(0.99, "npc"), y = unit(0.4, "npc"),
     just = c("right", "bottom"))






###气泡图
markers <- c('Ccr7',"Sell","Tcf7",'Igfbp4',
             'Ctla2a','Lmna','Gzmb',"Itgae",
             'Ccr8','Areg','Klrg1','Il1rl1',
             "Batf","Nfil3","Atf4","Blimp1",
             "Irf4")

markers <- c('Ccr7',"Sell","Tcf7",'Igfbp4',
             "Atf4",'Ctla2a','Lmna','Gzmb',"Itgae",
             "Batf",'Ccr8','Areg','Klrg1','Il1rl1')

Idents(Treg.intergrated)="seurat_clusters"

gene <- intersect(markers, rownames(Treg.intergrated))
p <- DotPlot(Treg.intergrated, features = rev(gene))
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
table(data$id)

ggplot(data, aes(x = id, y = features.plot)) +
  geom_point(
    aes(fill = avg.exp.scaled, size = pct.exp),
    color = 'black',
    shape = 21,
    stroke = 0.01) +
  labs(fill = "Z-score") +
  xlab("") + ylab("") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_size(
    range = c(0, 10),
    limits = c(0, 100),
    breaks = c(0,20,40,60,80,100)) +
  scale_x_discrete(labels = function(x) paste0("C", x)) +  # x轴改为C0, C1...
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 四周边框
    axis.line = element_blank(),  # 防止重复边框
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5,
                               size = 18, face = "bold"),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1,
                               size =18,face = "bold.italic" ),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),   # 图例刻度字体大小
    aspect.ratio = 3
  ) +
  guides(
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      ncol = 1,
      byrow = TRUE,
      override.aes = list(stroke = 0.4)
    ),
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5
    )
  )

ggsave("dotplot2.pdf", width = 6, height = 9, dpi = 600)


## scRNA_Vocalno
library(ggrepel)
library(tidyverse)

cells1 <- subset(Treg.intergrated@meta.data, seurat_clusters %in%  c("2"))  %>% rownames()
cells2 <- subset(Treg.intergrated@meta.data, seurat_clusters %in%  c("0"))  %>% rownames()
deg <- FindMarkers(Treg.intergrated, ident.1 = cells1, ident.2 = cells2)
write.csv(deg, "deg.csv", row.names = T)

deg <- read.csv('deg.csv',sep = ',',check.names = FALSE)
deg <- deg[,-1]
deg <- deg %>% 
  mutate(type = case_when(avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
                          abs(avg_log2FC) < 0.25 | p_val_adj > 0.05 ~ "no",
                          avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down"))
deg   <- deg %>%
  mutate(Difference = pct.1 - pct.2) 

deg$gene <- rownames(deg)

#markers <- deg %>% filter(gene %in% c("Plin2"))
markers <- deg %>% filter(gene %in% c("Ccr8",'Tnfrsf9','Pdcd1','Klrg1','Icos','Areg','Batf'))

ggplot(deg, aes(x=Difference, y=avg_log2FC, color = type)) + 
  geom_point(size=1) +
  theme_bw()+
  scale_color_manual(values=c("grey","grey","grey"))+
  # 图例
  theme(plot.title = element_text(hjust = 0.5,size=10),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
  )+
  geom_point(data = markers,color = "#CC3300",size = 1)+
  geom_text_repel(data=markers, aes(label=gene), color="black", 
                  size=5, fontface="italic", 
                  point.padding = 0.3, box.padding = 1,
                  segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  ylab('Log-Fold Change')+xlab('Percentage Difference')

ggsave("Vocalno.png",width = 6, height =4, dpi = 600)


## scRNA_Vocalno2

ggplot(deg,aes(avg_log2FC,-log10(p_val_adj),fill=type))+
  geom_point(shape=21,aes(size=-log10(p_val_adj)))+
  scale_fill_manual(values=c('seagreen','gray','orange'))+
  scale_color_manual(values=c('gray60','black'))+
  geom_vline(xintercept=c(-0.25,0.25),lty=2,col="gray30",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="gray30",lwd=0.6)+
  theme_bw(base_rect_size = 1)+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        plot.title = element_text(family = 'regular',hjust = 0.5),
        legend.position = c(0.5, 1),
        legend.justification = c(0.5, 1),
        legend.key.height = unit(0.5,'cm'),
        legend.background = element_rect(fill = NULL, colour = "black",size = 0.5))+
  geom_text_repel(data=markers, aes(label=rownames(markers)), color="black", 
                  size=6, fontface="italic", 
                  point.padding = 0.3, box.padding = 1,
                  segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)+
  xlim(-2,2)+
  guides(size=F,color=F)+
  ylab('-log10 (FDR)')+xlab('log2 (Fold Change)')



## scRNA_Vocalno3
library(scRNAtoolVis)
library(dplyr)
setwd("D:/4.ccr8/0719/")
diff.wilcox = FindAllMarkers(Treg.intergrated, test.use = "wilcox")
all.markers <- subset(diff.wilcox, p_val<0.05)

markers <- c("Klf2",'Ccr7','Igfbp4',"Sell",
             'Ccr2',"Cxcr3",'Gzmb','Ctla2a',
             'Klrg1','Ccr8','Il1rl1','Areg')
  
jjVolcano(diffData = all.markers,
          myMarkers = markers)

jjVolcano(diffData = all.markers,
          myMarkers = markers,
          tile.col = corrplot::COL2('RdBu', 15)[4:12],
          size  = 3.5,
          fontface = 'italic')
ggsave("Vocalno_all.png",width = 6, height =4, dpi = 600)





