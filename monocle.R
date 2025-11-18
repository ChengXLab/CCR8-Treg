BiocManager::install("monocle")
setwd("D:/4.ccr8/R/0902/monocle")
library(monocle)
#读入聚类好的object
rm(list = ls())
###准备数据
Treg.intergrated <- readRDS("Treg.intergrated.rds")
#选择感兴趣的细胞群,不运行
#Treg.intergrated <- subset(Treg.intergrated,ident = 0,1,2)

# convert a v5 assay to a v4 assay
Treg.intergrated[["RNA"]] <- as(object = Treg.intergrated[["RNA"]], Class = "Assay")

#提取表达矩阵
matrix <- Treg.intergrated@assays$RNA@data
matrix <- as.matrix(matrix)
#导入分群信息
pd <- new("AnnotatedDataFrame",data = Treg.intergrated@meta.data)
fdata <- data.frame(gene_short_name = rownames(matrix),row.names = row.names(matrix))
fd <- new("AnnotatedDataFrame",data = fdata)
#建立monocle object
monocle <- newCellDataSet(matrix,phenoData = pd,featureData = fd,
                          lowerDetectionLimit = 0.2,
                          expressionFamily = negbinomial.size())
#分析细胞间离散度
monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)
#分析群与群之间差异表达的基因，作为拟时依据
monocle <- detectGenes(monocle, min_expr = 0.05)
expressed_genes <- row.names(subset(fData(monocle),
                                    num_cells_expressed >= 500))
#saveRDS(monocle, "monocle_raw.rds")

#进行差异基因的选择
print(head(pData(monocle)))
diff_test_res <- differentialGeneTest(monocle[expressed_genes,],fullModelFormulaStr="~seurat_clusters")
write.csv(diff_test_res,file="diff_test_res_500.csv")
#diff_test_res <- read.csv('diff_test_res.csv',sep = ',',check.names = FALSE)
#dim(diff_test_res); table(duplicated(diff_test_res$name))
#diff_test_res <- diff_test_res[!duplicated(diff_test_res$name),]
#rownames(diff_test_res) <- diff_test_res[,1]

#选择p值最小的500个基因
monocle_ordering_genes <-
  row.names(diff_test_res)[order(diff_test_res$qval)][0:2000]

#monocle_ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

#根据差异基因运行拟时算法，总结进化规律，不需更改
monocle <- setOrderingFilter(monocle,
                             ordering_genes = monocle_ordering_genes)
#降维
monocle <- reduceDimension(monocle, method = 'DDRTree',
                           max_components = 2,
                           num_dim = 6,
                           reduction_formula = "~ orig.ident")
#排布细胞，慢
monocle <- orderCells(monocle)

plot_cell_trajectory(monocle, cell_size = 3 , color_by = "State")

ggsave("monocle_1000.png",width = 7 , height = 5, dpi=600)

#monocle <- orderCells(monocle, root_state = 3)

plot_cell_trajectory(monocle, cell_size = 3, color_by = "Pseudotime") 
ggsave("monocle_Pseu.png",width = 5 , height = 4, dpi=600)

plot_cell_trajectory(monocle, cell_size = 3 , color_by = "orig.ident")+
  facet_wrap(~orig.ident, nrow = 2)
#结果已出，先保存结果，下面开始作图
saveRDS(monocle, "monocle_500_2000.rds") 
colnames(monocle@phenoData@data)
monocle <- readRDS("E:/scRNA/HSLB-Treg/monocle.rds")

#colours_sample = c( "#82B0D2","#FFC5AC", "#77D822","#f88421"  )
#colours_sample = c( "#66B8FF","#7CAE00","#F8766D", "#C77CFF")
colours_sample = c("#F7B59F","#BCB7D5","#8DCCBF","#F7F6BC")
plot_cell_trajectory(monocle, cell_size = 3 , color_by = "orig.ident")+
  scale_color_manual(values=colours_sample) & NoAxes()
ggsave("monocle_sample.png",width = 5 , height = 4, dpi=600)
 
plot_cell_trajectory(monocle, cell_size = 3, color_by = "Pseudotime") & NoAxes()
ggsave("monocle_Pseu.png",width = 5 , height = 4, dpi=600)

plot_cell_trajectory(monocle, cell_size = 3 , color_by = "SCT_snn_res.0.5")+
  scale_color_manual(values=colours_cluster) & NoAxes()
ggsave("monocle_cluster.png",width = 4 , height = 3, dpi=600)


#分开按sample展示
meta <- monocle@phenoData@data
meta <- data.frame(meta,t(monocle@reducedDimS))
type <- paste(meta$orig.ident,meta$orig.ident,sep = "-")
meta$type <- type
ggplot(meta,aes(x=X1,y=X2,colour= orig.ident)) +
  geom_point(size=1,shape =19) +
  facet_wrap( ~ type,ncol = 1)+
  scale_color_manual(values=colours_sample) +
  theme_bw() + 
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text =element_blank())+
  theme(panel.border = element_rect(fill = NA,color = "black"))+
  theme(panel.grid=element_blank())
ggsave("monocle_sample_split.png", width = 3.2, height = 8, dpi=600)


#分开按cluster
meta <- monocle@phenoData@data
meta <- data.frame(meta,t(monocle@reducedDimS))
type <- paste(meta$orig.ident,meta$SCT_snn_res.0.5,sep = "-")
meta$type <- type
#colours = c( "#D6E4FF","#FFC5AC", "#0A7D58","#93790D"  )
#colours_cluster = c("#FF7869", "#A9C181", "#B77B13", "#84A9FF")
#colours_cluster = c("#FFA77C", "#A9C181", "#84A9FF","#C14F58","#60C2A6")
colours_cluster = c("#71ACD8", "#F4DA90", "#E9A66E", "#DC7656","#C14F58","#60C2A6","#DCE49F")
#write.csv(meta,file="meta.csv")
#meta <- read.csv('meta.csv',sep = ',',check.names = FALSE)
ggplot(meta,aes(x=X1,y=X2,colour= SCT_snn_res.0.5)) +
  geom_point(size=1,shape =19) +
  facet_wrap( ~ type, ncol = 3)+
  scale_color_manual(values=colours_cluster)+
  theme_bw() + 
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text =element_blank())+
  theme(panel.grid=element_blank())
ggsave("monocle_sample_split_cluster.png", width = 8, height = 8, dpi=600)

table(Treg.intergrated@meta.data$seurat_clusters)

#指定基因
colnames(monocle@phenoData@data)
markers <- c("Ccr8","Klrg1","Dgat2","rna_Pparg","Tnfrsf9",
             "Ccr5","Batf","Il2ra","Icos","Neb","Il1rl1",
             "Itgb8","Itgav","Areg","Ms4a6b","Itgae","Ctla4",
             "Bmyc","Tff1","Maf","Glrx","Lyn")

s.genes <- c("Ccr8","Itgae","Batf","Areg",
             "Klrg1","Tnfrsf9","Lyn","Itgav","Ccr5")

s.genes <- c("Ccr8")
plot_genes_in_Pseudotime(monocle[s.genes,], cell_size = 1,color_by = "seurat_clusters",ncol=2)+
  scale_color_manual(values=colours_cluster)+
  theme(legend.position = "none")+
  NoAxes()
ggsave("monocle_Genes_Igfbp4.png", width = 4, height = 4, dpi=600)

s.genes <- c("Ccr8","Tnfrsf9","Sell","Klf2")
plot_genes_in_Pseudotime(monocle[s.genes,], cell_size = 2,color_by = "seurat_clusters",ncol=2)+
  scale_color_manual(values=colours_cluster)& NoAxes()

ggsave("monocle_Genes_set.png", width = 8, height = 6, dpi=800)

#指定基因2
colnames(pData(monocle))
pData(monocle)$Ccr8 = log2(exprs(monocle)['Ccr8',]+1)
plot_cell_trajectory(monocle, color_by = "Ccr8")+
  scale_color_distiller(palette = "Spectral")
ggsave("Genes_Ccr8.png", width = 7 , height = 5, dpi=600)

#plot
meta_Ccr8 <- data.frame(meta,t(monocle@reducedDimS))
type <- paste(clusters$sample,clusters$sample,sep = "-")
meta_Ccr8$type <- type
meta_Ccr8$meta <- monocle@phenoData@data$Ccr8
#write.csv(meta,file="meta.csv")
#clusters <- read.csv('meta.csv',sep = ',',check.names = FALSE)
ggplot(meta_Ccr8,aes(x=X1,y=X2,colour= meta)) +
  geom_point(size=1,shape =19) +
  facet_wrap( ~ type,ncol = 4)+
  scale_color_distiller(palette = "Spectral")+
  theme_bw() + 
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text =element_blank())+
  theme(panel.grid=element_blank())
ggsave("plot_sample_split_cluster_Ccr8_.png", width = 10, height =2, dpi=600)


###山脊图
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)

plotdf=pData(monocle)

colours_cluster = c("#71ACD8", "#F4DA90", "#E9A66E", "#DC7656","#C14F58","#60C2A6","#DCE49F")
ggplot(plotdf, aes(x=Pseudotime,y=seurat_clusters,fill=seurat_clusters))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_fill_manual(values = colours_cluster)+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggsave("monocle_tmp.png",width = 13,height = 7,units = "cm", dpi=600)

write.csv(plotdf,file="plotdf.csv")
plotdf <- read.csv('plotdf.csv',sep = ',',check.names = FALSE)

colours_sample = c("#FB2D18","#FFC2D2","#929D9B")
ggplot(plotdf, aes(x=Pseudotime,y=orig.ident,fill=orig.ident))+
  geom_density_ridges(scale=1)+
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  scale_fill_manual(values = colours_sample)+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggsave("monocle_tmp_sample.png",width = 13,height = 7,units = "cm", dpi=600)


#0719fenqun
meta <- monocle@phenoData@data
meta <- data.frame(meta,t(monocle@reducedDimS))
type <- paste(meta$orig.ident,meta$SCT_snn_res.0.5,sep = "-")
meta$type <- type
write.csv(meta,file="meta0724.csv")
meta <- read.csv('meta0724.csv',sep = ',',check.names = FALSE)
ggplot(meta, aes(x=Pseudotime,y=type,fill=type))+
  geom_density_ridges(scale=1)+
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  scale_fill_manual(values = colours_cluster)+
  theme_minimal()+
  theme_bw() + 
  theme(panel.grid = element_blank())
ggsave("monocle_tmp_split_c2.png",width = 13,height = 9,units = "cm", dpi=600)

dev.off()



##heatmap
diff_test_res <- read.csv('diff_test_res.csv',sep = ',',check.names = FALSE)
rownames(diff_test_res) <- diff_test_res[,1]
diff_test_res <- diff_test_res[,-1]
monocle_ordering_genes <-row.names(diff_test_res)[order(diff_test_res$qval)][0:500]
Time_diff <- differentialGeneTest(monocle[monocle_ordering_genes,], cores = 1, 
                                     fullModelFormulaStr = "~sm.ns(Pseudotime)")
#Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
#write.csv(Time_diff, "monocle_Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p1 <- plot_pseudotime_heatmap(monocle[Time_genes,], 
                        num_clusters=2, show_rownames=F,return_heatmap=T)
#genes <- c("Ccr8","Tnfrsf9","Batf","Areg","Klf2","Sell","Ccr7","Klrg1","Il1rl1")

ggsave("monocle_Time_heatmapAll.pdf",p1, width = 3.1, height = 30, dpi=600)


plot_pseudotime_heatmap(monocle[Time_genes,],
                        cluster_rows = TRUE,
                        cores = 1,
                        show_rownames = TRUE)



library(ClusterGVis)
diff_test_res <- read.csv('diff_test_res_200.csv',sep = ',',check.names = FALSE)
rownames(diff_test_res) <- diff_test_res[,1]
diff_test_res <- diff_test_res[,-1]
monocle_ordering_genes <-row.names(diff_test_res)[order(diff_test_res$qval)][0:500]
df <- plot_pseudotime_heatmap2(monocle[monocle_ordering_genes[1:450],],
                               num_clusters = 2,
                               cores = 1)
svg("heatmap2.svg" ,          
    width = 7, height = 7, pointsize = 12,          
    onefile = FALSE, family = "sans", bg = "white",          
    antialias = c("default", "none", "gray", "subpixel"))
gene <- c("Ccr8","Areg","Klrg1","Il1rl1","Ccr7","Sell","Tcf7","Igfbp4")
#gene = sample(df$wide.res$gene,20,replace = F)
visCluster(object = df,plot.type = "heatmap",
           markGenes = gene,
           pseudotime_col = c("white", "#E41A1C"))
dev.off()
# test code
df <- plot_genes_branched_heatmap2(monocle[row.names(subset(diff_test_res,qval < 1e-4)),],                                   
                                   branch_point = 1,                                   
                                   num_clusters = 2,                                   
                                   cores = 1,                                   
                                   use_gene_short_name = T,                                   
                                   show_rownames = F)

#Treg.intergrated@meta.data <- Treg.intergrated@meta.data[, -16]

visCluster(object = df,plot.type = "heatmap")



library(Seurat)
library(monocle)
library(tidyverse)
#Treg.intergrated <- readRDS("D:/4.ccr8/R/0902/Treg.intergrated.rds")
Treg.intergrated <- readRDS("D:/4.ccr8/R/0802/Monocle/Treg.intergrated3.rds")
paper <- "Ccr8"
papermarker<-str_to_title(trimws(strsplit(paper,'/')[[1]]))
features <- list(papermarker)
#使用Seurat包自带的AddModuleScore函数打分
Treg.intergrated<- AddModuleScore(Treg.intergrated,
                           features = features,
                           ctrl = 100,
                           name = "features")
head(Treg.intergrated@meta.data)
#这里就得到了基因集评分结果，但是注意列名为 features1
colnames(Treg.intergrated@meta.data)[16] <- 'Score'

#subset monocle
pData(monocle)[,"CB"]=rownames(pData(monocle))
anno.from.monocle = pData(monocle)
anno.from.monocle = anno.from.monocle[,c("Pseudotime","State","CB")]
#anno.from.monocle$CB <- paste(substring(anno.from.monocle$CB,3,18))

#subset Seurat@meta.data
Treg.intergrated@meta.data$CB = rownames(Treg.intergrated@meta.data)
metadata = Treg.intergrated@meta.data

#plotdata
metadata = metadata %>% 
  inner_join(anno.from.monocle,by = "CB")
#metadata$CB <- paste(substring(metadata$CB,3,18))
metadata <- metadata[!duplicated(metadata$CB), ]
rownames(metadata) = metadata$CB
plotdata = metadata[,c("Pseudotime","Score","State","seurat_clusters")]

colours_cluster = c("#71ACD8","#E9A66E", "#F4DA90",  "#DC7656","#C14F58","#60C2A6","#DCE49F")

ggplot(plotdata, aes(x = Pseudotime, y = Score)) +
  geom_point(aes(color = seurat_clusters), size = 1.5) +  # 添加透明度让点更美观
  stat_smooth(se = FALSE, color = "black", method = "loess", span = 0.75, linewidth = .8) +
  scale_x_continuous("Pseudotime", limits = c(-0, NA)) +
  scale_y_continuous("Score", breaks = pretty(plotdata$Score, n = 5)) + 
  scale_color_manual("seurat_clusters", values = colours_cluster) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 11, face = "bold"),  # 坐标轴标题
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # 调整边框宽度
    panel.grid.major = element_line(color = "grey90", linewidth = 0.1),  # 主要网格线宽度
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),  # 次要网格线宽度
    legend.position = "none"  # 控制图例
  )
ggsave("monocle_Genes_Ccr8.pdf", width = 2.8, height = 2.8, dpi=600)

