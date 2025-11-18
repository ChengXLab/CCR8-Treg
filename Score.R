rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(msigdbr)
library(RColorBrewer)

Treg.intergrated <- readRDS("D:/4.ccr8/R/0902/Treg.intergrated.rds")
setwd("D:/4.ccr8/R/0902/Score")

df <- read.csv('Core.csv',sep = ',',check.names = FALSE)
lst1=list()  
for(i in 1:ncol(df)) {      
  lst1[[i]] <- df[ , i]    
}
names(lst1)=colnames(df)  
print(lst1)
papermarker = lst1

#Naive
paper <- "Ccr7/Igfbp4/Tcf7/Sell"
papermarker<-str_to_title(trimws(strsplit(paper,'/')[[1]]))

#Tissue Treg
paper<-"Ccr8,Klrg1,Pparg,Tnfrsf9,Ccr5,Batf,Areg,Pdcd1,Tnfrsf4,
Itgb8,Bmyc,Cst7,Cd44,Cd69,Ms4a6b,Wls,Il2ra,Sdf4,Manbal,
Ccl2,Ptrh2,Irf6,Trf,Tff1,Il1r1,Npnt,Lman1,Clybl,Ky,
Setdb1,Ccr3,Plxnc1,Creb3l1,Alox5ap,Neb,Iscu,Tlr2,Capn2,
Cd83,Tpd52,Pdcd1lg2,Plek,Lyn,Raph1,Lgmn,Ttc39c,Rbks,Pepd,
Mettl17,Sult2b1,Tra2a,Nav2,Prmt7,Ctsc,Zfp35,Borcs5,Ccl8"
papermarker<-str_to_title(trimws(strsplit(paper,',')[[1]]))
papermarker

#ECM_KEGG
paper <-"LAMC3,CHAD,COL1A1,COL1A2,COL2A1,COL3A1,COL4A1,COL4A2,COL4A4,COL4A6,
COL5A1,COL5A2,COL6A1,COL6A2,COL6A3,COL11A1,COL11A2,COMP,COL6A6,DAG1,LAMB4,
ITGA11,SV2C,FN1,GP1BA,GP1BB,GP5,GP9,LAMA1,HMMR,HSPG2,TNC,
IBSP,ITGA6,ITGA1,ITGA2,ITGA2B,ITGA3,ITGA4,ITGA5,ITGA7,ITGA9,ITGAV,ITGB1,ITGB3,ITGB4,ITGB5,ITGB6,ITGB7,ITGB8,
AGRN,LAMA2,LAMA3,LAMA4,LAMA5,LAMB1,LAMB2,LAMB3,LAMC1,LAMC2,COL5A3,GP6,RELN,SDC1,SDC2,SDC4,TNN,SPP1,
THBS1,THBS2,THBS3,THBS4,TNR,TNXB,VTN,VWF,
ITGA10,ITGA8,CD36,CD44,CD47,SDC3,SV2B,SV2A,GZMB,SELK,CCR2"
papermarker<-str_to_title(trimws(strsplit(paper,',')[[1]]))

#Proliferation_GO: 0045787
paper <-"Akt1,Anapc5,Anp32b,Apex1,Atad5,Aurkb,Birc5,Brca1,Bub1,
Calr,Ccnb1,Ccna2,Ccnd2,Cdc6,Cdca5,Cdca8,Cdk1,Cdk4,Cenpe,
Dtl,E2f8,Eif4e,Eif4g1,Ezh2,Fbxo5,Fen1,Incenp,Mad2l1,Ncapd2,Ncaph,Ncaph2,Npm1,Nup62,Nusap1,
Pebp1,Phb2,Racgap1,Rad23a,Rad51ap1,Ranbp1,Rdx,Rrm1,Rrm2,Smc2,Spag5,Tfdp1,Trp53,Ube2c"
papermarker<-str_to_title(trimws(strsplit(paper,',')[[1]]))

#Exhaustion_Freuchet et al. Nat Immunol 2023
paper <- "PRF1、IFNG、GNLY、NKG7、GZMB、GZMA、GZMH、KLRK1、KLRB1、KLRD1、CTSW、CST7"
papermarker<-str_to_title(trimws(strsplit(paper,'、')[[1]]))



paper <- "Nr4a1、Egr1、Egr2、Myc 、 Dusp2"
papermarker<-str_to_title(trimws(strsplit(paper,'、')[[1]]))

paper <- "Nr4a1、Nr4a2、Nr4a3"
papermarker<-str_to_title(trimws(strsplit(paper,'、')[[1]]))

# 将基因转为list 
features <- list(papermarker)

#使用Seurat包自带的AddModuleScore函数打分
sce_score<- AddModuleScore(Treg.intergrated,
                           features = features,
                           ctrl = 100,
                           name = "features")
head(sce_score@meta.data)
#这里就得到了基因集评分结果，但是注意列名为 features1
colnames(sce_score@meta.data)[12] <- 'Score'
#sce_score@meta.data[1:2,1:14]

#mydata<- FetchData(sce_score,vars = c("tSNE_1","tSNE_2","Score"))
#mydata<- FetchData(sce_score,vars = c("UMAP_1","UMAP_2","Score"))
mydata<- FetchData(sce_score,vars = c("umap_1","umap_2","Score"))
#colours = c("#76489C","#718993","#74C774","#A9D156","#E4E400")
#colours = c("#76489C","#718993","#BDBEC1",'#FFCC33','#CC3333')
#colours = c("#E2E2E2","#BDBEC1","#F9C0CA", "#FC4E2A", "#E31A1C", "#BD0026")
#colours <- brewer.pal(9, "Blues")
#colours = c("#E2E2E2","#BDBEC1","#FED976", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026")
#colours = rev(brewer.pal(11, "RdBu"))RdGy

colours = c("#2166AC", "#4393C3","#D1E5F0", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

ggplot(mydata,aes(x = umap_1,y =umap_2,colour = Score))+
  geom_point(size = 2)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = colours )+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
        legend.position = "none")& NoAxes()
ggsave("Score.pdf",width = 4 , height = 4, dpi=600)


ggplot(mydata,aes(x = umap_1,y =umap_2,colour = Score))+
  geom_point(size = 2)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = colours )+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "right")& NoAxes()
ggsave("Scorekk3.pdf",width = 4 , height = 4, dpi=600)

#write.csv(sce_score@meta.data,file="score.csv")

##Vlnplot
library(ggpubr)
vln.dat=FetchData(sce_score,c("Score","orig.ident","seurat_clusters"))
vln.dat$seurat_clusters <- factor(vln.dat$seurat_clusters, 
                                  levels = c("0", "1", "2"),
                                  labels = c("C0", "C1", "C2")) # 修改后的分类


ggviolin(vln.dat, "seurat_clusters", "Score",
         fill = "seurat_clusters", #小提琴内部颜色对应的数据列
         color = "seurat_clusters", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T, 
         palette = c("#71ACD8","#F4DA90", "#E9A66E", "#60C2A6","#DC7656","#DCE49F"),  #小提琴自定义颜
         font.y = 1,  #y轴标题字体大小
         font.tickslab = c(20,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(   
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  theme(legend.position = "none")+
labs(x = NULL, #设置x轴标题
     y = "Score", #设置y轴标题
     title = NULL)
ggsave("vln_score4.pdf",width = 4, height = 2, dpi=600)

#p value
vln.dat=FetchData(sce_score,c("Score","orig.ident","seurat_clusters"))
my_comparisons=list(c("0","1"),
                    c("0","2"),
                    c("1","2"))
ggviolin(vln.dat, "seurat_clusters", "Score",
         fill = "seurat_clusters", #小提琴内部颜色对应的数据列
         color = "seurat_clusters", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T, 
         palette = c("#71ACD8", "#F4DA90", "#E9A66E","#60C2A6","#C14F58","#DC7656","#DCE49F"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(20, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 30,  #y轴标题字体大小
         font.tickslab = c(30,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(   
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  stat_compare_means(comparisons = my_comparisons)
ggsave("vln_p_score3.png",width = 8, height = 4, dpi=600)

#orig.ident
vln.dat=FetchData(sce_score,c("Score","orig.ident","seurat_clusters"))

ggviolin(vln.dat, "orig.ident", "Score",
         fill = "orig.ident", #小提琴内部颜色对应的数据列
         color = "orig.ident", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T, 
         palette = c("#FB2D18","#FFC2D2","#929D9B"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(20, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 30,  #y轴标题字体大小
         font.tickslab = c(30,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(   
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "Score", #设置y轴标题
       title = NULL)
ggsave("vln_score4.png",width = 8, height = 4, dpi=600)

##boxplot
data<- FetchData(sce_score,vars = c("orig.ident","Score"))
#write.csv(data,file="Score_data.csv")
#data <- read.csv('Score_data.csv',sep = ',',check.names = FALSE)
ggplot(data, aes(x=orig.ident,y=`Score`)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL,title = "Score")+ geom_jitter(col="#00000033", pch=19,cex=2.5, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(orig.ident)))+
  NoLegend()+theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行,标题居中

ggsave("Score_Bar.png",width = 8 , height = 6, dpi=600)


##boxplot
data<- FetchData(sce_score,vars = c("SCT_snn_res.0.5","Score"))
#write.csv(data,file="Score_data.csv")
#data <- read.csv('Score_data.csv',sep = ',',check.names = FALSE)
ggplot(data, aes(x = SCT_snn_res.0.5,y=`Score`)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL,title = "Score")+ 
  geom_jitter(col="#00000033", pch=19,cex=2.5, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(SCT_snn_res.0.5)))+
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5))& NoAxes() 

ggsave("Score_Bar_clusters.png",width = 8 , height = 6, dpi=600)





##Heart
scRNAsub = Treg.intergrated[,Treg.intergrated@meta.data$orig.ident %in% c("Heart")]
sce_score<- AddModuleScore(scRNAsub,
                           features = features,
                           ctrl = 100,
                           name = "features")
head(sce_score@meta.data)
colnames(sce_score@meta.data)[9] <- 'Score'
sce_score@meta.data[1:2,1:9]

write.csv(sce_score@meta.data,file="score_H.csv")



mydata<- FetchData(sce_score,vars = c("UMAP_1","UMAP_2","Score"))
colours = c("#76489C","#718993","#BDBEC1",'#FFCC33','#CC3333')
colours = c("#d10c00", "#e5acac", "#ffb432", "#f9d089", "#3264fc", "#91acff", "#28c2ff", "#a3e4ff", "#8851ff", "#b08cff")
cl_allCells_pal = c("#d10c00","#e5acac","#ef8753","#217dff","#6db9ff","#599958","#7c7b7b","#c6c6c6")

a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Score))+
  geom_point(size = 2)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = colours )
a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("Score_H.png",width = 8 , height = 6, dpi=600)
