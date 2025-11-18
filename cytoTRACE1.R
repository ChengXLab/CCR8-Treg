####CytoTRACE2
#installing
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)
library(tidyverse)
library(Seurat)

# convert a v5 assay to a v4 assay
Treg.intergrated[["RNA"]] <- as(object = Treg.intergrated[["RNA"]], Class = "Assay")
#saveRDS(Treg.intergrated, "TregV4.rds") 

cytotrace2_result_sce <- cytotrace2(Treg.intergrated, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts",
                                    species = 'mouse',
                                    seed = 1234)
cytotrace2_result_sce

saveRDS(cytotrace2_result_sce, "Treg.intergrated_cytotrace2.rds") 


# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = Treg.intergrated@meta.data$orig.ident) %>% 
  set_rownames(., colnames(Treg.intergrated))

# plotting
FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1) + 
  scale_colour_gradientn(colours = (c("#5E4FA2","#66C2A5", "#E6F598", "#FEE08B", "#F46D43","#9E0142")),
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)& NoAxes()
ggsave("uamp_Cytotrace2.png",width = 4, height = 3, dpi=600)



FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1) + 
  scale_colour_gradientn(colours = (c("#5E4FA2","#66C2A5", "#E6F598", "#FEE08B", "#F46D43","#9E0142")),
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) + 
  theme(legend.position = "none",
        plot.title = element_blank())& NoAxes()
ggsave("CytoTRACE2.pdf",width = 3, height = 3, dpi=600)





####CytoTRACE
install.packages("devtools")
devtools::install_local("E:/R/R-4.3.1/CytoTRACE_0.3.3.tar.gz")

library(CytoTRACE)
Treg.intergrated <- readRDS("D:/4.ccr8/Treg/Treg.intergrated.rds")
exp1 <- as.matrix(Treg.intergrated@assays$RNA$counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)

phenot <- Treg.intergrated$seurat_clusters
phenot <- as.character(phenot)
names(phenot) <- rownames(Treg.intergrated@meta.data)
emb <- Treg.intergrated@reductions[["umap"]]@cell.embeddings

plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')
