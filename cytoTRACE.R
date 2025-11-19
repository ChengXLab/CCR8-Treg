library(CytoTRACE2)
library(tidyverse)
library(Seurat)

cytotrace2_result_sce <- cytotrace2(Treg.intergrated, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts",
                                    species = 'mouse',
                                    seed = 1234)
cytotrace2_result_sce

saveRDS(cytotrace2_result_sce, "Treg.intergrated_cytotrace2.rds") 

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

