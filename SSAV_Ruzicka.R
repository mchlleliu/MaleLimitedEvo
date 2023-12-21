###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                         Ruzicka et al. 2019 X SSAV genes
# 
# 
###################################

# Packages
#########
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(ggstatsplot)
#########


# Enrichment test for SA candidates in Ruzicka et al. genes
##########
all_genes_Ruz <- tibble(all.genes)
all_genes_Ruz <- all_genes_Ruz %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)

# no difference if only looking at Chr2 either.
all_genes_Ruz_Chr2 <- merge(all_genes_Ruz, Chrs, by = "FlyBaseID", all = T)
all_genes_Ruz_Chr2 <- all_genes_Ruz_Chr2[!is.na(all_genes_Ruz_Chr2$Sig) & 
                                           all_genes_Ruz_Chr2$Chr == "2",]

test <- fisher.test(x = all_genes_Ruz_Chr2$Sig, y = all_genes_Ruz_Chr2$IsRuz)

# enrichment for all genes in SSAV males and females, vs Ruzicka all genes
mos_plot_Ruz

# enrichment for Chr 2 genes in SSAV males and females, vs Ruzicka Chr 2 genes
mos_plot_Ruz_Chr2 <- ggbarstats(
  all_genes_Ruz_Chr2[!is.na(all_genes_Ruz_Chr2$Sig) & 
                  !is.na(all_genes_Ruz_Chr2$IsRuz),], IsRuz, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("in_Ruz", "not_in_Ruz"),
                    values = c("darkorchid4", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_x_discrete(labels = c("Background", "Candidates")) +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        #legend.justification = c("right", "bottom"),
        #legend.box.just = "left",
        #legend.box.background = element_rect(),
        legend.box.background = element_rect(),
        #legend.box.margin = margin(4, 6, 6, 6),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )


##########

pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests/Ruz_mosaic.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
mos_plot_Ruz
dev.off()
