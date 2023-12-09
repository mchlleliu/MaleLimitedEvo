###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                         Innocenti and Morrow X SSAV genes
# 
# 
###################################

# Enrichment test for SA candidates in Innocenti and Morrow
##########
all_genes_InnoMorrow <- tibble(all.genes)
all_genes_InnoMorrow <- all_genes_InnoMorrow %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)

test <- fisher.test(x = all_genes_InnoMorrow$Sig, y = all_genes_InnoMorrow$IsInnoMorr)

mos_plot_InnoMorr <- ggbarstats(
  all_genes_InnoMorrow, IsInnoMorr, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("SA_InnoMorr", "notSA_InnoMorr"),
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


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests//Inno_Morrow.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
mos_plot_InnoMorr
dev.off()
##########



SSAV_InnoMor <- all_genes_InnoMorrow[!is.na(all_genes_InnoMorrow$IsInnoMorr) 
                                     & all_genes_InnoMorrow$IsInnoMorr & all_genes_InnoMorrow$Sig,]
write_clip(SSAV_InnoMor$FlyBaseID) # write to clipboard for GO analysis.

## to do:
# look at sex-specific phenotype effects --> Connallon & Clark 2011