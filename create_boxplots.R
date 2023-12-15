library(ggplot2)
library(readr)
library(dplyr)
library(rstatix)

enrichment_blood_gwas <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/blood_gwas_enrichment_table.tsv")
enrichment_blood_eQTL <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/blood_eQTL_enrichment_table.tsv")
enrichment_blood_sQTL <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/blood_sQTL_enrichment_table.tsv")

enrichment_blood_gwas$Variant <- "GWAS"
enrichment_blood_eQTL$Variant <- "eQTL"
enrichment_blood_sQTL$Variant <- "sQTL"

enrichment_table = rbind(enrichment_blood_gwas, enrichment_blood_eQTL, enrichment_blood_sQTL)

enrichment_table$Category = factor(enrichment_table$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
enrichment_table$Variant <- factor(enrichment_table$Variant, levels = c("GWAS", "sQTL",  "eQTL"))

boxplot <- ggplot(enrichment_table, aes(x=Category, y=enrichment_ratio, fill=Variant)) + 
  geom_boxplot(outlier.shape = NA) +
  ggtitle("Enrichment in Lymphocyte DHS Footprints") + 
  xlab("\nLymphocyte DHS Footprints") + 
  ylab("Enrichment Ratio\n") + 
  scale_fill_manual(name = "Variants", labels = c("White Blood Cell \nCount GWAS", "Whole Blood sQTL", "Whole Blood eQTL"), values = c("#BCAB79", "#A53F2B", "#e17f3c")) +
  theme_classic() + 
  ylim(0,16) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))

####### ENRICHMENT LUNG PLOT

enrichment_lung_gwas <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/lung_gwas_enrichment_table.tsv")
enrichment_lung_eQTL <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/lung_eQTL_enrichment_table.tsv")
enrichment_lung_sQTL <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/lung_sQTL_enrichment_table.tsv")

enrichment_lung_gwas$Variant <- "GWAS"
enrichment_lung_eQTL$Variant <- "eQTL"
enrichment_lung_sQTL$Variant <- "sQTL"

enrichment_table = rbind(enrichment_lung_gwas, enrichment_lung_eQTL, enrichment_lung_sQTL)

enrichment_table$Category = factor(enrichment_table$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
enrichment_table$Variant <- factor(enrichment_table$Variant, levels = c("GWAS", "sQTL", "eQTL"))

boxplot <- ggplot(enrichment_table, aes(x=Category, y=enrichment_ratio, fill=Variant)) + 
  geom_boxplot(outlier.shape = NA) +
  ggtitle("Enrichment in Lung DHS Footprints") + 
  xlab("\nLung DHS Footprints") + 
  ylab("Enrichment Ratio\n") + 
  scale_fill_manual(name = "Variants", labels = c("Pulmonary\nFunction GWAS", "Lung sQTL", "Lung eQTL"), values = c("#5d8ee4", "#bbc1c6", "#7757c5")) +
  theme_classic() + 
  ylim(0,16) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))

dunn_test(enrichment_ratio ~ Category, data = enrichment_blood_eQTL)


####### BLOOD eQTL vs. GWAS boxplot 

enrichment_lung_eQTL$Category = factor(enrichment_lung_eQTL$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))

eQTL_boxplot <- ggplot(enrichment_lung_eQTL, aes(x=Category, y=enrichment_ratio, fill=Variant)) + 
  geom_boxplot(outlier.shape = NA) +
  ggtitle("Enrichment in Lung DHS Footprints") + 
  xlab("\nLung DHS Footprints") + 
  ylab("Enrichment Ratio\n") + 
  scale_fill_manual(name = "Variants", labels = c("Lung eQTL"), values = c("#C0AFE2")) +
  theme_classic() + 
  ylim(0,16) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))

enrichment_lung_gwas$Category = factor(enrichment_lung_gwas$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
gwas_boxplot <- ggplot(enrichment_lung_gwas, aes(x=Category, y=enrichment_ratio, fill=Variant)) + 
  geom_boxplot(outlier.shape = NA) +
  ggtitle("Enrichment in Lung DHS Footprints") + 
  xlab("\nLung DHS Footprints") + 
  ylab("Enrichment Ratio\n") + 
  scale_fill_manual(name = "Variants", labels = c("Lung eQTL"), values = c("#5d8ee4")) +
  theme_classic() + 
  ylim(0,16) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text = element_text(size = 11, color = "black"), 
        axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))

dunn_test(enrichment_ratio ~ Category, data = enrichment_lung_eQTL)


