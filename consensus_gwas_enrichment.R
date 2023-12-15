library(readr)
library(dplyr)
library(regioneR)
library(ggplot2)

thousand_genome_panel_hg38 = read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/thousand_genome_panel_hg38.tsv.gz")
gwas_catalog <- read_tsv("/Users/ashvinravi/Desktop/bmi206/group_project/gwas_catalog_v1.0-associations_e110_r2023-11-08.tsv")
european_gwas_catalog <- gwas_catalog %>% filter(grepl("European", `INITIAL SAMPLE SIZE`)) %>% select(CHR_ID, CHR_POS, SNPS) 

gwas_snps <- european_gwas_catalog %>% distinct(CHR_ID, CHR_POS, .keep_all = T)
gwas_snps = gwas_snps[!grepl("x", gwas_snps$CHR_POS),]
gwas_snps = gwas_snps[!grepl(";", gwas_snps$CHR_POS),]
gwas_snps = gwas_snps[!is.na(gwas_snps$CHR_POS),]
colnames(gwas_snps) = c("chrom", "start", "SNP")

gwas_snps$start = as.numeric(gwas_snps$start)
gwas_snps$end = gwas_snps$start + 1 
gwas_snps$chrom = paste0("chr", gwas_snps$chrom)
gwas_snps_maf <- inner_join(gwas_snps, thousand_genome_panel_hg38, by = c("chrom", "start", "end", "SNP"))

# Calculate Distribution Function of MAFs of GWAS gwas_snps
dist.fun = approxfun(density.default(gwas_snps_maf$MAF))

# Load DHS sites without footprints 
load("/Users/ashvinravi/Desktop/bmi206/individual_project/dhs.rda")
dhs_wout_footprints = dhs
colnames(dhs_wout_footprints) = c("chrom", "start", "end")

# Load consensus footprints
consensus_footprints = read.table("/Users/ashvinravi/Desktop/bmi206/individual_project/consensus_footprints_and_collapsed_motifs_hg38.bed.gz", sep='\t')
colnames(consensus_footprints) = c('contig', 'start', 'stop', 'identifier', 'mean_signal', 'num_samples', 'num_fps', 'width', 'summit_pos', 'core_start', 'core_end', 'motif_clusters')
consensus_footprints$posterior_prob <- 1 - (2^(-consensus_footprints$mean_signal))

consensus_0.001 = consensus_footprints[(1 - consensus_footprints$posterior_prob) <= 0.001,]
consensus_0.01 = consensus_footprints[(1 - consensus_footprints$posterior_prob) <= 0.01,]
consensus_0.05 = consensus_footprints[(1 - consensus_footprints$posterior_prob) <= 0.05,]

dhs_wout_footprints_GR = toGRanges(data.frame(chr=dhs_wout_footprints$chrom,
                                          start=dhs_wout_footprints$start, 
                                          end=dhs_wout_footprints$end))

consensus_0.001_GR = toGRanges(data.frame(chr=consensus_0.001$contig, 
                                          start=consensus_0.001$start, 
                                          end=consensus_0.001$stop))
consensus_0.01_GR = toGRanges(data.frame(chr=consensus_0.01$contig, 
                                         start=consensus_0.01$start, 
                                         end=consensus_0.01$stop)) 
consensus_0.05_GR = toGRanges(data.frame(chr=consensus_0.05$contig, 
                                         start=consensus_0.05$start, 
                                         end=consensus_0.05$stop)) 

gwas_GR = toGRanges(data.frame(chr=gwas_snps_maf$chrom, 
                               start=gwas_snps_maf$start,
                               end=(gwas_snps_maf$end)))

gwas_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, gwas_GR)
gwas_overlap_0.05 = numOverlaps(consensus_0.05_GR, gwas_GR)
gwas_overlap_0.01 = numOverlaps(consensus_0.01_GR, gwas_GR)
gwas_overlap_0.001 = numOverlaps(consensus_0.001_GR, gwas_GR)

randomized_overlap_0.05 = c()
randomized_overlap_0.01 = c()
randomized_overlap_0.001 = c()
randomized_overlap_dhs = c()

print("Generating random regions...")
for (x in 1:100) {
  # Sample Random gwas_snps from panel from 1000G panel. 
  randomized_gwas_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 172272, weight_by = dist.fun(MAF), replace=T) %>%
    ungroup()
  print(x)
  randomized_gwas_snps_GR = toGRanges(data.frame(chr=randomized_gwas_snps$chrom,
                                            start=randomized_gwas_snps$start,
                                            end=as.numeric(randomized_gwas_snps$end)))
  
  randomized_overlap_dhs = c(randomized_overlap_dhs, numOverlaps(randomized_gwas_snps_GR, dhs_wout_footprints_GR))
  randomized_overlap_0.05 = c(randomized_overlap_0.05, numOverlaps(randomized_gwas_snps_GR, consensus_0.05_GR))
  randomized_overlap_0.01 = c(randomized_overlap_0.01, numOverlaps(randomized_gwas_snps_GR, consensus_0.01_GR))
  randomized_overlap_0.001 = c(randomized_overlap_0.001, numOverlaps(randomized_gwas_snps_GR, consensus_0.001_GR))
}

enrichment_dhs = data.frame(gwas_overlap_dhs / randomized_overlap_dhs)
enrichment_dhs$category = "DHS w/out footprints"
colnames(enrichment_dhs) = c("enrichment_ratio", "Category")
enrichment_0.05 = data.frame(gwas_overlap_0.05 / randomized_overlap_0.05)
enrichment_0.05$category = "0.05"
colnames(enrichment_0.05) = c("enrichment_ratio", "Category")
enrichment_0.01 = data.frame(gwas_overlap_0.01 / randomized_overlap_0.01)
enrichment_0.01$category = "0.01"
colnames(enrichment_0.01) =  c("enrichment_ratio", "Category")
enrichment_0.001 = data.frame(gwas_overlap_0.001 / randomized_overlap_0.001)
enrichment_0.001$category = "0.001"
colnames(enrichment_0.001) = c("enrichment_ratio", "Category")
enrichment_table = rbind(enrichment_dhs, enrichment_0.05, enrichment_0.01, enrichment_0.001)
enrichment_table$Category = factor(enrichment_table$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
# write_tsv(enrichment_table, "/Users/ashvinravi/Desktop/bmi206/individual_project/enrichment_table.tsv")

boxplot <- ggplot(enrichment_table, aes(x=Category, y=enrichment_ratio, fill=Category)) + geom_boxplot() + theme_classic()

############ eQTL/sQTL analysis

thousand_genome_panel_hg38 = read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/thousand_genome_panel_hg38.tsv.gz")

randomized_overlap_0.05 = c()
randomized_overlap_0.01 = c()
randomized_overlap_0.001 = c()
randomized_overlap_dhs = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random gwas_snps from panel from 1000G panel. 
  randomized_gwas_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 172272, weight_by = dist.fun(MAF), replace=T) %>%
    ungroup()
  print(x)
  randomized_gwas_snps_GR = toGRanges(data.frame(chr=randomized_gwas_snps$chrom,
                                                 start=randomized_gwas_snps$start,
                                                 end=as.numeric(randomized_gwas_snps$end)))
  
  randomized_overlap_dhs = c(randomized_overlap_dhs, numOverlaps(randomized_gwas_snps_GR, dhs_wout_footprints_GR))
  randomized_overlap_0.05 = c(randomized_overlap_0.05, numOverlaps(randomized_gwas_snps_GR, consensus_0.05_GR))
  randomized_overlap_0.01 = c(randomized_overlap_0.01, numOverlaps(randomized_gwas_snps_GR, consensus_0.01_GR))
  randomized_overlap_0.001 = c(randomized_overlap_0.001, numOverlaps(randomized_gwas_snps_GR, consensus_0.001_GR))
}

enrichment_dhs = data.frame(gwas_overlap_dhs / randomized_overlap_dhs)
enrichment_dhs$category = "DHS w/out footprints"
colnames(enrichment_dhs) = c("enrichment_ratio", "Category")
enrichment_0.05 = data.frame(gwas_overlap_0.05 / randomized_overlap_0.05)
enrichment_0.05$category = "0.05"
colnames(enrichment_0.05) = c("enrichment_ratio", "Category")
enrichment_0.01 = data.frame(gwas_overlap_0.01 / randomized_overlap_0.01)
enrichment_0.01$category = "0.01"
colnames(enrichment_0.01) =  c("enrichment_ratio", "Category")
enrichment_0.001 = data.frame(gwas_overlap_0.001 / randomized_overlap_0.001)
enrichment_0.001$category = "0.001"
colnames(enrichment_0.001) = c("enrichment_ratio", "Category")
enrichment_table = rbind(enrichment_dhs, enrichment_0.05, enrichment_0.01, enrichment_0.001)
enrichment_table$Category = factor(enrichment_table$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table, "/Users/ashvinravi/Desktop/bmi206/individual_project/enrichment_results/consensus_gwas_enrichment.tsv")

enrichment_table = read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/enrichment_results/consensus_gwas_enrichment.tsv")
enrichment_table$Category = factor(enrichment_table$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
boxplot_gwas <- ggplot(enrichment_table, aes(x=Category, y=enrichment_ratio, fill="GWAS")) + geom_boxplot() + theme_classic() + scale_fill_manual(values = "#72bf6a") + 
  ylim(1, 2.5) + ylab("Enrichment Ratio\n") + xlab("\nConsensus DHS Footprints")



                  