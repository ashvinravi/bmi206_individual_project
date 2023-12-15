library(readr)
library(dplyr)
library(regioneR)
library(ggplot2)

set.seed(252)
thousand_genome_panel_hg38 = read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/thousand_genome_panel_hg38.tsv.gz")

list_of_eQTL_files = list.files("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_eQTL/", pattern = "egenes.txt.gz")
setwd("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_eQTL/")
consensus_eQTL_variants = data.frame()
for (file in list_of_eQTL_files) {
  eGenes = read_tsv(file)
  print(file)
  significant_variants = eGenes %>% 
    filter(qval <= 0.05) %>% 
    select(chr, variant_pos, rs_id_dbSNP151_GRCh38p7) %>% 
    distinct()
  consensus_eQTL_variants = rbind(consensus_eQTL_variants, significant_variants)
}

consensus_eQTL_variants <- consensus_eQTL_variants %>% distinct()
colnames(consensus_eQTL_variants) = c("chrom", "start", "SNP")
consensus_eQTL_variants$end <- consensus_eQTL_variants$start + 1
### Join by minor allele frequency
eQTL_snps_maf <- inner_join(consensus_eQTL_variants, thousand_genome_panel_hg38, by = c("chrom", "start", "end", "SNP"))
dist.fun.consensus.eQTL = approxfun(density.default(eQTL_snps_maf$MAF))

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

eQTL_GR = toGRanges(data.frame(chr=eQTL_snps_maf$chrom, 
                               start=eQTL_snps_maf$start,
                               end=(eQTL_snps_maf$end)))

eQTL_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, eQTL_GR)
eQTL_overlap_0.05 = numOverlaps(consensus_0.05_GR, eQTL_GR)
eQTL_overlap_0.01 = numOverlaps(consensus_0.01_GR, eQTL_GR)
eQTL_overlap_0.001 = numOverlaps(consensus_0.001_GR, eQTL_GR)

randomized_overlap_0.05_eQTL = c()
randomized_overlap_0.01_eQTL = c()
randomized_overlap_0.001_eQTL = c()
randomized_overlap_dhs_eQTL = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random eQTL_snps from panel from 1000G panel. 
  randomized_eQTL_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 217813, weight_by = dist.fun.consensus.eQTL(MAF), replace=T) %>%
    ungroup()
  print(x)
  randomized_eQTL_snps_GR = toGRanges(data.frame(chr=randomized_eQTL_snps$chrom,
                                                 start=randomized_eQTL_snps$start,
                                                 end=as.numeric(randomized_eQTL_snps$end)))
  
  randomized_overlap_dhs_eQTL = c(randomized_overlap_dhs_eQTL, numOverlaps(randomized_eQTL_snps_GR, dhs_wout_footprints_GR))
  randomized_overlap_0.05_eQTL = c(randomized_overlap_0.05_eQTL, numOverlaps(randomized_eQTL_snps_GR, consensus_0.05_GR))
  randomized_overlap_0.01_eQTL = c(randomized_overlap_0.01_eQTL, numOverlaps(randomized_eQTL_snps_GR, consensus_0.01_GR))
  randomized_overlap_0.001_eQTL = c(randomized_overlap_0.001_eQTL, numOverlaps(randomized_eQTL_snps_GR, consensus_0.001_GR))
}

enrichment_dhs_eQTL = data.frame(eQTL_overlap_dhs / randomized_overlap_dhs_eQTL)
enrichment_dhs_eQTL$category = "DHS w/out footprints"
colnames(enrichment_dhs_eQTL) = c("enrichment_ratio", "Category")
enrichment_0.05_eQTL = data.frame(eQTL_overlap_0.05 / randomized_overlap_0.05_eQTL)
enrichment_0.05_eQTL$category = "0.05"
colnames(enrichment_0.05_eQTL) = c("enrichment_ratio", "Category")
enrichment_0.01_eQTL = data.frame(eQTL_overlap_0.01 / randomized_overlap_0.01_eQTL)
enrichment_0.01_eQTL$category = "0.01"
colnames(enrichment_0.01_eQTL) =  c("enrichment_ratio", "Category")
enrichment_0.001_eQTL = data.frame(eQTL_overlap_0.001 / randomized_overlap_0.001_eQTL)
enrichment_0.001_eQTL$category = "0.001"
colnames(enrichment_0.001_eQTL) = c("enrichment_ratio", "Category")
enrichment_table_eQTL = rbind(enrichment_dhs_eQTL, enrichment_0.05_eQTL, enrichment_0.01_eQTL, enrichment_0.001_eQTL)
enrichment_table_eQTL$Category = factor(enrichment_table_eQTL$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table_eQTL, "/Users/ashvinravi/Desktop/bmi206/individual_project/enrichment_results/enrichment_table_eQTL_consensus.tsv")

boxplot_eQTL <- ggplot(enrichment_table_eQTL, aes(x=Category, y=enrichment_ratio, fill="eQTL")) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = "#ff5252") + 
  ylim(1, 3) +  ylab("Enrichment Ratio\n") + xlab("\nConsensus DHS Footprints")

############ sQTL enrichment 

list_of_sQTL_files = list.files("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_sQTL/", pattern = "sgenes.txt.gz")
setwd("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_sQTL/")
consensus_sQTL_variants = data.frame()
for (file in list_of_sQTL_files) {
  eGenes = read_tsv(file)
  print(file)
  significant_variants = eGenes %>% 
    filter(qval <= 0.05) %>% 
    select(chr, variant_pos, rs_id_dbSNP151_GRCh38p7) %>% 
    distinct()
  consensus_sQTL_variants = rbind(consensus_sQTL_variants, significant_variants)
}

consensus_sQTL_variants <- consensus_sQTL_variants %>% distinct()
colnames(consensus_sQTL_variants) = c("chrom", "start", "SNP")
consensus_sQTL_variants$end <- consensus_sQTL_variants$start + 1
### Join by minor allele frequency
sQTL_snps_maf <- inner_join(consensus_sQTL_variants, thousand_genome_panel_hg38, by = c("chrom", "start", "end", "SNP"))
dist.fun.consensus.sQTL = approxfun(density.default(sQTL_snps_maf$MAF))

sQTL_GR = toGRanges(data.frame(chr=sQTL_snps_maf$chrom, 
                               start=sQTL_snps_maf$start,
                               end=(sQTL_snps_maf$end)))

sQTL_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, sQTL_GR)
sQTL_overlap_0.05 = numOverlaps(consensus_0.05_GR, sQTL_GR)
sQTL_overlap_0.01 = numOverlaps(consensus_0.01_GR, sQTL_GR)
sQTL_overlap_0.001 = numOverlaps(consensus_0.001_GR, sQTL_GR)

randomized_overlap_0.05_sQTL = c()
randomized_overlap_0.01_sQTL = c()
randomized_overlap_0.001_sQTL = c()
randomized_overlap_dhs_sQTL = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random sQTL_snps from panel from 1000G panel. 
  randomized_sQTL_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 55805, weight_by = dist.fun.consensus.sQTL(MAF), replace=T) %>%
    ungroup()
  print(x)
  randomized_sQTL_snps_GR = toGRanges(data.frame(chr=randomized_sQTL_snps$chrom,
                                                 start=randomized_sQTL_snps$start,
                                                 end=as.numeric(randomized_sQTL_snps$end)))
  
  randomized_overlap_dhs_sQTL = c(randomized_overlap_dhs_sQTL, numOverlaps(randomized_sQTL_snps_GR, dhs_wout_footprints_GR))
  randomized_overlap_0.05_sQTL = c(randomized_overlap_0.05_sQTL, numOverlaps(randomized_sQTL_snps_GR, consensus_0.05_GR))
  randomized_overlap_0.01_sQTL = c(randomized_overlap_0.01_sQTL, numOverlaps(randomized_sQTL_snps_GR, consensus_0.01_GR))
  randomized_overlap_0.001_sQTL = c(randomized_overlap_0.001_sQTL, numOverlaps(randomized_sQTL_snps_GR, consensus_0.001_GR))
}

enrichment_dhs_sQTL = data.frame(sQTL_overlap_dhs / randomized_overlap_dhs_sQTL)
enrichment_dhs_sQTL$category = "DHS w/out footprints"
colnames(enrichment_dhs_sQTL) = c("enrichment_ratio", "Category")
enrichment_0.05_sQTL = data.frame(sQTL_overlap_0.05 / randomized_overlap_0.05_sQTL)
enrichment_0.05_sQTL$category = "0.05"
colnames(enrichment_0.05_sQTL) = c("enrichment_ratio", "Category")
enrichment_0.01_sQTL = data.frame(sQTL_overlap_0.01 / randomized_overlap_0.01_sQTL)
enrichment_0.01_sQTL$category = "0.01"
colnames(enrichment_0.01_sQTL) =  c("enrichment_ratio", "Category")
enrichment_0.001_sQTL = data.frame(sQTL_overlap_0.001 / randomized_overlap_0.001_sQTL)
enrichment_0.001_sQTL$category = "0.001"
colnames(enrichment_0.001_sQTL) = c("enrichment_ratio", "Category")
enrichment_table_sQTL = rbind(enrichment_dhs_sQTL, enrichment_0.05_sQTL, enrichment_0.01_sQTL, enrichment_0.001_sQTL)
enrichment_table_sQTL$Category = factor(enrichment_table_sQTL$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table_sQTL, "/Users/ashvinravi/Desktop/bmi206/individual_project/enrichment_results/enrichment_table_sQTL_consensus.tsv")

boxplot_sQTL <- ggplot(enrichment_table_sQTL, aes(x=Category, y=enrichment_ratio, fill="sQTL")) + geom_boxplot() + theme_classic() + scale_fill_manual(values = "#6497b1") + 
  ylim(1, 2.5) + ylab("Enrichment Ratio\n") + xlab("\nConsensus DHS Footprints")
