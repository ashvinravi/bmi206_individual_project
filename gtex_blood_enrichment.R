library(readr)
library(dplyr)
library(regioneR)
library(ggplot2)
library(tidyr)

# In this script, we will calculate the enrichment of lymphocyte DHS peaks/footprints in the lymphocyte count GWAS. 

thousand_genome_panel_hg38 = read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/thousand_genome_panel_hg38.tsv.gz")

gwas_catalog <- read_tsv("/Users/ashvinravi/Desktop/bmi206/group_project/gwas_catalog_v1.0-associations_e110_r2023-11-08.tsv")
european_gwas_catalog <- gwas_catalog %>% 
  filter(grepl("European", `INITIAL SAMPLE SIZE`)) %>% 
  filter(grepl("White blood", `DISEASE/TRAIT`)) %>% 
  select(CHR_ID, CHR_POS, SNPS) 

european_gwas_catalog$CHR_ID[is.na(european_gwas_catalog$CHR_ID)] = sub("(.*?)[\\.|:].*", "\\1", european_gwas_catalog$SNPS[is.na(european_gwas_catalog$CHR_ID)])
european_gwas_catalog$CHR_POS[is.na(european_gwas_catalog$CHR_POS)] = sub(".*:", "", european_gwas_catalog$SNPS[is.na(european_gwas_catalog$CHR_POS)]) 

european_gwas_catalog$CHR_ID[grepl("chr", european_gwas_catalog$CHR_ID)] <- gsub("^.{0,3}", "", european_gwas_catalog$CHR_ID[grepl("chr", european_gwas_catalog$CHR_ID)])
european_gwas_catalog <- european_gwas_catalog[!grepl("rs", european_gwas_catalog$CHR_ID),]  

gwas_snps <- european_gwas_catalog %>% distinct(CHR_ID, CHR_POS, .keep_all = T)
gwas_snps = gwas_snps[!grepl("x", gwas_snps$CHR_POS),]
gwas_snps = gwas_snps[!grepl(";", gwas_snps$CHR_POS),]
gwas_snps = gwas_snps[!is.na(gwas_snps$CHR_POS),]
colnames(gwas_snps) = c("chrom", "start", "SNP")

gwas_snps$start = as.numeric(gwas_snps$start)
gwas_snps$end = gwas_snps$start + 1 
gwas_snps$chrom = paste0("chr", gwas_snps$chrom)
gwas_snps_maf <- inner_join(gwas_snps, thousand_genome_panel_hg38, by = c("chrom", "start", "end", "SNP"))

dist.fun.gwas = approxfun(density.default(gwas_snps_maf$MAF))

######### Lymphocyte Samples 

lymphocyte_dhs_wout_footprints = read.table("/Users/ashvinravi/Desktop/bmi206/individual_project/lymphocyte_data/lymphocyte_dhs_peaks_wout_footprints.bed", sep='\t')
lymphocyte_0.001 = read.table("/Users/ashvinravi/Desktop/bmi206/individual_project/lymphocyte_data/lymphocyte_footprints_0.001.bed")
lymphocyte_0.01 = read.table("/Users/ashvinravi/Desktop/bmi206/individual_project/lymphocyte_data/lymphocyte_footprints_0.01.bed")
lymphocyte_0.05 = read.table("/Users/ashvinravi/Desktop/bmi206/individual_project/lymphocyte_data/lymphocyte_footprints_0.05.bed")

colnames(lymphocyte_dhs_wout_footprints) = c("chr", "start", "end")
colnames(lymphocyte_0.001) = c("chr", "start", "end")
colnames(lymphocyte_0.01) = c("chr", "start", "end")
colnames(lymphocyte_0.05) = c("chr", "start", "end")


dhs_wout_footprints_GR = toGRanges(data.frame(chr=lymphocyte_dhs_wout_footprints$chr,
                                              start=lymphocyte_dhs_wout_footprints$start, 
                                              end=lymphocyte_dhs_wout_footprints$end))

consensus_0.001_GR = toGRanges(data.frame(chr=lymphocyte_0.001$chr, 
                                          start=lymphocyte_0.001$start, 
                                          end=lymphocyte_0.001$end))
consensus_0.01_GR = toGRanges(data.frame(chr=lymphocyte_0.01$chr, 
                                         start=lymphocyte_0.01$start, 
                                         end=lymphocyte_0.01$end)) 
consensus_0.05_GR = toGRanges(data.frame(chr=lymphocyte_0.05$chr, 
                                         start=lymphocyte_0.05$start, 
                                         end=lymphocyte_0.05$end)) 

gwas_GR = toGRanges(data.frame(chr=gwas_snps_maf$chrom, 
                               start=gwas_snps_maf$start,
                               end=(gwas_snps_maf$end)))

gwas_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, gwas_GR)
gwas_overlap_0.05 = numOverlaps(consensus_0.05_GR, gwas_GR)
gwas_overlap_0.01 = numOverlaps(consensus_0.01_GR, gwas_GR)
gwas_overlap_0.001 = numOverlaps(consensus_0.001_GR, gwas_GR)

randomized_overlap_gwas_0.05 = c()
randomized_overlap_gwas_0.01 = c()
randomized_overlap_gwas_0.001 = c()
randomized_overlap_gwas_dhs = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random gwas_snps from panel from 1000G panel. 
  randomized_gwas_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 1600, 
                 weight_by = dist.fun.gwas(MAF),
                 replace=T) %>%
    ungroup()
  print(x)
  randomized_gwas_snps_GR = toGRanges(data.frame(chr=randomized_gwas_snps$chrom,
                                                 start=randomized_gwas_snps$start,
                                                 end=as.numeric(randomized_gwas_snps$end)))
  
  randomized_overlap_gwas_dhs = c(randomized_overlap_gwas_dhs, numOverlaps(randomized_gwas_snps_GR, dhs_wout_footprints_GR) + 1)
  randomized_overlap_gwas_0.05 = c(randomized_overlap_gwas_0.05, numOverlaps(randomized_gwas_snps_GR, consensus_0.05_GR) + 1)
  randomized_overlap_gwas_0.01 = c(randomized_overlap_gwas_0.01, numOverlaps(randomized_gwas_snps_GR, consensus_0.01_GR) + 1)
  randomized_overlap_gwas_0.001 = c(randomized_overlap_gwas_0.001, numOverlaps(randomized_gwas_snps_GR, consensus_0.001_GR) + 1)
}

enrichment_dhs_gwas = data.frame(gwas_overlap_dhs / randomized_overlap_gwas_dhs)
enrichment_dhs_gwas$category = "DHS w/out footprints"
colnames(enrichment_dhs_gwas) = c("enrichment_ratio", "Category")
enrichment_0.05_gwas = data.frame(gwas_overlap_0.05 / randomized_overlap_gwas_0.05)
enrichment_0.05_gwas$category = "0.05"
colnames(enrichment_0.05_gwas) = c("enrichment_ratio", "Category")
enrichment_0.01_gwas = data.frame(gwas_overlap_0.01 / randomized_overlap_gwas_0.01)
enrichment_0.01_gwas$category = "0.01"
colnames(enrichment_0.01_gwas) =  c("enrichment_ratio", "Category")
enrichment_0.001_gwas = data.frame(gwas_overlap_0.001 / randomized_overlap_gwas_0.001)
enrichment_0.001_gwas$category = "0.001"
colnames(enrichment_0.001_gwas) = c("enrichment_ratio", "Category")
enrichment_table_gwas = rbind(enrichment_dhs_gwas, enrichment_0.05_gwas, enrichment_0.01_gwas, enrichment_0.001_gwas)
enrichment_table_gwas$Category = factor(enrichment_table_gwas$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table_gwas, "/Users/ashvinravi/Desktop/bmi206/individual_project/blood_gwas_enrichment_table.tsv")

boxplot_gwas <- ggplot(enrichment_table_gwas, aes(x=Category, y=enrichment_ratio, fill=Category)) + geom_boxplot() + theme_classic()

####### eQTL enrichment 
gtex_blood_egenes <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")
significant_variant_pairs <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")

significant_variant_pairs <- significant_variant_pairs %>% separate(variant_id, into = c("chrom", "pos", "A1", "A2", "genome_build"), sep = "_")
significant_variants <- significant_variant_pairs %>% select(chrom, pos) %>% mutate(end = as.numeric(pos) + 1)

lead_blood_eQTL_variants <- gtex_blood_egenes %>% 
  filter(qval <= 0.05) %>% 
  select(chr, variant_pos, rs_id_dbSNP151_GRCh38p7) %>% 
  distinct()

colnames(lead_blood_eQTL_variants) = c("chrom", "start", "SNP")
lead_blood_eQTL_variants <- inner_join(lead_blood_eQTL_variants, thousand_genome_panel_hg38, by = c("chrom", "start", "SNP"))

dist.fun.eQTL <- approxfun(density.default(lead_blood_eQTL_variants$MAF))

eQTL_GR = toGRanges(data.frame(chr=lead_blood_eQTL_variants$chrom, 
                               start=lead_blood_eQTL_variants$start,
                               end=as.numeric(lead_blood_eQTL_variants$start) + 1))

eQTL_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, eQTL_GR)
eQTL_overlap_0.05 = numOverlaps(consensus_0.05_GR, eQTL_GR)
eQTL_overlap_0.01 = numOverlaps(consensus_0.01_GR, eQTL_GR)
eQTL_overlap_0.001 = numOverlaps(consensus_0.001_GR, eQTL_GR)

randomized_overlap_eQTL_dhs = c()
randomized_overlap_eQTL_0.05 = c()
randomized_overlap_eQTL_0.01 = c()
randomized_overlap_eQTL_0.001 = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random gwas_snps from panel from 1000G panel. 
  randomized_eQTL_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 9566, 
                 weight_by = dist.fun.eQTL(MAF),
                 replace=T) %>%
    ungroup()
  print(x)
  randomized_eQTL_snps_GR = toGRanges(data.frame(chr=randomized_eQTL_snps$chrom,
                                                 start=randomized_eQTL_snps$start,
                                                 end=as.numeric(randomized_eQTL_snps$end)))
  
  randomized_overlap_eQTL_dhs = c(randomized_overlap_eQTL_dhs, numOverlaps(randomized_eQTL_snps_GR, dhs_wout_footprints_GR) + 1)
  randomized_overlap_eQTL_0.05 = c(randomized_overlap_eQTL_0.05, numOverlaps(randomized_eQTL_snps_GR, consensus_0.05_GR) + 1)
  randomized_overlap_eQTL_0.01 = c(randomized_overlap_eQTL_0.01, numOverlaps(randomized_eQTL_snps_GR, consensus_0.01_GR) + 1)
  randomized_overlap_eQTL_0.001 = c(randomized_overlap_eQTL_0.001, numOverlaps(randomized_eQTL_snps_GR, consensus_0.001_GR) + 1)
}

enrichment_dhs_eQTL = data.frame(eQTL_overlap_dhs / randomized_overlap_eQTL_dhs)
enrichment_dhs_eQTL$category = "DHS w/out footprints"
colnames(enrichment_dhs_eQTL) = c("enrichment_ratio", "Category")
enrichment_0.05_eQTL = data.frame(eQTL_overlap_0.05 / randomized_overlap_eQTL_0.05)
enrichment_0.05_eQTL$category = "0.05"
colnames(enrichment_0.05_eQTL) = c("enrichment_ratio", "Category")
enrichment_0.01_eQTL = data.frame(eQTL_overlap_0.01 / randomized_overlap_eQTL_0.01)
enrichment_0.01_eQTL$category = "0.01"
colnames(enrichment_0.01_eQTL) =  c("enrichment_ratio", "Category")
enrichment_0.001_eQTL = data.frame(eQTL_overlap_0.001 / randomized_overlap_eQTL_0.001)
enrichment_0.001_eQTL$category = "0.001"
colnames(enrichment_0.001_eQTL) = c("enrichment_ratio", "Category")
enrichment_table_eQTL = rbind(enrichment_dhs_eQTL, enrichment_0.05_eQTL, enrichment_0.01_eQTL, enrichment_0.001_eQTL)
enrichment_table_eQTL$Category = factor(enrichment_table_eQTL$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table_eQTL, "/Users/ashvinravi/Desktop/bmi206/individual_project/blood_eQTL_enrichment_table.tsv")

boxplot_eQTL <- ggplot(enrichment_table_eQTL, aes(x=Category, y=enrichment_ratio, fill=Category)) + geom_boxplot() + theme_classic()

####### sQTL enrichment 
gtex_blood_sgenes <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_sQTL/Whole_Blood.v8.sgenes.txt.gz")
significant_splicing_variant_pairs <- read_tsv("/Users/ashvinravi/Desktop/bmi206/individual_project/GTEx_Analysis_v8_sQTL/Whole_Blood.v8.sqtl_signifpairs.txt.gz")

significant_splicing_variant_pairs <- significant_splicing_variant_pairs %>% separate(variant_id, into = c("chrom", "pos", "A1", "A2", "genome_build"), sep = "_")
significant_splicing_variants <- significant_splicing_variant_pairs %>% select(chrom, pos) %>% mutate(end = as.numeric(pos) + 1)

lead_blood_sQTL_variants <- gtex_blood_sgenes %>% 
  filter(qval <= 0.05) %>% 
  select(chr, variant_pos, rs_id_dbSNP151_GRCh38p7) %>% 
  distinct()

colnames(lead_blood_sQTL_variants) = c("chrom", "start", "SNP")
lead_blood_sQTL_variants <- inner_join(lead_blood_sQTL_variants, thousand_genome_panel_hg38, by = c("chrom", "start", "SNP"))

dist.fun.sQTL <- approxfun(density.default(lead_blood_sQTL_variants$MAF))

sQTL_GR = toGRanges(data.frame(chr=lead_blood_sQTL_variants$chrom, 
                               start=lead_blood_sQTL_variants$start,
                               end=as.numeric(lead_blood_sQTL_variants$start) + 1))

sQTL_overlap_dhs = numOverlaps(dhs_wout_footprints_GR, sQTL_GR)
sQTL_overlap_0.05 = numOverlaps(consensus_0.05_GR, sQTL_GR)
sQTL_overlap_0.01 = numOverlaps(consensus_0.01_GR, sQTL_GR)
sQTL_overlap_0.001 = numOverlaps(consensus_0.001_GR, sQTL_GR)

randomized_overlap_sQTL_0.05 = c()
randomized_overlap_sQTL_0.01 = c()
randomized_overlap_sQTL_0.001 = c()
randomized_overlap_sQTL_dhs = c()

print("Generating random regions...")
for (x in 1:1000) {
  # Sample Random gwas_snps from panel from 1000G panel. 
  randomized_sQTL_snps <- thousand_genome_panel_hg38 %>%
    slice_sample(n = 2846, 
                 weight_by = dist.fun.sQTL(MAF),
                 replace=T) %>%
    ungroup()
  print(x)
  randomized_sQTL_snps_GR = toGRanges(data.frame(chr=randomized_sQTL_snps$chrom,
                                                 start=randomized_sQTL_snps$start,
                                                 end=as.numeric(randomized_sQTL_snps$end)))
  
  randomized_overlap_sQTL_dhs = c(randomized_overlap_sQTL_dhs, numOverlaps(randomized_sQTL_snps_GR, dhs_wout_footprints_GR) + 1)
  randomized_overlap_sQTL_0.05 = c(randomized_overlap_sQTL_0.05, numOverlaps(randomized_sQTL_snps_GR, consensus_0.05_GR) + 1)
  randomized_overlap_sQTL_0.01 = c(randomized_overlap_sQTL_0.01, numOverlaps(randomized_sQTL_snps_GR, consensus_0.01_GR) + 1)
  randomized_overlap_sQTL_0.001 = c(randomized_overlap_sQTL_0.001, numOverlaps(randomized_sQTL_snps_GR, consensus_0.001_GR) + 1)
}

enrichment_dhs_sQTL = data.frame(sQTL_overlap_dhs / randomized_overlap_sQTL_dhs)
enrichment_dhs_sQTL$category = "DHS w/out footprints"
colnames(enrichment_dhs_sQTL) = c("enrichment_ratio", "Category")
enrichment_0.05_sQTL = data.frame(sQTL_overlap_0.05 / randomized_overlap_sQTL_0.05)
enrichment_0.05_sQTL$category = "0.05"
colnames(enrichment_0.05_sQTL) = c("enrichment_ratio", "Category")
enrichment_0.01_sQTL = data.frame(sQTL_overlap_0.01 / randomized_overlap_sQTL_0.01)
enrichment_0.01_sQTL$category = "0.01"
colnames(enrichment_0.01_sQTL) =  c("enrichment_ratio", "Category")
enrichment_0.001_sQTL = data.frame(sQTL_overlap_0.001 / randomized_overlap_sQTL_0.001)
enrichment_0.001_sQTL$category = "0.001"
colnames(enrichment_0.001_sQTL) = c("enrichment_ratio", "Category")
enrichment_table_sQTL = rbind(enrichment_dhs_sQTL, enrichment_0.05_sQTL, enrichment_0.01_sQTL, enrichment_0.001_sQTL)
enrichment_table_sQTL$Category = factor(enrichment_table_sQTL$Category, levels=c("DHS w/out footprints", "0.05", "0.01", "0.001"))
write_tsv(enrichment_table_sQTL, "/Users/ashvinravi/Desktop/bmi206/individual_project/blood_sQTL_enrichment_table.tsv")

boxplot_sQTL <- ggplot(enrichment_table_sQTL, aes(x=Category, y=enrichment_ratio, fill=Category)) + geom_boxplot() + theme_classic()



