library(dplyr)
library(readr)

thousand_genome_panel <- read_tsv("Desktop/bmi206/group_project/1000G_Phase3_annotation_full.tsv.gz")

frequency_files <- list.files("/Users/ashvinravi/Desktop/bmi206/individual_project/1000G_Phase3_frq/")

frequency_panel = data.frame(matrix(ncol=6, nrow=0))
for (f in frequency_files) {
  print(f)
  frq = read.table(f, header = T)
  colnames(frq) = c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")
  colnames(frequency_panel) = c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")
  frequency_panel = rbind(frq, frequency_panel)
  
}

frequency_panel = as.data.frame(frequency_panel)
colnames(frequency_panel) = c("chrom", "SNP", "A1", "A2", "MAF", "NCHROBS")
frequency_panel$chrom = paste0("chr", frequency_panel$chrom)
total_panel = inner_join(thousand_genome_panel, frequency_panel, by=c("chrom", "SNP"))

write_tsv(total_panel, "../1000G_Phase3_panel_frq.tsv")
