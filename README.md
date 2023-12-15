# BMI 206: Individual Project
### Ashvin Ravi 
### December 15, 2023 

This repository contains all scripts for anlaysis + tissue-type specific source data for DHS footprinting data. 

## Raw Data: 

**Lung DHS data** can be found in ```lung_footprint_data/lung_union_data```. \
**Lymphocyte DHS data** can be found in ```lymphocyte_data```. \
**Consensus DHS data** is available in the ```dhs.rdata``` file. \

**Consensus DHS footprints** are available for download here (https://www.vierstra.org/resources/dgf#downloads). 

**GTEx v8 eQTL and sQTL** eGenes and sGenes are available for download here (https://www.gtexportal.org/home/downloads/adult-gtex/qtl). 

**NHGRI GWAS catalog** is available for download here (https://www.ebi.ac.uk/gwas/). 

## Analysis Scripts: 
```gtex_blood_enrichment.R``` and ```gtex_lung_enrichment.R``` corresponds to analysis scripts for Figure 6a and 6b respectively. 

```consensus_gwas_enrichment.R```corresponds to analysis scripts for Figure 6c.

```all_eQTL_sQTL.R```corresponds to analysis scripts for Figure 6d-6e.

The 1000G Phase 3 EUR panel is split by chromosome in the ```thousand_genome_panel_hg38_by_chr``` folder. These files can be concatenated to create a comprehensive panel. 

## Enrichment Results: 
All bootstrap test results are available in the ```enrichment_results``` folder. The ```create_boxplot.R``` provides the code for creating the figures. The produced figures are all available in ```enrichment_figures```. 





