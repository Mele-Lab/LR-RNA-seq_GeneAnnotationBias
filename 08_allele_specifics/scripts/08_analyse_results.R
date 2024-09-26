## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/08_allele_specifics"
##
##    load setup sources
mn5source <- "/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
localsource <-"/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
ifelse(file.exists(mn5source), source(mn5source), source(localsource))
rm(mn5source, localsource)
##    sets wd, load tidyverse, data.table, sets DT threads
setup_script(relative_path, 3, 48)
##
##    catch arguments from standard input 
##      catch_args(number_of_args, "obj_name1", "obj_name2")
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------0
data <- fread("data/04_calc_ase/41_EU1_GM12273_ase_annotated.tsv")

data <- fread("data/04_calc_ase/10_NI5_GM19117_ase_annotated.tsv")
data$fdr <- p.adjust(data$BINOM_P, method="BH")
data[,geneid.v:=tstrsplit(GENE_ID,",")[[1]]]

length(unique(data[fdr<0.05, geneid.v]))
length(unique(data[GENOTYPE_WARNING==0 & BLACKLIST==0 & MULTI_MAPPING==0 & OTHER_ALLELE_WARNING==0 & HIGH_INDEL_WARNING==0 & rawDepth>=20, geneid.v]))
nrow(data[GENOTYPE_WARNING==0 & BLACKLIST==0 & MULTI_MAPPING==0 & OTHER_ALLELE_WARNING==0 & HIGH_INDEL_WARNING==0 & rawDepth>=20, ])
hist(data[fdr<0.05, BINOM_P])
