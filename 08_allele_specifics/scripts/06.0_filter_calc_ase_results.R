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
catch_args(3, "SAMPLE", "CELLINE", "TYPE")
##
## 0----------------------------END OF HEADER----------------------------------0
data <- fread(paste0("data/",TYPE,"/04_calc_ase/", SAMPLE, "_ase_annotated.tsv"))
data[, filter:=ifelse(GENOTYPE_WARNING==0 & BLACKLIST==0 & MULTI_MAPPING==0 & OTHER_ALLELE_WARNING==0 & HIGH_INDEL_WARNING==0 & rawDepth>=20, "pass", "fail")]
filtered_data <- data[filter=="pass",][,BINOM_P_ADJUSTED:=p.adjust(BINOM_P, method="BH")][, filter:=NULL]
fwrite(filtered_data, paste0("data/",TYPE,"/04_calc_ase/", SAMPLE,"_ase_annotated_filtered.tsv"), sep="\t", quote = F, row.names = F)
