## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell"
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
data <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_data_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo.tsv")

## FILTER dataATION
gtf <- fread("../../noveldataations/merged/240917_merge_geneEntry_correctedScaffolds_nochrEBV.sorted.gtf")
gtf$V2 <- "Cerberus"
data[, filter:=ifelse((structural_category=="FSM" |
                          structural_category=="ISM" & existsFSMinTranscript==F |
                          structural_category=="NIC" & sample_sharing>1 & flair_max_counts>=3 & length >=300|
                          structural_category%in%c("NNC", "Genic", "Fusion", "Antisense", "Intergenic") & sample_sharing>1 & flair_max_counts>=3 & sj_less_recountsupported_counts>=25 & tool_sharing>1 & length >=300) &
                         associated_gene_biotype%in%unique(data$associated_gene_biotype)[grepl("lncRNA|^$|protein_coding", unique(data$associated_gene_biotype))], "pass", "fail")]

filtereddata <- data[filter=="pass" & structural_category!="ISM"]