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
library(data.table)
library(tidyverse)
data <- fread("~/mounts/mn5/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf", header=F)
data <-data[V3=="exon"]
data <- data[grepl("gene_type \"lncRNA", V9) | grepl("gene_type \"protein_coding", V9)]
data[, transcriptid.v :=tstrsplit(V9, "\"")[[4]]]
data[, exonlength:=abs(V5-V4)]
data[, length := sum(exonlength), by="transcriptid.v"]

gencode <- unique(data[, .(transcriptid.v, length)])
ggplot(gencode, aes(x=length))+
  geom_density()+
  xlim(c(0, 5000))+
  geom_vline(xintercept=300, linetype="dashed", color="darkgrey")+
  mytheme
nrow(gencode[length>=300])/nrow(gencode)
