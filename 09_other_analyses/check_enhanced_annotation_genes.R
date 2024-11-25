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


gencode <- fread("../../../Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf")
gencode <- gencode[V3!="exon"]
gencode[V3=="gene", geneid:=tstrsplit(V9,"\"")[[2]]]
enhanced <- fread("../novelannotations/241018_v47_poder_merge.gtf")
enhanced <- enhanced[V3!="exon"]
enhanced[V3=="gene", geneid:=tstrsplit(V9,"\"")[[2]]]
poder <- fread("../novelannotations/merged/poder_v1.gtf")
poder <- poder[V3!="exon"]
poder[V3=="gene", geneid:=tstrsplit(V9,"\"")[[2]]]

# are poder genes inside enhanced?
sum(unique(poder$geneid)%in%unique(enhanced$geneid))/length(unique(poder$geneid))


data <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
length(unique(data[structural_category%in%c("Intergenic")]$geneid.v))

