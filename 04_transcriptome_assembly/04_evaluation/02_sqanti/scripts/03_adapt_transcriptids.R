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

sqanti <- fread("04_transcriptome_assembly/04_evaluation/02_sqanti/data/240909merge/240909merge_classification.txt")
gencode <- fread("/home/pclavell/mounts/mn5/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.all_genes.geneid.info.txt", header = F)
mytrx <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/novelannotations/merged/240910_transcripts_built_novel_genes.gtf")
sqanti <- sqanti[, .(isoform, structural_category, associated_gene, associated_transcript)]
colnames(gencode) <- c("geneid", "info")
gencode[, `:=`(gene_biotype=tstrsplit(info, "\"")[[4]], gene_symbol=tstrsplit(info, "\"")[[6]], info=NULL)]



sqanti[, one_associated_gene:=ifelse(structural_category=="intergenic", "novel",
                                     ifelse(structural_category=="antisense", tstrsplit(associated_gene, "_")[[2]],associated_gene))]
          
bigger <-gencode[sqanti, on=c("gene_symbol"="one_associated_gene")]
problematicfsm <-  unique(bigger[is.na(geneid)&structural_category!="intergenic" & structural_category=="full-splice_match",])
problematicgenic <-  unique(bigger[is.na(geneid)&structural_category!="intergenic" & structural_category=="genic",])

problematicfsm[, `:=`(gene1=tstrsplit(gene_symbol,"_")[[1]],
                      gene2=ifelse(grepl("_",gene_symbol), tstrsplit(gene_symbol,"_")[[2]], "singlegene"))]
problems1 <-gencode[, gene_bio1:=gene_biotype][, .(gene_bio1, gene_symbol)][problematicfsm, on=c("gene_symbol"="gene1")]
problems2 <-gencode[, gene_bio2:=gene_biotype][, .(gene_bio2, gene_symbol)][problems1, on=c("gene_symbol"="gene2")]

heatmap(table(problems2$gene_bio2, problems2$gene_bio1))
