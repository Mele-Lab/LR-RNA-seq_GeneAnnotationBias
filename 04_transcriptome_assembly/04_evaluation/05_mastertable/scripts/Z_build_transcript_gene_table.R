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
## 0----------------------------END OF HEADER----------------------------------
# load data
data <- fread("04_transcriptome_assembly/04_evaluation/02_sqanti/data/240909merge/240909merge_classification.txt")

fwrite(unique(data[,.(isoform, associated_gene)]), "04_transcriptome_assembly/04_evaluation/02_sqanti/data/240909merge/240909merge_associatedgene.tsv", 
       quote = F, row.names = F, col.names = F, sep="\t")


# prepare Enhanced annotation gene-transcript-biotype table
podermastertable <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
enhanced <- fread("../novelannotations/241018_v47_poder_merge.gtf")
poder <- fread("../novelannotations/merged/poder_v1.gtf")

gencode <- fread("../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
gencode <- gencode[, .(geneid.v, transcriptid.v, gene_biotype)]

#given these genes, add them the transcripts from gencode and the novel transcripts from poder
enhanced <- enhanced[V3=="transcript"]
enhanced[, `:=`(geneid.v=tstrsplit(V9, "\"")[[2]])]
enhanced <- enhanced[, .(geneid.v)]

enhanced <- gencode[enhanced, on="geneid.v"]
enhanced <- enhanced[!is.na(transcriptid.v)]

podersub <-unique(podermastertable[structural_category!="FSM", .(structural_category, geneid.v, isoform, associated_gene_biotype)][, `:=`(transcriptid.v=isoform, isoform=NULL, gene_biotype=associated_gene_biotype, associated_gene_biotype=NULL, structural_category=NULL)])
enhanced2 <- rbind.data.frame(enhanced, podersub)
enhanced2[gene_biotype=="protein_coding", gene_biotype:="Protein Coding"]
fwrite(enhanced2, "../novelannotations/241018_v47_poder_merge.gene_transcript_plusBiotypes.tsv", row.names = F, quote = F, sep="\t")
