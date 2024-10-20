## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
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

TYPE <- "pantrx"

# LOAD DATA
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  annot <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene_with_biotypes.tsv")
}

dge <- fread(paste0("data/02_DEGres_", TYPE,".tsv"))
dva <- fread(paste0("data/03_DVGres_", TYPE,".tsv"))
dtu <- fread(paste0("data/04_DTUres_", TYPE,".tsv"))



#### OVERLAP BETWEEN PAIRWISE COMPARISONS

for(analysis in c("dge", "dva", "dtu")){
  dgesig <- dge[fdr<0.05, .(contrast_name, geneid.v)][, sig:=1]
  dgesigwide <- dcast(dgesig, geneid.v~contrast_name, fill=0)
  UpSetR::upset(dgesigwide, nsets=28)
}


dgeupcon <- dge[fdr<0.05 & logFC>0.5, .(contrast, geneid.v)]
colnames(dgeupcon)[1] <- "pop"
dgeupred <- dge[fdr<0.05 & logFC< -0.5, .(reference, geneid.v)]
colnames(dgeupred)[1] <- "pop"
dgeup <- rbind.data.frame(dgeupcon, dgeupred)
# in How many contrasts does a single gene appear by ancestry

dgeup[, num_contrasts :=.N, by=c("pop", "geneid.v")]
ggplot(unique(dgeup), aes(x=num_contrasts, fill=pop))+
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"))+
  scale_fill_manual(values =popcols )+
  mytheme+
  labs(y="# Upregulated Genes", x="# Of Contrasts where a Gene is upregulated")



dgeup <- unique(dgeup)[, sig:=1]
dgeupwide <- dcast(dgeup, geneid.v~pop, fill=0)
UpSetR::upset(dgeupwide, nsets=8)


#### OVERLAP BETWEEN ANALYSES
colnames(dge)[colnames(dge%in%c("pval", "fdr", "logFC"))] <- paste0(c("pval", "fdr","logFC"), "_DGE")
colnames(dva)[colnames(dva%in%c("pval", "fdr", "LogVarRatio"))] <- paste0(c("pval", "fdr","LogVarRatio"), "_DVE")
colnames(dtu)[colnames(dtu%in%c("pval", "fdr"))] <- paste0(c("pval", "fdr"), "_DTU")
#### WHAT ARE THE BIOTYPES OF THE GENES


