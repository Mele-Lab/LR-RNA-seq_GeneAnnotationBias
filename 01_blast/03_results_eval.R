## 0---------HEADER----------0
##
## AUTHOR: Pau Clavell-Revelles
## EMAIL: pauclavellrevelles@gmail.com
## CREATION DATE: 2023-12-20
##
## NOTES:
## 
##   
##
## ---------------------------

## DEFINE RELATIVE PATH (to projects)
relative_path <- "Projects/gencode_diversity/deduplication"

## DEFINE ARGUMENTS
args <- commandArgs(trailingOnly=TRUE)

if(length(args)>0){
  cat("CAPTURING ARGUMENTS...\n\n", sep= "")
  FILENAME <- args[1]
  NUMREADS <- args[2]
}


## set working directory for local and mn4
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")

machine <- ifelse(Sys.info()[7]=="pclavell", "local",
                  ifelse(Sys.info()[7]=="bsc83549", "mn4",
                         "ERROR"))
if(machine=="ERROR"){
  stop("ERROR: User not recognized, this is not your laptop nor your MN4 login")
}
cat(paste0("YOU ARE WORKING IN ", machine, "... \n\n"))

if(machine=="local"){
  projects_path_root <- "/home/pclavell/mounts/projects/"
  scratch_path <- "/home/pclavell/mounts/scratch/"
  home_path <- "/home/pclavell/mounts/home/"
}else if(machine=="mn4"){
  projects_path_root <- "/gpfs/projects/bsc83/"
  scratch_path <- "/gpfs/scratch/bsc83/"
  home_path <- "~/"
}

wd <- paste0(projects_path_root, relative_path, "/")
if(!dir.exists(wd)){dir.create(wd, recursive = T)}
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, "... \n\n",sep = "")
## ---------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 6, digits = 4) # non-scientific notation

## ---------------------------

## load up the packages we will need:
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
library(tidyverse)
library(data.table)
if(machine=="mn4"){setDTthreads(threads=48)}

## 0---------END OF HEADER-------------------------------------------------------------------0 ##
FILENAME <- "fiveLink_results.tsv"
NUMREADS <- 4438695
alignment_res <- fread(paste0("01_blast/01_results/",FILENAME))
colnames(alignment_res) <- c("qseqid", 
                             "qlen", "slen", 
                             "qstart", "qend", 
                             "sstart", "send", 
                             "qseq", "sseq", 
                             "evalue", "bitscore", 
                             "length", 
                             "pident", "nident", "mismatch", "positive", "gapopen", "gaps", 
                             "sstrand")

# add ranking of multimappings
alignment_res[, order_of_evalue := seq_len(.N), by = qseqid][,order_of_evalue:=factor(order_of_evalue)]
alignment_res[, umi := ifelse(grepl("NNNNNNNNNNNNNNN",sseq),"yes","no" )]
alignment_res[, num_mappings := .N, by=qseqid][, mapping:=ifelse(num_mappings>1, "multimapping", "unique")]

# plot number of mappings
uq <- unique(alignment_res, by="qseqid")
uq_mapcount <- table(uq[, "mapping"])
uq_mapcount<- setDT(as.data.frame(unlist(uq_mapcount)))
uq_mapcount <- rbind(uq_mapcount,data.table(mapping="unmapped", Freq=(NUMREADS-sum(uq_mapcount$Freq))), fill=T)
ggplot(uq_mapcount)+
  geom_col(aes(x=mapping, y=Freq, fill=mapping))+
  theme_bw()+
  scale_fill_manual(values=c("#ebb434", "#358a11","#8a1111"))+
  ylab("# reads")

# how many alignments contain UMI
ggplot(alignment_res)+
  geom_bar(aes(x=umi, fill=mapping), position="stack")+
  theme_bw()+
  ylab("# alignments")+
  scale_fill_manual(values=c("#ebb434", "#358a11"))

# how many reads have an alignment containing UMI
uu <- unique(alignment_res, by=c("qseqid", "umi"))
uu[, atleast1umi := ifelse(sum(umi=="yes")>0, "yes", "no"), by="qseqid"]

ggplot(unique(uu, by="qseqid"))+
  geom_bar(aes(x=atleast1umi, fill=atleast1umi))+
  theme_bw()+
  ylab("# reads")+
  xlab("At least 1 UMI?")

# check pident
ggplot(alignment_res)+
  geom_density(aes(x=pident, 
                   fill=umi),alpha=0.5)+
  facet_wrap(~umi)+
  theme_bw()

# check length
ggplot(alignment_res)+
  geom_density(aes(x=length, 
                   fill=umi),alpha=0.5)+
  theme_bw()

# check gaps
ggplot(alignment_res)+
  geom_bar(aes(x=gaps, fill=umi), position="dodge")+
  theme_bw()

# mismatch
ggplot(alignment_res[1:100000,])+
  geom_point(aes(x=nident, y=mismatch, color=umi, alpha=0.3))

ggplot(alignment_res)+
  geom_bar(aes(x=nident, fill=umi))
ggplot(alignment_res)+
  geom_bar(aes(x=mismatch, fill=umi))
