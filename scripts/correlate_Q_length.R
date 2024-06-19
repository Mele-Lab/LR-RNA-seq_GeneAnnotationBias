## 0---------HEADER----------0
##
## AUTHOR: Pau Clavell-Revelles
## EMAIL: pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-14
##
## NOTES:
## 
##   
##
## ---------------------------

## DEFINE RELATIVE PATH (to projects/bsc83)
relative_path <- "Projects/pantranscriptome/pclavell/ONT_preprocessing"

## DEFINE ARGUMENTS
args <- commandArgs(trailingOnly=TRUE)

if(length(args)>0){
  cat("CAPTURING ARGUMENTS...\n\n", sep= "")
  FILENAME <- args[1]
  SAMPLE <- args[2]
}


## set working directory for local and mn5
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")

machine <- ifelse(Sys.info()[7]=="pclavell", "local",
                  ifelse(grepl("bsc", Sys.info()[7]), "mn5",
                         "ERROR"))
if(machine=="ERROR"){
  stop("ERROR: User not recognized, this is not your laptop nor your MN5 login")
}
cat(paste0("YOU ARE WORKING IN ", machine, "... \n\n"))

if(machine=="local"){
  mn5projects <- "/home/pclavell/mounts/mn5/"
}else if(machine=="mn5"){
  mn5projects <- "/gpfs/projects/bsc83/"
}

wd <- paste0(mn5projects, relative_path, "/")
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
if(machine=="mn5"){setDTthreads(threads=48)}

## 0---------END OF HEADER-------------------------------------------------------------------0 ##
# nanostats <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/nanoplotstats/20240212_HS_30_PE2_HG01928_preprocessed.fastq.gz")

nanostats <- fread(FILENAME)
correlation <- cor(nanostats$quals, nanostats$lengths)
cat("Cor_Q_length", correlation)
sorted <- nanostats$lengths[order(nanostats$length, decreasing=T)]
cumsumm <- cumsum(as.numeric(sorted))
sum <- sum(sorted)
n50<-sorted[min(which(cumsumm>(sum/2)))]
cat("N50", n50)

thousandbp <-sum(nanostats$lengths>1000)/nrow(nanostats)*100


fwrite(data.frame("corr_Qlength"=correlation, "n50"=n50, "reads_longer1k"=thousandbp), paste0("nanoplotstats/", SAMPLE, "_corrQlength_n50.tsv"))
