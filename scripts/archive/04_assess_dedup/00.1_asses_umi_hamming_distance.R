## 0---------HEADER----------0
##
## AUTHOR: Pau Clavell-Revelles
## EMAIL: pauclavellrevelles@gmail.com
## CREATION DATE: 2024-01-18
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
  SEP <- args[2]
  WIN <- args[3]
  MIN <- args[4]
  CHR <- args[5]
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

# load data
data <- fread("02_extract_UMI/01_extracted_UMI/FAX06702_extracted_UMI.tsv", nrows=10000)
data <- fread("random_sample.txt", nrows=10000)


# Convert each string into a list of characters
char_list <- strsplit(unlist(data[,4]), "")

# Find the maximum length of characters in any string
min(sapply(char_list, length))
max(sapply(char_list, length))

# IF THE LENGTH IS THE SAME:
result_dataframe <- data.frame(matrix(unlist(char_list), ncol = length(char_list), byrow = TRUE))

test2<-cultevo::hammingdists(t(result_dataframe))

hist(test2)
