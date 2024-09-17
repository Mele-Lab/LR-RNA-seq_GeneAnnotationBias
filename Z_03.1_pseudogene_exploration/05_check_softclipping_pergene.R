## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/03_mapping"
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

# load data
data <- fread("data_unmasked_filtered/1_PY1_GM10492_overlapped.bed", skip=18000000)


extract_gene_type <- function(gene_info) {
  matches <- regmatches(gene_info, regexpr('gene_type "([^"]+)"', gene_info))
  sub('gene_type "([^"]+)"', '\\1', matches)
}

data[, gene_biotype := sapply(V17, extract_gene_type)]

# Function to extract and sum the soft-clipped bases
extract_soft_clipped <- function(cigar_string) {
  matches <- gregexpr("([0-9]+)S", cigar_string)
  if (matches[[1]][1] == -1) {
    return(0)
  }
  soft_clips <- regmatches(cigar_string, matches)[[1]]
  soft_clip_numbers <- as.numeric(gsub("S", "", soft_clips))
  return(sum(soft_clip_numbers))
}

# Apply the function to column 7 and create a new column with the total number of soft-clipped bases
data[, softclip := sapply(V7, extract_soft_clipped)]

filt <- data[gene_biotype%in%c("protein_coding", "lncRNA", "processed_pseudogene")]
ggplot(filt, aes(x=softclip, fill=gene_biotype))+
  geom_density(alpha=0.5)+
  xlim(c(0,57))+mytheme
