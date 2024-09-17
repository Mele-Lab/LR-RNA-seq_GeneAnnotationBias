# ## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:
##    write path relative to
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/03_recount3"
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
catch_args(1, MYPATH)

#MINSUPPORT=50
JUNCTIONS_FILE="../../04_evaluation/02_sqanti/data/pseudomasked_genomic_espresso/pseudomasked_genomic_espresso_junctions.txt"
##
ASSEMBLY<-gsub("_junctions.txt","",gsub(".*/", "",JUNCTIONS_FILE))

# ## 0----------------------------END OF HEADER----------------------------------0
args <- commandArgs(trailingOnly = T)
print("ARGS###################################################")
print(args)
MYPATH <- as.character(args[1])
cat(MYPATH)
library(recount3)
library(data.table)
#library(snapcount)


# load my junction file from sqanti
JUNCTIONS_FILE="../../04_evaluation/02_sqanti/data/pseudomasked_genomic_espresso/pseudomasked_genomic_espresso_junctions.txt"
myjunctions <- fread(JUNCTIONS_FILE)
myjunctions <-myjunctions[, .(isoform, chrom, strand, junction_number, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, canonical)]
myjunctions[, query := paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand)]
# load recount
recount <- fread("data/recount3_srv3h_morethan10counts.tsv")
colnames(recount) <- c("contig", "start", "end", "strand", "novel", "donor", "acceptor", "recount_samples", "recount_counts", "junction")



# PARSE BOTH DATAFRAMEs
myjunctions <-myjunctions[, .(isoform, chrom, strand, junction_number, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, canonical)]
myjunctions[, junction := paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand)]

mergedsj <-recount[, .(junction, recount_samples, recount_counts)][myjunctions, on="junction"]
mergedsj$recount_counts <- replace_na(mergedsj$recount_counts, 0)
mergedsj$recount_samples <- replace_na(mergedsj$recount_samples, 0)


# # plot recount_samples
# ggplot(mergedsj, aes(x=recount_samples, fill=junction_category)) +
#   geom_density(alpha=0.5)+
#   geom_vline(xintercept=25000, linetype="dashed", color="grey")+
#   mytheme
# 
# # plot recount_counts
# ggplot(mergedsj, aes(x=recount_counts+1, fill=junction_category)) +
#   geom_density(alpha=0.5)+
#   mytheme+
#   geom_vline(xintercept=500, linetype="dashed", color="grey")+
#   scale_x_continuous(
#     trans = "log10")

# # plot recount_samples and recount_counts
# ggplot(mergedsj, aes(x=log10(recount_samples+1), y=log10(recount_counts+1)) ) +
#   geom_hex(bins = 150) +
#   scale_fill_continuous(type = "viridis", trans = "log10") +
#   mytheme+
#   facet_wrap(~junction_category)
# 
# # plot ratio of recount_counts by recount_samples
# ggplot(mergedsj, aes(x=recount_counts/recount_samples, fill=junction_category)) +
#   geom_density(alpha=0.5)+
#   mytheme+
#   geom_vline(xintercept=1.6, linetype="dashed", color="grey")+
#   scale_x_continuous(
#     trans = "log10")

# decide on the thresholds
mergedsj[, recountsupported:=ifelse(recount_counts>=500  & recount_counts/recount_samples>=1.6, "recount_supp", "non_recount_supp")]

## COMPUTE PERCENTAGES
# nrow(mergedsj[recountsupported=="recount_supp" & junction_category=="novel"])/nrow(mergedsj[junction_category=="novel"])*100
# nrow(mergedsj[recountsupported=="non_recount_supp" & junction_category=="known"])/nrow(mergedsj[junction_category=="known"])*100

# # plot recount_samples and recount_counts
# ggplot(mergedsj, aes(x=log10(recount_samples+1), y=log10(recount_counts+1)) ) +
#   geom_hex(bins = 150) +
#   scale_fill_continuous(type = "viridis", trans = "log10") +
#   mytheme+
#   facet_wrap(~junction_category+recountsupported)

# Add column to state if we trust a novel SJ (unsupported known will be kept)
mergedsj[, reliability:=ifelse(junction_category=="novel" & recountsupported=="non_recount_supp", "unreliable", "reliable")]

# Convert SJ data.table to transcript-wise SJ info
mergedsj[, sj_reliability := paste(reliability, collapse=","), by=isoform]
mergedsj[, sj_support := paste(recount_counts, collapse=","), by=isoform]
mergedsj[, sj_supported := paste(recountsupported, collapse=","), by=isoform]
mergedsj[, sj_novelty := paste(junction_category, collapse=","), by=isoform]
mergedsj[, sj_start_site_novelty := paste(start_site_category, collapse=","), by=isoform]
mergedsj[, sj_end_site_novelty := paste(end_site_category, collapse=","), by=isoform]
mergedsj[, sj_canonical := paste(canonical, collapse=","), by=isoform]
mergedsj[, exons := sapply(strsplit(sj_canonical, ","), length)+1] # count SJ and add 1

# create the datatable to merge with the sqanti output
mylogical <- c(colnames(mergedsj)%in%c("isoform", "chrom", "strand") | grepl("sj_", colnames(mergedsj)))
sj_trxwise <- unique(mergedsj[,..mylogical])
nrow(sj_trxwise[!grepl("unreliable", sj_reliability),])/nrow(sj_trxwise)*100
nrow(sj_trxwise[!grepl("unreliable", sj_reliability) & grepl("novel", sj_novelty) ,])/nrow(sj_trxwise[grepl("novel", sj_novelty),])*100
