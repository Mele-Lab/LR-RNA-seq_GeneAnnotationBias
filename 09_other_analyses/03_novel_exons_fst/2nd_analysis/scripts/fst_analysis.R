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

data <- data.frame()
paths <- list.files("09_other_analyses/03_novel_exons_fst/2nd_analysis/data/")
for(file in paths){
  predata <- fread(paste0("09_other_analyses/03_novel_exons_fst/2nd_analysis/data/", file))
  predata[, comp:=gsub(".FST.txt","",file)]
  data <- rbind.data.frame(data, predata)
}

colnames(data) <- c("contig", "pos", "fst", "pops")
data <- data[contig!="CHROM"]
data[fst=="-nan", fst:=NA]
data[, fst:=as.numeric(fst)]

data[, contains_ceu:=grepl("CEU", pops)]


data[contains_ceu==TRUE, meanfst:=mean(fst, na.rm=T), by=c("contig", "pos")]
fst <- data[contains_ceu==TRUE, .(meanfst, contig, pos)]
fst <- unique(fst[!is.na(meanfst)])

# poder <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
# poder[, ceufound:=fifelse(CEU>0 & population_sharing==1, "Only CEU",
#                           fifelse(CEU>0 & (CEU+ITU+LWK+ PEL+ YRI)>1, "Population Shared", 
#                                   fifelse(CEU==0 & (CEU+ITU+LWK+ PEL+ YRI)>1, "Only non-CEU", "Other")))]


sub <- fread("../novelannotations/analysis_tables/241103_exons_novelty_pop_det.tsv")
sub <- sub[, V1:=NULL]
sub[, `:=`(contig=tstrsplit(eid, "_")[[1]],
           start=as.numeric(tstrsplit(eid, "_")[[3]])+1, ## adding one because the table was generated from pyranges which converts to 0-based
           end=tstrsplit(eid, "_")[[4]])][, eid:=NULL]
sub[, `:=`(newstart=as.numeric(start)-2, newend=as.numeric(end)+2)]
### intersect

library(GenomicRanges)
# Convert `sub` to a GRanges object
gr_sub <- GRanges(
  seqnames = sub$contig,
  ranges = IRanges(start = sub$newstart, end = sub$newend),
  mcols= sub[, .(novelty, AJI, CEU, HAC, ITU, LWK, MPC, PEL, YRI)]
)

fst[, pos:=as.numeric(pos)]
# Convert `fst` to a GRanges object
gr_fst <- GRanges(
  seqnames = fst$contig,
  ranges = IRanges(start = fst$pos, end = fst$pos),
  mcols=fst[, .(meanfst)] # keep the 'fst' column as metadata
)


# Find overlaps between `gr_sub` and `gr_fst`
overlaps <- findOverlaps(gr_fst, gr_sub)

# Extract overlapping rows from `fst` and `sub`
matched_fst <- gr_fst[queryHits(overlaps)]
matched_sub <- gr_sub[subjectHits(overlaps)]

# Combine the results into a single data.table
result <- data.table(
  contig = as.character(seqnames(matched_fst)),
  position = start(matched_fst),
  fst = mcols(matched_fst)$mcols.meanfst,
  novelty = mcols(matched_sub)$mcols.novelty,
  AJI = mcols(matched_sub)$mcols.AJI,
  CEU = mcols(matched_sub)$mcols.CEU,
  HAC = mcols(matched_sub)$mcols.HAC,
  ITU = mcols(matched_sub)$mcols.ITU,
  LWK = mcols(matched_sub)$mcols.LWK,
  MPC = mcols(matched_sub)$mcols.MPC,
  PEL = mcols(matched_sub)$mcols.PEL,
  YRI = mcols(matched_sub)$mcols.YRI
)

#### FIND IF VARIANTS ARE IN SPLICE SITES-
# Calculate the exon boundaries (first 2 and last 2 bp)
gr_sub$first_2bp_start <- start(gr_sub)
gr_sub$first_2bp_end <- pmin(end(gr_sub), start(gr_sub) + 1) # First 2 bp
gr_sub$last_2bp_start <- pmax(start(gr_sub), end(gr_sub) - 1) # Last 2 bp
gr_sub$last_2bp_end <- end(gr_sub)

# Check if each position in `result` is within these boundaries
result[, within_2bp := {
  # Extract the matched ranges from `gr_sub`
  exon_start <- start(matched_sub)
  exon_end <- end(matched_sub)
  
  # Check if the position is within the first or last 2 bp
  is_within_first_2bp <- position >= exon_start & position <= (exon_start + 1)
  is_within_last_2bp <- position >= (exon_end - 1) & position <= exon_end
  
  # Combine the two conditions
  is_within_first_2bp | is_within_last_2bp
}]

rm(data, predata, fst)
gc()


# Final arrangements
result[fst<0, fst:=0]
result[, highfst:=fifelse(fst>=0.25, "fst>=0.25", "fst<0.25")]
colnames(result)[grep("within_2bp",colnames(result))] <- "splice_site_variant"


n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

ggplot(result[CEU==0 & AJI==0], aes(x=novelty, y=fst, fill=novelty))+
  geom_violin(alpha=0.7, trim = TRUE)+
  geom_boxplot(outliers = F, width=0.05)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  ggpubr::stat_compare_means(comparisons = list(c("Known", "Novel"), c("Known", "Novel 5'/3'"), c("Novel", "Novel 5'/3'")),method = "wilcox.test",
                             method.args = list(alternative = "two.sided"))+
  labs(x="Exon", y="FST of Variants within exon+-2 bp")+
  stat_summary( fun=mean, geom = "point", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  guides(fill="none")+
  scale_fill_manual(values=c("darkgrey", "red","orange"))

table(result[CEU==0 & AJI==0, .(novelty, highfst)])
