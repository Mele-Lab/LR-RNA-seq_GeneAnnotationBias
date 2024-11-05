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
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
# data <- fread("../novelannotations/analysis_tables/241031_exon_novelty_associated_gene.tsv", header=T)[, V1:=NULL]
data <- fread("../novelannotations/analysis_tables/241103_exons_novelty_pop_det.tsv")
sub <- data[, V1:=NULL]
sub[, `:=`(contig=tstrsplit(eid, "_")[[1]],
           start=tstrsplit(eid, "_")[[3]],
           end=tstrsplit(eid, "_")[[4]])][, eid:=NULL]
sub[, `:=`(newstart=as.numeric(start)-2, newend=as.numeric(end)+2)]
# fwrite(sub[, .(contig, newstart, newend)], "09_other_analyses/03_novel_exons_fst/data/01_expanded_novel&novel53prime_exons.v241031.bed", sep="\t", row.names = F, quote = F, col.names = F)
# load fst
fst <- fread("09_other_analyses/03_novel_exons_fst/data/CEUvsYRI_fst_noNA.fst", header=F)
colnames(fst) <- c("contig", "position", "FST")

library(GenomicRanges)
# Convert `sub` to a GRanges object
gr_sub <- GRanges(
  seqnames = sub$contig,
  ranges = IRanges(start = sub$newstart, end = sub$newend),
  mcols= sub[, .(novelty, AJI, CEU, HAC, ITU, LWK, MPC, PEL, YRI)]
)

# Convert `fst` to a GRanges object
gr_fst <- GRanges(
  seqnames = fst$contig,
  ranges = IRanges(start = fst$position, end = fst$position),
  mcols=fst[, .(FST)] # keep the 'fst' column as metadata
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
  fst = mcols(matched_fst)$mcols.FST,
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

result[fst<0, fst:=0]

# PLOT
result[, highfst:=fifelse(fst>=0.25, "fst>=0.25", "fst<0.25")]
ggplot(result, aes(x=fst, color=novelty))+
  geom_density()+
  facet_wrap(~highfst, scales="free_y")+
  mytheme
ggplot(result, aes(x=novelty, y=fst, fill=novelty))+
  geom_violin(alpha=0.7, trim = TRUE)+
  geom_boxplot(outliers = F, width=0.05)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  ggpubr::stat_compare_means(comparisons = list(c("Known", "Novel"), c("Known", "Novel 5'/3'")),method = "t.test",
                     method.args = list(alternative = "two.sided"))+
  labs(x="Exon", y="FST of Variants within exon+-2 bp")+
  stat_summary( fun=mean, geom = "point", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  guides(fill="none")+
  scale_fill_manual(values=c("darkgrey", "red","orange"))
ggplot(result[YRI>0], aes(x=novelty, y=fst, fill=novelty))+
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
p<-ggplot(result[YRI>0], aes(fst, col=novelty)) + 
  stat_ecdf(geom = "step")+
  mytheme+
  scale_color_manual(values=c("darkgrey", "red","orange"))+
  labs(x=expression(F[ST] ~ "(CEU-YRI)"), color="Exon")+
  ggmagnify::geom_magnify(from = c(xmin = 0.05, xmax = 0.2, ymin = 0.75, ymax = 1), 
               to = c(xmin = 0.3, xmax = 0.9, ymin = 0.1, ymax = 0.75))
ggsave(plot=p,"10_figures/suppfig/ECDF_fstCEUYRI.YRIdiscoveredExons.pdf", dpi=700, width = 18, height = 15,  units = "cm")

# zoomin <- ggplot(result[YRI>0], aes(fst, col=novelty)) + 
#   stat_ecdf(geom = "step")+
#   mytheme+
#   scale_color_manual(values=c("darkgrey", "red","orange"))+
#   labs(x=expression(F[ST]), color="Exon")+
#   coord_cartesian(xlim= c(0, 0.15), ylim=c(0.5, 1))
# main + 
#   ggmapinset::geom_inset(inset = zoomin, xmin = 0.25, xmax = 1, ymin = 0, ymax = 0.75)

ggplot(result, aes(x=novelty, y=fst, fill=novelty))+
  geom_boxplot(outliers = F)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  stat_summary( fun=mean, geom = "point", position = position_dodge(0.9),fun.args = list(y=-0.05))+
  ggpubr::stat_compare_means(comparisons = list(c("Known", "Novel"), c("Known", "Novel 5'/3'")),method = "t.test",
                             method.args = list(alternative = "two.sided"),   tip.length = 0.01, label.y = c(0.1, 0.12))+
  labs(x="Exon", y="FST of Variants within exon+-2 bp")+
  guides(fill="none")
