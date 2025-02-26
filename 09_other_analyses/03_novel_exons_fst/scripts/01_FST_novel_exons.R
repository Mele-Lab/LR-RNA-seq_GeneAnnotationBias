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



data <- fread("../novelannotations/analysis_tables/250219_novel_exon_fsts.tsv")

novelpval <- format(wilcox.test(data[novelty=="Novel", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)
novelfivethree <- format(wilcox.test(data[novelty=="Novel 5'/3'", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)


svg("10_figures/01_plots/supp/10_fst/ECDF_fstCEUYRI.YRIdiscoveredExons.svg", width = 6, height = 4) 
ggplot(data, aes(mean_mean_fst, col=novelty)) + 
  stat_ecdf(geom = "step")+
  mytheme+
  scale_color_manual(values=c("darkgrey", "#AD5113","#519D8F"))+
  labs(x=expression("Mean " ~ F[ST] ~ "(CEU-All Populations)"), color="Exonic Region", y="Cumulative Proportion")+
  ggmagnify::geom_magnify(from = c(xmin = 0.025, xmax = 0.1, ymin = 0.65, ymax = 0.99), 
                          to = c(xmin = 0.15, xmax = 0.65, ymin = 0.1, ymax = 0.75))+
  annotate(geom="text", label=paste0("p=",novelpval), x=0.5, y=0.15, color="#AD5113")+
  annotate(geom="text", label=paste0("p=",novelfivethree), x=0.5, y=0.25, color="#519D8F")
dev.off()

# library(GenomicRanges)
# 
# 
# n_fun <- function(x, y){
#   return(data.frame(y = y, label = paste0("n = ",length(x))))
# }
# # data <- fread("../novelannotations/analysis_tables/241031_exon_novelty_associated_gene.tsv", header=T)[, V1:=NULL]
# data <- fread("../novelannotations/analysis_tables/241103_exons_novelty_pop_det.tsv")
# sub <- data[, V1:=NULL]
# sub[, `:=`(contig=tstrsplit(eid, "_")[[1]],
#            start=as.numeric(tstrsplit(eid, "_")[[3]])+1, ## adding one because the table was generated from pyranges which converts to 0-based
#            end=tstrsplit(eid, "_")[[4]])][, eid:=NULL]
# sub[, `:=`(newstart=as.numeric(start)-2, newend=as.numeric(end)+2)]
# # fwrite(sub[, .(contig, newstart, newend)], "09_other_analyses/03_novel_exons_fst/data/01_expanded_novel&novel53prime_exons.v241031.bed", sep="\t", row.names = F, quote = F, col.names = F)
# 
# # load fst
# fst <- fread("09_other_analyses/03_novel_exons_fst/data/CEUvsYRI_fst_noNA.fst", header=F)
# colnames(fst) <- c("contig", "position", "FST")
# 
# # Convert `sub` to a GRanges object
# gr_sub <- GRanges(
#   seqnames = sub$contig,
#   ranges = IRanges(start = sub$newstart, end = sub$newend),
#   mcols= sub[, .(novelty, AJI, CEU, HAC, ITU, LWK, MPC, PEL, YRI)]
# )
# 
# # Convert `fst` to a GRanges object
# gr_fst <- GRanges(
#   seqnames = fst$contig,
#   ranges = IRanges(start = fst$position, end = fst$position),
#   mcols=fst[, .(FST)] # keep the 'fst' column as metadata
# )
# 
# 
# # Find overlaps between `gr_sub` and `gr_fst`
# overlaps <- findOverlaps(gr_fst, gr_sub)
# 
# # Extract overlapping rows from `fst` and `sub`
# matched_fst <- gr_fst[queryHits(overlaps)]
# matched_sub <- gr_sub[subjectHits(overlaps)]
# 
# # Combine the results into a single data.table
# result <- data.table(
#   contig = as.character(seqnames(matched_fst)),
#   position = start(matched_fst),
#   fst = mcols(matched_fst)$mcols.FST,
#   novelty = mcols(matched_sub)$mcols.novelty,
#   AJI = mcols(matched_sub)$mcols.AJI,
#   CEU = mcols(matched_sub)$mcols.CEU,
#   HAC = mcols(matched_sub)$mcols.HAC,
#   ITU = mcols(matched_sub)$mcols.ITU,
#   LWK = mcols(matched_sub)$mcols.LWK,
#   MPC = mcols(matched_sub)$mcols.MPC,
#   PEL = mcols(matched_sub)$mcols.PEL,
#   YRI = mcols(matched_sub)$mcols.YRI
# )
# 
# #### FIND IF VARIANTS ARE IN SPLICE SITES-
# # Calculate the exon boundaries (first 2 and last 2 bp)
# gr_sub$first_2bp_start <- start(gr_sub)
# gr_sub$first_2bp_end <- pmin(end(gr_sub), start(gr_sub) + 1) # First 2 bp
# gr_sub$last_2bp_start <- pmax(start(gr_sub), end(gr_sub) - 1) # Last 2 bp
# gr_sub$last_2bp_end <- end(gr_sub)
# 
# # Check if each position in `result` is within these boundaries
# result[, within_2bp := {
#   # Extract the matched ranges from `gr_sub`
#   exon_start <- start(matched_sub)
#   exon_end <- end(matched_sub)
#   
#   # Check if the position is within the first or last 2 bp
#   is_within_first_2bp <- position >= exon_start & position <= (exon_start + 1)
#   is_within_last_2bp <- position >= (exon_end - 1) & position <= exon_end
#   
#   # Combine the two conditions
#   is_within_first_2bp | is_within_last_2bp
# }]
# 
# #3--
# 
# 
# 
# # Final arrangements
# result[fst<0, fst:=0]
# result[, highfst:=fifelse(fst>=0.25, "fst>=0.25", "fst<0.25")]
# colnames(result)[grep("within_2bp",colnames(result))] <- "splice_site_variant"
# 
# # Save data
# fwrite(result, "09_other_analyses/03_novel_exons_fst/data/variants_withinexonsOrSS.tsv", row.names = F, quote = F, sep="\t")
# 
# 
# # PLOT
# ggplot(result[YRI>0], aes(x=novelty, y=fst, fill=novelty))+
#   geom_violin(alpha=0.7, trim = TRUE)+
#   geom_boxplot(outliers = F, width=0.05)+
#   mytheme+
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-0.05))+
#   ggpubr::stat_compare_means(comparisons = list(c("Known", "Novel"), c("Known", "Novel 5'/3'"), c("Novel", "Novel 5'/3'")),method = "wilcox.test",
#                              method.args = list(alternative = "two.sided"))+
#   labs(x="Exon", y="FST of Variants within exon+-2 bp")+
#   stat_summary( fun=mean, geom = "point", position = position_dodge(0.9),fun.args = list(y=-0.05))+
#   guides(fill="none")+
#   scale_fill_manual(values=c("darkgrey", "red","orange"))
# novelpval <- 0.00022
# novelfivethree <- 0.0039
# p<-ggplot(result[YRI>0], aes(fst, col=novelty)) + 
#   stat_ecdf(geom = "step")+
#   mytheme+
#   scale_color_manual(values=c("darkgrey", "#AD5113","#519D8F"))+
#   labs(x=expression(F[ST] ~ "(CEU-YRI)"), color="Exon", y="Cumulative Proportion")+
#   ggmagnify::geom_magnify(from = c(xmin = 0.05, xmax = 0.2, ymin = 0.75, ymax = 1), 
#                to = c(xmin = 0.3, xmax = 0.9, ymin = 0.1, ymax = 0.75))+
#   annotate(geom="text", label=paste0("p=",novelpval), x=0.65, y=0.3, color="#AD5113")+
#   annotate(geom="text", label=paste0("p=",novelfivethree), x=0.65, y=0.37, color="#519D8F")
# 
# ggsave(plot=p,"10_figures/01_plots/supp/10_fst/ECDF_fstCEUYRI.YRIdiscoveredExons.pdf", dpi=700, width = 6, height = 4,  units = "in")
# 
# 


#####-------- REPEAT ANALYSIS TO INCLUDE INFO ABOUT SAMPLE DISCOVERY
# data <- fread("../novelannotations/analysis_tables/241126_exons_novelty_sample_det.tsv")
# sub <- data[, V1:=NULL]
# sub[, `:=`(contig=tstrsplit(eid, "_")[[1]],
#            start=as.numeric(tstrsplit(eid, "_")[[3]])+1, ## adding one because the table was generated from pyranges which converts to 0-based
#            end=tstrsplit(eid, "_")[[4]])][, eid:=NULL]
# sub[, `:=`(newstart=as.numeric(start)-2, newend=as.numeric(end)+2)]
# # fwrite(sub[, .(contig, newstart, newend)], "09_other_analyses/03_novel_exons_fst/data/01_expanded_novel&novel53prime_exons.v241031.bed", sep="\t", row.names = F, quote = F, col.names = F)
# 
# # load fst
# fst <- fread("09_other_analyses/03_novel_exons_fst/data/CEUvsYRI_fst_noNA.fst", header=F)
# colnames(fst) <- c("contig", "position", "FST")
# library(GenomicRanges)
# library(data.table)
# 
# # Convert `sub` to a GRanges object
# gr_sub <- GRanges(
#   seqnames = sub$contig,
#   ranges = IRanges(start = sub$newstart, end = sub$newend),
#   mcols = sub[, .(novelty, ITU5, ITU4, ITU3, ITU2, ITU1, PEL6, PEL5, PEL4, PEL3, PEL2, PEL1,
#                   HAC6, AJI6, AJI5, AJI4, AJI3, AJI2, AJI1, LWK5, LWK4, LWK3, LWK2, LWK1,
#                   YRI7, YRI6, YRI5, YRI3, HAC5, HAC4, HAC3, HAC2, HAC1, YRI2, YRI1,
#                   CEU5, CEU4, CEU3, CEU2, CEU1, MPC4, MPC3, MPC2, MPC1)]
# )
# 
# # Convert `fst` to a GRanges object
# gr_fst <- GRanges(
#   seqnames = fst$contig,
#   ranges = IRanges(start = fst$position, end = fst$position),
#   mcols = fst[, .(FST)] # Keep the 'fst' column as metadata
# )
# 
# # Find overlaps between `gr_sub` and `gr_fst`
# overlaps <- findOverlaps(gr_fst, gr_sub)
# 
# # Extract overlapping rows from `fst` and `sub`
# matched_fst <- gr_fst[queryHits(overlaps)]
# matched_sub <- gr_sub[subjectHits(overlaps)]
# 
# # Combine the results into a single data.table
# result <- data.table(
#   contig = as.character(seqnames(matched_fst)),
#   position = start(matched_fst),
#   fst = mcols(matched_fst)$mcols.FST,
#   novelty = mcols(matched_sub)$mcols.novelty,
#   ITU5 = mcols(matched_sub)$mcols.ITU5,
#   ITU4 = mcols(matched_sub)$mcols.ITU4,
#   ITU3 = mcols(matched_sub)$mcols.ITU3,
#   ITU2 = mcols(matched_sub)$mcols.ITU2,
#   ITU1 = mcols(matched_sub)$mcols.ITU1,
#   PEL6 = mcols(matched_sub)$mcols.PEL6,
#   PEL5 = mcols(matched_sub)$mcols.PEL5,
#   PEL4 = mcols(matched_sub)$mcols.PEL4,
#   PEL3 = mcols(matched_sub)$mcols.PEL3,
#   PEL2 = mcols(matched_sub)$mcols.PEL2,
#   PEL1 = mcols(matched_sub)$mcols.PEL1,
#   HAC6 = mcols(matched_sub)$mcols.HAC6,
#   AJI6 = mcols(matched_sub)$mcols.AJI6,
#   AJI5 = mcols(matched_sub)$mcols.AJI5,
#   AJI4 = mcols(matched_sub)$mcols.AJI4,
#   AJI3 = mcols(matched_sub)$mcols.AJI3,
#   AJI2 = mcols(matched_sub)$mcols.AJI2,
#   AJI1 = mcols(matched_sub)$mcols.AJI1,
#   LWK5 = mcols(matched_sub)$mcols.LWK5,
#   LWK4 = mcols(matched_sub)$mcols.LWK4,
#   LWK3 = mcols(matched_sub)$mcols.LWK3,
#   LWK2 = mcols(matched_sub)$mcols.LWK2,
#   LWK1 = mcols(matched_sub)$mcols.LWK1,
#   YRI7 = mcols(matched_sub)$mcols.YRI7,
#   YRI6 = mcols(matched_sub)$mcols.YRI6,
#   YRI5 = mcols(matched_sub)$mcols.YRI5,
#   YRI3 = mcols(matched_sub)$mcols.YRI3,
#   HAC5 = mcols(matched_sub)$mcols.HAC5,
#   HAC4 = mcols(matched_sub)$mcols.HAC4,
#   HAC3 = mcols(matched_sub)$mcols.HAC3,
#   HAC2 = mcols(matched_sub)$mcols.HAC2,
#   HAC1 = mcols(matched_sub)$mcols.HAC1,
#   YRI2 = mcols(matched_sub)$mcols.YRI2,
#   YRI1 = mcols(matched_sub)$mcols.YRI1,
#   CEU5 = mcols(matched_sub)$mcols.CEU5,
#   CEU4 = mcols(matched_sub)$mcols.CEU4,
#   CEU3 = mcols(matched_sub)$mcols.CEU3,
#   CEU2 = mcols(matched_sub)$mcols.CEU2,
#   CEU1 = mcols(matched_sub)$mcols.CEU1,
#   MPC4 = mcols(matched_sub)$mcols.MPC4,
#   MPC3 = mcols(matched_sub)$mcols.MPC3,
#   MPC2 = mcols(matched_sub)$mcols.MPC2,
#   MPC1 = mcols(matched_sub)$mcols.MPC1
# )
# 
# # Add splice site logic
# gr_sub$first_2bp_start <- start(gr_sub)
# gr_sub$first_2bp_end <- pmin(end(gr_sub), start(gr_sub) + 1)  # First 2 bp
# gr_sub$last_2bp_start <- pmax(start(gr_sub), end(gr_sub) - 1)  # Last 2 bp
# gr_sub$last_2bp_end <- end(gr_sub)
# 
# result[, splice_site_variant := {
#   exon_start <- start(matched_sub)
#   exon_end <- end(matched_sub)
#   within_first_2bp <- position >= exon_start & position <= (exon_start + 1)
#   within_last_2bp <- position >= (exon_end - 1) & position <= exon_end
#   within_first_2bp | within_last_2bp
# }]
# 
# # Final adjustments
# result[fst < 0, fst := 0]
# result[, highfst := fifelse(fst >= 0.25, "fst>=0.25", "fst<0.25")]
# colnames(result)[grep("splice_site_variant", colnames(result))] <- "splice_site_variant"
# 
# fwrite(result, "09_other_analyses/03_novel_exons_fst/data/variants_withinexonsOrSS_sampleinfo_included.tsv", row.names = F, quote = F, sep="\t")
# 
