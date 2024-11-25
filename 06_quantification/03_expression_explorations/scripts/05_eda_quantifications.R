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
library(edgeR)
library(ggbeeswarm)
library(ggh4x)
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
library(ggpubr)

gcounts <- fread("../novelannotations/quantifications/kallisto_quant/matrix.abundance.genelevel.tsv")
tcounts <- fread("../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")

# Load data
data <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
metadata[, total_pop_throughput:=sum(map_reads_generalmap), by="population"]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,4,2,6,7,3,5)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
data[, associated_gene_biotype := factor(associated_gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene"))]
data[associated_gene_biotype=="Novel/Ambiguous Gene", associated_gene_biotype:="Novel/Ambiguous\nGene"]
# Compute cpm of gene counts
# gcounts <- column_to_rownames(gcounts, var="geneid.v")
# gcounts <- as.matrix(gcounts)
# gcpm <- cpm(gcounts)
# 
# # Transform to long format data.table
# gcpm <- gcpm[rowSums(gcpm>=0.1)>0,]
# gcpm <- rownames_to_column(as.data.frame(gcpm), var="geneid.v")
gcounts_long <- melt(gcounts, id.vars="geneid.v", variable.name="sample", value.name = "count")
setDT(gcounts_long)
gcounts_long <- unique(data[, .(associated_gene_biotype, geneid.v)])[gcounts_long, on=c("geneid.v")]
gcounts_long <- metadata[, .(sample, population, map_reads_generalmap,total_pop_throughput)][gcounts_long, on="sample"]
gcounts_long[, Expressed:=fifelse(count>=1, "Expressed", "Not Expressed")]
gcounts_long[is.na(associated_gene_biotype), associated_gene_biotype:="Novel/Ambiguous Gene"]
gcounts_long[, count_expressed:=uniqueN(geneid.v), by=c("sample", "Expressed", "associated_gene_biotype")]

ggplot(unique(gcounts_long[, .(population,sample, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed,associated_gene_biotype)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed, col=Expressed))+
  geom_quasirandom(width = 0.3, aes(size=map_reads_generalmap/10^6, alpha=Expressed))+
  mytheme+
  ylim(c(0,max(gcounts_long$count_expressed)))+
  labs(x="", y="# Expressed PODER genes", color="", size="Reads (M)", alpha="")+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  scale_alpha_manual(values=c(0.6, 0.2))+
  scale_size_continuous( range = c(0.5, 3))+
  theme(
    legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
    legend.margin = margin(0, 0, 0, 0),
    legend.position ="top"
  )+
  guides(
    col = guide_legend(ncol = 2, byrow = TRUE),  # Two columns for the color legend
    size = guide_legend(ncol = 3)  # Single column for the size legend
  )+
  facet_wrap(~associated_gene_biotype)+
  stat_summary(data = unique(gcounts_long[Expressed=="Expressed", .(population,sample, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed,associated_gene_biotype)]), 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed,group = associated_gene_biotype),
               fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 1)
ggsave("10_figures/suppfig_quantification2/jitter_ExpressedGenes_PerPopulation.pdf", dpi=700, width = 6.5, height = 2.5,  units = "in")



#### Do the same with trx
# # Compute cpm of gene counts
# tcounts <- column_to_rownames(tcounts, var="transcript_id")
colnames(tcounts) <- gsub("_1","", colnames(tcounts))
# tcounts <- as.matrix(tcounts)
# tcpm <- cpm(tcounts)

# Transform to long format data.table
# tcpm <- tcpm[rowSums(tcpm>0)>0,]
# tcpm <- rownames_to_column(as.data.frame(tcpm), var="transcriptid.v")
tcounts_long <- melt(tcounts, id.vars="transcript_id", variable.name="cell_line_id", value.name = "count")
setDT(tcounts_long)
tcounts_long <- unique(data[, .(associated_gene_biotype, geneid.v, isoform,structural_category)])[tcounts_long, on=c("isoform"="transcript_id")]
tcounts_long <- metadata[, .(sample, cell_line_id, population, map_reads_generalmap,total_pop_throughput)][tcounts_long, on="cell_line_id"]
tcounts_long[, Expressed:=fifelse(count>=1000, "Expressed", "Not Expressed")]
tcounts_long[, count_expressed:=uniqueN(isoform), by=c("sample", "Expressed")]
tcounts_long[, count_expressed_cat_biotype:=uniqueN(isoform), by=c("sample", "Expressed", "associated_gene_biotype", "structural_category")]
tcounts_long[, count_expressed_cat:=uniqueN(isoform), by=c("sample", "Expressed", "structural_category")]
tcounts_long[, eur:=fifelse(population%in%c("CEU", "AJI"), "European", "non-European")]
tcounts_long[, afr:=fifelse(population%in%c("YRI", "MPC", "LWK"), "African", "OOA")]
tcounts_long[, expressed_transcripts_per_gene_persample:=uniqueN(isoform)/map_reads_generalmap*10^6, by=c("sample", "geneid.v")]


ggplot(unique(tcounts_long[, .(population, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed, col=Expressed))+
  stat_summary(data = tcounts_long[Expressed == "Expressed"], 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed),
               fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
  geom_quasirandom(width = 0.3, aes(size=map_reads_generalmap/10^6, alpha=Expressed))+
  mytheme+
  ylim(c(0,max(tcounts_long$count_expressed)))+
  labs(x="", y="# Expressed PODER transcripts", color="", size="Reads (M)", alpha="")+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  scale_alpha_manual(values=c(0.6, 0.2))+
  theme(
    legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
    legend.margin = margin(0, 0, 0, 0),
    legend.position ="top"
  )+
  scale_size_continuous( range = c(0.5, 3))+
  guides(
    col = guide_legend(ncol = 1, byrow = TRUE),  # Two columns for the color legend
    size = guide_legend(ncol = 3)  # Single column for the size legend
  )
ggsave("10_figures/suppfig_quantification2/jitter_ExpressedTrx_PerPopulation.pdf", dpi=700, width = 3.25, height = 2.5,  units = "in")


ggplot(unique(tcounts_long[, .(population, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed_cat,structural_category)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed_cat, col=Expressed))+
  geom_quasirandom(width = 0.3, aes(size=map_reads_generalmap/10^6, alpha=Expressed))+
  mytheme+
  labs(x="", y="# Expressed PODER transcripts", color="", size="Reads (M)", alpha="")+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  scale_alpha_manual(values=c(0.6, 0.2))+
  scale_size_continuous( range = c(0.5, 3))+
  theme(
    legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
    legend.margin = margin(0, 0, 0, 0),
    legend.position ="top"
  )+
  guides(
    col = guide_legend(ncol = 2, byrow = TRUE),  # Two columns for the color legend
    size = guide_legend(ncol = 3)  # Single column for the size legend
  )+
  facet_wrap(~structural_category, scales="free_y")+
  stat_summary(data = unique(tcounts_long[Expressed=="Expressed", .(population, total_pop_throughput,map_reads_generalmap,structural_category, Expressed, count_expressed_cat)]), 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed_cat, group=structural_category),
               fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 1)
ggsave("10_figures/suppfig_quantification2/jitter_ExpressedTrx_PerPopulationPerSQANTI.pdf", dpi=700, width = 6, height = 6,  units = "in")





# PLOT DETECTION VS DEPTH
discoverylines <-unique(tcpm_long[Expressed == "Expressed", .(sample, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)])
discoverylines <-unique(discoverylines[, totalcount :=sum(count_expressed_cat), by="sample"][, .(sample, population,map_reads_generalmap,totalcount)])
ggplot(discoverylines, 
       aes(x=map_reads_generalmap/10^6, y=totalcount))+
  geom_smooth(method = "lm", color="darkgrey")+
  geom_point(aes( col=population))+
  mytheme+
  labs( y = "# Expressed PODER transcripts",  x = "Reads (M)", color="Population")+
  scale_color_manual(values = popcol) +
  stat_cor(label.y = 75000, size=6*0.35) +
  stat_regline_equation(label.y = 72500, size=6*0.35)+
  theme(
    legend.key.size = unit(0.05, "cm"),  # Reduce size of legend keys
    legend.margin = margin(0, 0, 0, 0)  # Remove extra margin around legend
  )
ggsave("10_figures/suppfig_quantification1/jitter_ExpressedTrx_PerThroughput.pdf", dpi=700, width = 2.75, height = 2,  units = "in")


n_fun <- function(x) {
  min_y <- min(x) - 0.1 * abs(min(x))  # Place label just below the minimum value with a small margin
  return(data.frame(y = min_y, label = paste0("n = ", length(x))))
}
p <-ggplot(unique(tcounts_long[Expressed == "Expressed", .(eur, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)]), 
           aes(x = eur, y = count_expressed_cat/map_reads_generalmap*10^6, fill = eur)) +
  geom_violin(alpha = 0.75) +
  geom_quasirandom(width = 0.3, alpha = 0.6, aes(size = map_reads_generalmap / 10^6, color = population)) +
  geom_boxplot(outliers = FALSE, width = 0.05, show.legend = F) +
  mytheme +
  labs(x = "", y = "# Expressed PODER transcripts / Million Reads", color = "Population", size = "Reads (M)", fill = "") +
  facet_wrap(~structural_category, scales = "free_y") +
  scale_color_manual(values = popcol) +
  scale_fill_manual(values = c("#466995", "#A53860"),
                    guide = guide_legend(override.aes = list(size = 3) )) +
  geom_pwc(ref.group = "European", method = "t_test", label.size=6*0.35) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9), size=6*0.35) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  scale_size_continuous( range = c(0.5, 3))+
  theme(legend.position = c(0.65, 0.1),legend.box = "horizontal",
        legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
        legend.margin = margin(0, 0, 0, 0))+
  guides(fill = guide_legend(override.aes = list(shape = NA, color=NA),
                             keywidth = unit(0.2, "cm"),          # Adjust width of the legend squares
                             keyheight = unit(0.2, "cm")))  # Remove shape and fill from legend

ggadjust_pvalue(p)
ggsave("10_figures/suppfig_quantification3/jitter_ExpressedTrx_PerSQANTI&EURnonEUR.pdf", dpi=700, width = 6, height = 6,  units = "in")

p <-ggplot(unique(tcounts_long[Expressed == "Expressed", .(afr, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)]), 
           aes(x = afr, y = count_expressed_cat/map_reads_generalmap*10^6, fill = afr)) +
  geom_violin(alpha = 0.75) +
  geom_quasirandom(width = 0.3, alpha = 0.6, aes(size = map_reads_generalmap / 10^6, color = population)) +
  geom_boxplot(outliers = FALSE, width = 0.05, show.legend = F) +
  mytheme +
  labs(x = "", y = "# Expressed PODER transcripts / Million Reads", color = "Population", size = "Reads (M)", fill = "") +
  facet_wrap(~structural_category, scales = "free_y") +
  scale_color_manual(values = popcol) +
  scale_fill_manual(values = c("#466995", "#A53860"),
                    guide = guide_legend(override.aes = list(size = 3) )) +
  geom_pwc(ref.group = "African", method = "t_test", label.size=6*0.35) +
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9), size=6*0.35) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  scale_size_continuous( range = c(0.5, 3))+
  theme(legend.position = c(0.65, 0.1),legend.box = "horizontal",
        legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
        legend.margin = margin(0, 0, 0, 0))+
  guides(fill = guide_legend(override.aes = list(shape = NA, color=NA),
                             keywidth = unit(0.2, "cm"),          # Adjust width of the legend squares
                             keyheight = unit(0.2, "cm")))  # Remove shape and fill from legend

ggadjust_pvalue(p)

# PLOTS ABOUT TRANSCCRIPTS X GENE
ggplot(unique(tcounts_long[, .(population, sample, total_pop_throughput, afr, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
              aes(x=expressed_transcripts_per_gene_persample, color=afr))+
  stat_ecdf(geom="step")+
  scale_color_manual(values=c("#F7D257", "#496F5D"))+
  mytheme+
  ggmagnify::geom_magnify(from = c(xmin = 0.05, xmax = 5, ymin = 0.75, ymax = 1), 
                          to = c(xmin =10, xmax = 37, ymin = 0.1, ymax = 0.75))
ggplot(unique(tcounts_long[, .(population, sample, total_pop_throughput, afr, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
       aes(x=afr, y=expressed_transcripts_per_gene_persample, fill=afr))+
  geom_violin(alpha=0.7, adjust=2 )+
  geom_boxplot(width=0.1, outliers = F)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=0))+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "wilcox.test",
                   method.args = list(alternative = "two.sided"))+
  scale_y_continuous(trans="log10")+
  labs(x="", y="# Expressed Transcripts/Sample\nnormalized by M reads")+
  annotation_logticks(sides = "l")+
  guides(fill="none")
# ggplot(unique(tcounts_long[, .(population, sample, total_pop_throughput, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
#        aes(x=eur, y=expressed_transcripts_per_gene_persample, fill=eur))+
#   geom_violin(alpha=0.7, adjust=2 )+
#   geom_boxplot(width=0.1, outliers = F)+
#   scale_fill_manual(values = c("#466995", "#A53860")) +
#   mytheme+
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-2))+
#   stat_compare_means(comparisons = list(c("European", "non-European")),method = "wilcox.test",
#                      method.args = list(alternative = "two.sided"))+
#   scale_y_continuous(trans="log10")+
#   labs(x="", y="# Expressed Transcripts/Sample normalized by M reads")+
#   annotation_logticks(sides = "l")+
#   guides(fill="none")


# 
# 
# ggplot(unique(tcpm_long[Expressed == "Expressed", .(eur, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)]), 
#        aes(x = eur, y = count_expressed_cat/map_reads_generalmap*10^6, fill = eur)) +
#   geom_violin(alpha = 0.75) +
#   geom_boxplot(outliers = FALSE, width = 0.05) +
#   geom_quasirandom(width = 0.3, alpha = 0.6, aes(size = map_reads_generalmap / 10^6, color = population)) +
#   mytheme +
#   labs(x = "", y = "# Expressed PODER transcripts / Million Reads", color = "Population", size = "Reads (M)", fill = "") +
#   facet_wrap(~structural_category, scales = "free_y") +
#   scale_color_manual(values = popcol) +
#   scale_fill_manual(values = c("#466995", "#A53860")) +
#   geom_pwc(ref.group = "European", method = "t_test") +
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9), size=6*0.35) +
#   guides(fill = guide_legend(override.aes = list(shape = NA)))+
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
#   theme(legend.position = c(0.65, 0.1),legend.box = "horizontal")
# ggsave("10_figures/suppfig_09/jitter_ExpressedTrx_PerBiotype&EURnonEUR.pdf", dpi=700, width = 4, height = 4,  units = "in")
# ggplot(unique(tcpm_long[Expressed == "Expressed", .(afr, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)]), 
#        aes(x = afr, y = count_expressed_cat/map_reads_generalmap*10^6, fill = afr)) +
#   geom_violin(alpha = 0.75) +
#   geom_boxplot(outliers = FALSE, width = 0.05) +
#   geom_quasirandom(width = 0.3, alpha = 0.6, aes(size = map_reads_generalmap / 10^6, color = population)) +
#   mytheme +
#   labs(x = "", y = "# Expressed PODER transcripts / Million Reads", color = "Population", size = "Reads (M)", fill = "") +
#   facet_wrap(~structural_category, scales = "free_y") +
#   scale_color_manual(values = popcol) +
#   scale_fill_manual(values =c("#F7D257", "#496F5D")) +
#   geom_pwc(ref.group = "African", method = "t_test",method.args = list(alternative="greater")) +
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9), size=6*0.35) +
#   guides(fill = guide_legend(override.aes = list(shape = NA)))+
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
#   theme(legend.position = c(0.65, 0.1),legend.box = "horizontal")
# 
# # check at distribution of tpm per category
# ggplot(tcpm_long[structural_category%in%c("FSM", "NIC", "NNC")], aes(x=population, y=CPM, fill=population))+
#   geom_violin(trim = T)+
#   facet_wrap(~structural_category)+
#   labs(x="")+
#   scale_fill_manual(values=popcol)+
#   guides(fill="none")+
#   mytheme+
#   scale_y_continuous(trans="log10")


# 
# 
# #### REPEAT THE SAME BUT IN THIS CASE PUT THRESHOLD ON COUNTS----------------------------------------------------------------------------
# #### Do the same with trx
# # Compute cpm of gene counts
# tcounts <- column_to_rownames(tcounts, var="transcript_id")
# colnames(tcounts) <- gsub("_1","", colnames(tcounts))
# tcpm <- tcounts
# 
# # Transform to long format data.table
# tcpm <- rownames_to_column(as.data.frame(tcpm), var="transcriptid.v")
# tcpm_long <- melt(tcpm, id.vars="transcriptid.v", variable.name="cell_line_id", value.name = "count")
# setDT(tcpm_long)
# tcpm_long <- unique(data[, .(associated_gene_biotype, geneid.v, isoform,structural_category)])[tcpm_long, on=c("isoform"="transcriptid.v")]
# tcpm_long <- metadata[, .(sample, cell_line_id, population, map_reads_generalmap,total_pop_throughput)][tcpm_long, on="cell_line_id"]
# 
# 
# # plot distribution of counts to known where to put the threshold
# ggplot(tcpm_long, aes(x=count, colour = structural_category))+
#   stat_ecdf(geom = "step")+
#   mytheme+
#   facet_wrap(~associated_gene_biotype)+
#   scale_color_manual(values=colsqanti)+
#   scale_x_continuous(trans="log10", limits = c(0.01, 10e7))+
#   geom_vline(xintercept=1, linetype="dashed", color="darkgrey")+
#   labs(x="# Transcript Counts", y="ECFD", color="Structural\nCategory")+
#   theme(
#     legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
#     legend.margin = margin(0, 0, 0, 0)  # Remove extra margin around legend
#   )
# ggsave("10_figures/suppfig_quantification1/ecdf_counts_perGeneBiotype.pdf", dpi=700, width = 6, height = 2.5,  units = "in")
# 
# 
# 
# tcpm_long[, Expressed:=fifelse(count>=1, "Expressed", "Not Expressed")]
# tcpm_long[, count_expressed:=uniqueN(isoform), by=c("sample", "Expressed")]
# tcpm_long[, count_expressed_cat_biotype:=uniqueN(isoform), by=c("sample", "Expressed", "associated_gene_biotype", "structural_category")]
# tcpm_long[, count_expressed_cat:=uniqueN(isoform), by=c("sample", "Expressed", "structural_category")]
# tcpm_long[, eur:=fifelse(population%in%c("CEU", "AJI"), "European", "non-European")]
# 
# discoverylines <-unique(tcpm_long[Expressed == "Expressed", .(sample, population, structural_category, total_pop_throughput, map_reads_generalmap, count_expressed_cat)])
# discoverylines <-unique(discoverylines[, totalcount :=sum(count_expressed_cat), by="sample"][, .(sample, population,map_reads_generalmap,totalcount)])
# ggplot(discoverylines, 
#        aes(x=map_reads_generalmap/10^6, y=totalcount))+
#   geom_smooth(method = "lm", color="darkgrey")+
#   geom_point(aes( col=population))+
#   mytheme+
#   labs( y = "# Expressed PODER transcripts",  x = "Reads (M)", color="Population")+
#   scale_color_manual(values = popcol) +
#   stat_cor(label.y = 75000, size=6*0.35) +
#   stat_regline_equation(label.y = 72500, size=6*0.35)+
#   theme(
#     legend.key.size = unit(0.05, "cm"),  # Reduce size of legend keys
#     legend.margin = margin(0, 0, 0, 0)  # Remove extra margin around legend
#   )
# ggsave("10_figures/suppfig_quantification1/jitter_ExpressedTrx_PerThroughput.pdf", dpi=700, width = 2.75, height = 2,  units = "in")
# 
# 
# 
# 
# ggplot(unique(tcounts_long[, .(population, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed)]), 
#        aes(x=reorder(population, total_pop_throughput), y=count_expressed, col=Expressed))+
#   stat_summary(data = tcounts_long[Expressed == "Expressed"], 
#                aes(x = reorder(population, total_pop_throughput), y = count_expressed),
#                fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
#   geom_quasirandom(width = 0.3, aes(size=map_reads_generalmap/10^6, alpha=Expressed))+
#   mytheme+
#   ylim(c(0,max(tcounts_long$count_expressed)))+
#   labs(x="", y="# Expressed PODER transcripts", color="", size="Reads (M)", alpha="")+
#   scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
#   scale_alpha_manual(values=c(0.6, 0.2))+
#   theme(
#     legend.key.size = unit(0.1, "cm"),  # Reduce size of legend keys
#     legend.margin = margin(0, 0, 0, 0),
#     legend.position ="top"
#   )+
#   guides(
#     col = guide_legend(ncol = 1, byrow = TRUE),  # Two columns for the color legend
#     size = guide_legend(ncol = 3)  # Single column for the size legend
#   )
# ggsave("10_figures/suppfig_quantification2/jitter_ExpressedTrx_PerPopulation.pdf", dpi=700, width = 3.25, height = 2.5,  units = "in")
# 
# tcpm_long[, afr:=fifelse(population%in%c("YRI", "MPC", "LWK"), "African", "OOA")]
# tcpm_long[, expressed_transcripts_per_gene_persample:=uniqueN(isoform)/map_reads_generalmap*10^6, by=c("sample", "geneid.v")]
# ggplot(unique(tcpm_long[, .(population, sample, total_pop_throughput, afr, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
#        aes(x=afr, y=expressed_transcripts_per_gene_persample, fill=afr))+
#   geom_violin(alpha=0.7, adjust=2 )+
#   geom_boxplot(width=0.1, outliers = F)+
#   scale_fill_manual(values=c("#F7D257", "#496F5D"))+
#   mytheme+
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-2))+
#   stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
#                      method.args = list(alternative = "two.sided"))+
#   scale_y_continuous(trans="log10")+
#   labs(x="", y="# Expressed Transcripts/Sample normalized by M reads")+
#   annotation_logticks(sides = "l")+
#   guides(fill="none")
# ggplot(unique(tcpm_long[, .(population, sample, total_pop_throughput, afr, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
#        aes(x=eur, y=expressed_transcripts_per_gene_persample, fill=eur))+
#   geom_violin(alpha=0.7, adjust=2 )+
#   geom_boxplot(width=0.1, outliers = F)+
#   scale_fill_manual(values = c("#466995", "#A53860")) +
#   mytheme+
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=-2))+
#   stat_compare_means(comparisons = list(c("European", "non-European")),method = "t.test",
#                      method.args = list(alternative = "two.sided"))+
#   scale_y_continuous(trans="log10")+
#   labs(x="", y="# Expressed Transcripts/Sample normalized by M reads")+
#   annotation_logticks(sides = "l")+
#   guides(fill="none")
# ggplot(unique(tcpm_long[, .(population, sample, total_pop_throughput, eur,geneid.v, expressed_transcripts_per_gene_persample)]),
#        aes(x=population, y=expressed_transcripts_per_gene_persample, fill=population))+
#   geom_violin(alpha=0.7)+
#   geom_boxplot(width=0.1, outliers = F)+
#   scale_fill_manual(values =popcol) +
#   mytheme+
#   stat_summary(fun.data = n_fun, geom = "text", position = position_dodge(0.9),fun.args = list(y=0))
# # design <- c(
# #   "
# # AABBCC
# # #DDEE#
# # #FFGG#
# # "
# # )
# ggplot(unique(tcpm_long[, .(population,structural_category,total_pop_throughput,map_reads_generalmap, Expressed, count_expressed_cat)]), 
#        aes(x=reorder(population, total_pop_throughput), y=count_expressed_cat, col=Expressed))+
#   stat_summary(data = tcpm_long[Expressed == "Expressed"], 
#                aes(x = reorder(population, total_pop_throughput), y = count_expressed_cat),
#                fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
#   geom_quasirandom(width = 0.3, alpha=0.6, aes(size=map_reads_generalmap/10^6))+
#   mytheme+
#   labs(x="", y="# PODER transcripts", color="", size="Reads (M)")+
#   facet_wrap(~structural_category, scales="free_y")+
#   #ggh4x::facet_manual(~structural_category, scales="free_y", design = design)+
#   scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
#   theme(legend.position = c(0.65, 0.15),legend.direction="vertical",legend.box = "horizontal")
# ggsave("10_figures/suppfig_09/jitter_ExpressedTrx_PerBiotype&Population.pdf", dpi=700, width = 30, height = 13,  units = "cm")

