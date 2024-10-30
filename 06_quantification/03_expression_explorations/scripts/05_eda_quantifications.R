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
gcounts <- column_to_rownames(gcounts, var="geneid.v")
gcounts <- as.matrix(gcounts)
gcpm <- cpm(gcounts)

# Transform to long format data.table
gcpm <- gcpm[rowSums(gcpm>=0.1)>0,]
gcpm <- rownames_to_column(as.data.frame(gcpm), var="geneid.v")
gcpm_long <- melt(gcpm, id.vars="geneid.v", variable.name="sample", value.name = "CPM")
setDT(gcpm_long)
gcpm_long <- unique(data[, .(associated_gene_biotype, geneid.v)])[gcpm_long, on=c("geneid.v")]
gcpm_long <- metadata[, .(sample, population, map_reads_generalmap,total_pop_throughput)][gcpm_long, on="sample"]
gcpm_long[, Expressed:=fifelse(CPM>=0.1, "Expressed", "Not Expressed")]
gcpm_long[, count_expressed:=uniqueN(geneid.v), by=c("sample", "Expressed", "associated_gene_biotype")]

ggplot(unique(gcpm_long[, .(population,associated_gene_biotype, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed, col=Expressed))+
  stat_summary(data = gcpm_long[Expressed == "Expressed"], 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed),
               fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
  geom_quasirandom(width = 0.3, alpha=0.6, aes(size=map_reads_generalmap/10^6))+
  mytheme+
  ylim(c(0,max(gcpm_long$count_expressed)))+
  labs(x="", y="# PODER genes", color="", size="Reads (M)")+
  facet_wrap(~associated_gene_biotype)+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  theme(legend.position = "top")
ggsave("10_figures/suppfig_09/jitter_ExpressedGenes_PerBiotype&Population.pdf", dpi=700, width = 30, height = 13,  units = "cm")



#### Do the same with trx
# Compute cpm of gene counts
tcounts <- column_to_rownames(tcounts, var="transcript_id")
colnames(tcounts) <- gsub("_1","", colnames(tcounts))
tcounts <- as.matrix(tcounts)
tcpm <- cpm(tcounts)

# Transform to long format data.table
tcpm <- tcpm[rowSums(tcpm>=0.1)>0,]
tcpm <- rownames_to_column(as.data.frame(tcpm), var="transcriptid.v")
tcpm_long <- melt(tcpm, id.vars="transcriptid.v", variable.name="cell_line_id", value.name = "CPM")
setDT(tcpm_long)
tcpm_long <- unique(data[, .(associated_gene_biotype, geneid.v, isoform,structural_category)])[tcpm_long, on=c("isoform"="transcriptid.v")]
tcpm_long <- metadata[, .(sample, cell_line_id, population, map_reads_generalmap,total_pop_throughput)][tcpm_long, on="cell_line_id"]
tcpm_long[, Expressed:=fifelse(CPM>=0.1, "Expressed", "Not Expressed")]
tcpm_long[, count_expressed:=uniqueN(isoform), by=c("sample", "Expressed")]
tcpm_long[, count_expressed_cat_biotype:=uniqueN(isoform), by=c("sample", "Expressed", "associated_gene_biotype", "structural_category")]
tcpm_long[, count_expressed_cat:=uniqueN(isoform), by=c("sample", "Expressed", "structural_category")]
tcpm_long[, eur:=fifelse(population%in%c("CEU", "AJI"), "European", "non-European")]

ggplot(unique(tcpm_long[, .(population, total_pop_throughput,map_reads_generalmap, Expressed, count_expressed)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed, col=Expressed))+
  stat_summary(data = tcpm_long[Expressed == "Expressed"], 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed),
               fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
  geom_quasirandom(width = 0.3, alpha=0.6, aes(size=map_reads_generalmap/10^6))+
  mytheme+
  ylim(c(0,max(tcpm_long$count_expressed)))+
  labs(x="", y="# PODER transcripts", color="", size="Reads (M)")+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  theme(legend.position = "top")
ggsave("10_figures/suppfig_09/jitter_ExpressedTrx_PerPopulation.pdf", dpi=700, width = 15, height = 13,  units = "cm")

# design <- c(
#   "
# AABBCC
# #DDEE#
# #FFGG#
# "
# )
ggplot(unique(tcpm_long[, .(population,structural_category,total_pop_throughput,map_reads_generalmap, Expressed, count_expressed_cat)]), 
       aes(x=reorder(population, total_pop_throughput), y=count_expressed_cat, col=Expressed))+
  stat_summary(data = tcpm_long[Expressed == "Expressed"], 
               aes(x = reorder(population, total_pop_throughput), y = count_expressed_cat),
               fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 1)+
  geom_quasirandom(width = 0.3, alpha=0.6, aes(size=map_reads_generalmap/10^6))+
  mytheme+
  labs(x="", y="# PODER transcripts", color="", size="Reads (M)")+
  facet_wrap(~structural_category, scales="free_y")+
  #ggh4x::facet_manual(~structural_category, scales="free_y", design = design)+
  scale_color_manual(values=rev(c("#c44536", "#457b9d")))+
  theme(legend.position = c(0.65, 0.15),legend.direction="vertical",legend.box = "horizontal")
ggsave("10_figures/suppfig_09/jitter_ExpressedTrx_PerBiotype&Population.pdf", dpi=700, width = 30, height = 13,  units = "cm")





ggplot(unique(tcpm_long[Expressed=="Expressed", .(eur, population,structural_category,total_pop_throughput,map_reads_generalmap, count_expressed_cat)]), 
       aes(x=eur, y=count_expressed_cat, fill=eur))+
  geom_violin(alpha=0.75)+
  geom_boxplot(outliers = F, width=0.05)+
  geom_quasirandom(width = 0.3, alpha=0.6, aes(size=map_reads_generalmap/10^6, color=population))+
  mytheme+
  labs(x="", y="# PODER transcripts", color="", size="Reads (M)")+
  facet_wrap(~structural_category)+
  scale_color_manual(values=popcol)+
  scale_fill_manual(values=c("#466995", "#A53860"))+
  geom_pwc(ref.group="European",
           method="t_test")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = 30000), 
               position = position_dodge(0.9))
  