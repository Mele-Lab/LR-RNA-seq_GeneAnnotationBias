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

data <- fread("../novelannotations/analysis_tables/241029_num_exons_by_novelty.tsv")
data[, novelty:=factor(novelty, levels=c("Known", "Novel 5'/3'", "Novel"))]
ggplot(data, aes(x=novelty, y=n_e, fill=novelty))+
  geom_col()+
  geom_text(aes(label=paste(round(perc, digits=1), "%", sep=" ")), vjust=-0.5, size=6*0.35)+
  mytheme+
  labs(x="Internal Exon", y="# Internal PODER Exons")+
  scale_fill_manual(values=c("darkgrey", "#519D8F", "#AD5113"))+
  guides(fill="none")+
  ylim(c(0,165000))
ggsave("10_figures/fig_02/barplot_internalPODERexonsNovelty.pdf", dpi=700, width = 2, height = 2.5,  units = "in")


data <- fread("../novelannotations/analysis_tables/241108_exon_novelty_associated_gene.tsv", header = T)
data[, associated_gene_biotype:=factor(associated_gene_biotype, levels=c("Protein Coding", "lncRNA", "Intergenic", "Genic", "Fusion", "Antisense"))]
data[, novelty:=factor(novelty, levels=c("Known", "Novel 5'/3'", "Novel"))]
ggplot(data[novelty!="Known"], aes(x=associated_gene_biotype, fill=novelty))+
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  mytheme+
  labs(x="", y="# Internal PODER Exons", fill="Exon")+
  theme(legend.position = c(0.85, 0.45))+
  scale_fill_manual(values=c("#519D8F", "#AD5113"))+
  geom_text(aes(label=after_stat(count)), stat="count", position = position_dodge2(width = 0.9, preserve = "single"),vjust=-0.5, size=6*0.35)+
  geom_vline(xintercept=2.60, linetype="dashed")+
  annotate("text", x=2.75, y=6100, label="Novel/Ambiguous\nGenes", hjust=0, vjust=1, size=7*0.35, fontface = "bold", family="Helvetica")+
  annotate("text", x=2.45, y=6100, label="Annotated\nGenes", hjust=1, vjust=1, size=7*0.35, fontface = "bold", family="Helvetica")
ggsave("10_figures/01_plots/main/fig_02/barplot_internalPODERexonsNovelty_Perbiotype.pdf", dpi=700, width = 3.5, height = 2.5,  units = "in")




data <- fread("../novelannotations/analysis_tables/241031_n_supported_ics_by_novelty.tsv")
data[, structural_category:=factor(structural_category, levels=c("FSM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense"))]
data[, supported_by_external:=factor(supported_by_external, levels=rev(c("TRUE", "FALSE")))]
# ggplot(data, aes(x=structural_category, y=n_t,  fill=structural_category))+
#   geom_col(aes(alpha=supported_by_external), position = position_dodge2(width = 0.9, preserve = "single"))+
#   geom_text(aes(label=n_t, group=supported_by_external), vjust=0.5,hjust=0,position = position_dodge2(width = 0.9, preserve = "single"), size=6*0.35)+
#   mytheme+
#   labs(x="", y="# PODER Transcripts", alpha="Externally\nValidated")+
#   scale_fill_manual(values=colsqanti)+
#   scale_alpha_manual(values=rev(c(0.5, 1)))+
#   guides(fill="none",alpha = guide_legend(override.aes = list(size = 3)))+
#   theme(legend.position=c(0.8,0.25))+
#   ylim(c(0, 71000))+
#   coord_flip()
# ggsave("10_figures/01_plots/main/fig_02/barplot_PODER_validation_bysqanticat.pdf", dpi=700, width = 2.35, height = 3,  units = "in")

sumdata <- unique(data[, sumcount:=sum(n_t), by="structural_category"][, maximum:=max(n_t), by="structural_category"][, .(structural_category, maximum, sumcount)])
data[, supported_by_external:=factor(fifelse(supported_by_external==TRUE, "Supported", "Not Supported"), levels=c("Supported", "Not Supported"))]
ggplot(data, aes(x=structural_category, y=n_t,  fill=structural_category))+
  geom_col(aes(alpha=supported_by_external), position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label=n_t, group=supported_by_external), vjust=0.5,hjust=0, angle=90,position = position_dodge2(width = 0.9, preserve = "single"), size=6*0.35)+
  mytheme+
  labs(x="Structural Category compared to GENCODE v47", y="# PODER Transcripts", alpha="Support in\nOther Catalogs")+
  scale_fill_manual(values=colsqanti)+
  scale_alpha_manual(values=rev(c(0.5, 1)))+
  guides(fill="none",alpha = guide_legend(override.aes = list(size = 3)))+
  theme(legend.position=c(0.8,0.8))+
  ylim(c(0, 85000))+
  geom_text(data=sumdata, aes(x=structural_category, y=maximum+15000, label = sumcount), size=6*0.35, fontface="bold")
ggsave("10_figures/01_plots/main/fig_02/barplot_PODER_validation_bysqanticat.pdf", dpi=700, width = 3, height = 3,  units = "in")


# ggplot(data, aes(x=structural_category, y=n_t,  fill=structural_category))+
#   geom_col(aes(alpha=supported_by_external))+
#   geom_text(aes(label=n_t, group=supported_by_external),position = position_stack(vjust=0.5), size=5*0.35)+
#   mytheme+
#   labs(x="", y="# PODER Transcripts", alpha="Externally\nValidated")+
#   scale_fill_manual(values=colsqanti)+
#   scale_alpha_manual(values=rev(c(0.5, 1)))+
#   guides(fill="none",alpha = guide_legend(override.aes = list(size = 3)))+
#   theme(legend.position=c(0.8,0.25))
# ggsave("10_figures/01_plots/main/fig_02/barplot_PODER_validation_bysqanticat.pdf", dpi=700, width = 2.35, height = 3,  units = "in")



data <- fread("../novelannotations/analysis_tables/241110_tau_1_min_cpm_0.8_min_tau.tsv")

wilcoxon <-wilcox.test(data[pop_spec_t==TRUE, tau], data[pop_spec_t==FALSE, tau])
fisher.test(table(data$pop_spec_t, data$tau_pop_spec_t))
ggplot(data, aes(x = tau, color = pop_spec_t)) +
  geom_density(show.legend=FALSE) +
  stat_density(
               geom="line",position="identity", linewidth = 0.5) + 
  guides(colour = guide_legend(override.aes=list(size=1)))+
  mytheme +
  xlim(c(-0.15, 1.15)) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#393939") + #832161#5CACB1#75B9BE
  scale_color_manual(name="Transcript Discovery",values = c("darkgrey", "black"), labels=c("TRUE"="Population-Specific", "FALSE"="Across Populations")) +
  labs(x = bquote(tau ~ "\n(Transcript Expression\nPopulation-Specificity)"), y="Density")+
  theme(legend.position = "top",
        legend.key.height = unit(0.3, "lines"))+
  annotate(geom="text", label=paste0("p=",format(wilcoxon$p.value, scientific = TRUE, digits = 3)), y=1.85, x=0.1, size=6*0.35)+
  guides(color=guide_legend(nrow=2))
ggsave("10_figures/01_plots/main/fig_03/density_tauPopSpecificity_1mincpm.pdf", dpi=700, width = 1.25, height = 2,  units = "in")





data <- fread("../novelannotations/analysis_tables/241122_perc_unsupported_ics_per_annot.tsv")
data[, dataset:=factor(dataset, levels=c("PODER", "ENCODE4", "GTEx", "CHESS3"))]
data$tech <- c("Long-Reads", "Short-Reads", "Long-Reads", "Long-Reads")
data$poder <- factor(c(0,0,0,1))
data[, total:=n_ic/(perc/100)]
data[, n_ic_supported:=n_ic/(perc/100)*(1-perc/100)]
data[, V1:=NULL]
mydata <- melt(data, id.vars = c("dataset", "poder", "total", "tech", "perc"), variable.name = "Support", value.name = "Number")
mydata[, perc:=round(Number/total*100, 1)]
mydata[, Support:=fifelse(Support=="n_ic", "Not Supported", "Supported")]


ggplot(mydata, aes(x=dataset, y=Number, fill=tech)) +
  geom_col(aes(alpha=Support)) +
  geom_text(aes(label=paste(perc, "%", sep=" "), group = Support), position=position_stack(vjust=0.35), size=6*0.35, show.legend = F) +
  geom_text(aes(label=Number, group = Support), position=position_stack(vjust=0.5), size=6*0.35, show.legend = F) +
  mytheme +
  labs(x="", y="# Novel Transcripts", alpha="External Catalog Support", fill="Technology") +
  scale_fill_manual(values=c("purple", "#4A90E2")) +
  scale_alpha_manual(values=c(0.6, 1)) +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        legend.position="top") +
  guides(fill = guide_legend(nrow = 1, title.position = "top"),
         alpha = guide_legend(nrow = 1, title.position = "top"))

ggsave("10_figures/01_plots/supp/12_external_supp/barplot_PODER_ENCODE_GTEX_CHESS_nonvalidation.pdf", dpi=700, width = 3, height = 2.25,  units = "in")


chess <- fread("../novelannotations/analysis_tables/chess_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="CHESS3"]
encode <- fread("../novelannotations/analysis_tables/enc_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="ENCODE4"]
gtex <- fread("../novelannotations/analysis_tables/gtex_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="GTExv9"]
poder <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")[, .(isoform, structural_category)][, effort:="PODER"]

mydata <- rbind.data.frame(rbind.data.frame(chess, encode), gtex)
categories <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")
names(categories) <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "intergenic", "genic", "fusion", "antisense")
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")

mydata[, structural_category:=categories[structural_category]]
mydata[, structural_category:=factor(structural_category, levels=categories)]
mydata <- rbind.data.frame(mydata, poder)
mydata$effort <- factor(mydata$effort,  levels = c("PODER", "GTExv9","ENCODE4", "CHESS3"))



ggplot(mydata[!is.na(structural_category)], aes(x=structural_category, fill=structural_category))+
  facet_wrap(~effort, nrow=1)+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  guides(fill="none")+
  labs(x="", y="# Transcripts")+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))+
  geom_text(aes(label = after_stat(count)), stat="count", angle=90, size=6*0.35, hjust=0)+
  ylim(c(0, 180000))
ggsave("10_figures/01_plots/supp/12_external_supp/barplot_ENCODE_GTEX_CHESS_sqanti.pdf", dpi=700, width = 7, height = 3,  units = "in")






data <- fread("../novelannotations/analysis_tables/241121_aa_per_gene_per_sample.tsv")
data[, ooa:=fifelse(ooa=="AFR", "African", "OOA")]
p<-ggplot(unique(data[, .(gid, ooa, population, n_aa_norm)]),
          aes(x=n_aa_norm, color=ooa))+
  stat_ecdf(geom="step")+
  scale_color_manual(values=c("#F7D257", "#496F5D"))+
  mytheme+
  ggmagnify::geom_magnify(from = c(xmin = 0.25, xmax = 1.25, ymin = 0.75, ymax = 1),
                          to = c(xmin =1.75, xmax = 6.4, ymin = 0.1, ymax = 0.85))+
  labs(x="# Predicted Proteins from Expressed Transcripts/Gene (and sample)\nNormalized by Million Reads", y="Cumulative Proportion", color="")+
  theme(legend.position="top",
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  annotate(geom="text", label=paste0("p=",wilcox.test(unique(data[, .(gid, ooa, population, n_aa_norm)])[ooa=="OOA", n_aa_norm],
                                          unique(data[, .(gid, ooa, population, n_aa_norm)])[ooa=="African", n_aa_norm])$p.val),
           x=4, y=0.4, size=6*0.35) # pvalue obtained from below plot
ggsave(plot=p,"10_figures/01_plots/supp/22_afr_trx/ecdf_Expressed_aa_PerGene_AFROOA_sample.pdf", dpi=700, width = 3, height = 3,  units = "in")

data <- fread("../novelannotations/analysis_tables/241121_aa_per_gene_per_pop.tsv")
data[, ooa:=fifelse(ooa=="AFR", "African", "OOA")]
p<-ggplot(unique(data[, .(gid, ooa, population, n_aa_norm)]),
          aes(x=n_aa_norm, color=ooa))+
  stat_ecdf(geom="step")+
  scale_color_manual(values=c("#F7D257", "#496F5D"))+
  mytheme+
  ggmagnify::geom_magnify(from = c(xmin = 0.25, xmax = 1.25, ymin = 0.75, ymax = 1),
                          to = c(xmin =1.75, xmax = 4.75, ymin = 0.1, ymax = 0.85))+
  labs(x="# Predicted Proteins from Expressed Transcripts/Gene (and population)\nNormalized by Million Reads", y="Cumulative Proportion", color="")+
  theme(legend.position="top",
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  annotate(geom="text", label=paste0("p=",wilcox.test(unique(data[, .(gid, ooa, population, n_aa_norm)])[ooa=="OOA", n_aa_norm],
                                                      unique(data[, .(gid, ooa, population, n_aa_norm)])[ooa=="African", n_aa_norm])$p.val),
           x=4, y=0.4, size=6*0.35) # pvalue obtained from below plot
ggsave(plot=p,"10_figures/01_plots/supp/22_afr_trx/ecdf_Expressed_aa_PerGene_AFROOA_pop.pdf", dpi=700, width = 3, height = 3,  units = "in")


# personalized GRCh38
metadata <- fread("00_metadata/data/samples_metadata_online.tsv")
popcol <- unique(metadata$color_pop)
names(popcol) <- unique(metadata$population)

data <- fread("10_figures/data/250210_perc_novel_hg38_absent_sjs_w_variant_per_cell_line.tsv")
data[, cell_line_id:=gsub("^..", "", cell_line_id)]
metadata[, cell_line_id:=gsub("^..", "", cell_line_id)]

data <- unique(metadata[, .(cell_line_id, population)])[data, on="cell_line_id"]

ggplot(data, aes(x="", y=perc))+
  geom_violin(fill="#61814B",alpha=0.5)+
  geom_boxplot(width=0.3,fill="#61814B")+
  ggbeeswarm::geom_quasirandom(size=0.5, aes(color=population))+
  mytheme+
  labs(x="", y="% Novel Splice Junctions\nfound only in personalized-GRCh38\nexplained by tools filters", color="")+
  theme(axis.ticks.x=element_blank(), legend.position = "top")+
  scale_color_manual(values=popcol)
ggsave("10_figures/01_plots/main/personalizedhg38/violin_percentage_NovelSJ_explained.pdf", dpi=700, width = 1.5, height = 2.25,  units = "in")


data <- fread("../novelannotations/analysis_tables/250221_personal_hg38_unique_novel_explainability.tsv")
data[, explained:=fifelse(`Exonic variant\n+-10bp from SJ` | `SS variant`, TRUE, FALSE)]
data[, total:=sum(n_sj), by=.(cell_line_id)]
data[, sum:=sum(n_sj), by=.(cell_line_id, explained)]
data[, per:=sum/total*100]
ggplot(unique(data[explained==TRUE]), aes(x="", y=per))+
  geom_violin(fill="#61814B",alpha=0.5)+
  geom_boxplot(width=0.3,fill="#61814B")+
  ggbeeswarm::geom_quasirandom(size=0.5)+
  mytheme+
  labs(x="", y="% Novel Splice Junctions\nfound only in personalized-GRCh38\nexplained by tools filters")+
  theme(axis.ticks.x=element_blank())

#### MAGE tau

data <- fread("../novelannotations/analysis_tables/250213_mage_tau_biotype_pop_spec.tsv")
ggplot(data[!associated_gene_biotype %in% c("Other", "Pseudogene")], 
       aes(x = tau, color = pop_spec_t)) +
  geom_density() +
  facet_wrap(~associated_gene_biotype) +
  mytheme +
  scale_color_manual(values = c("black", "grey", "darkred"), 
                     labels = c("Absent from PODER" = "Absent from\nPODER", 
                                "False" = "Population-Shared", 
                                "True" = "Population-Specific")) +
  labs(x = bquote(tau ~ "\n(Transcript Expression\nPopulation-Specificity)"), 
       y = "Density", color = "") +
  guides(color = guide_legend(override.aes = list(linetype = 1, shape = 2, size = 5)))+
  theme(legend.position = "top")
ggsave("10_figures/01_plots/supp/MAGE_tau.pdf", dpi=700, width = 4, height = 2.25,  units = "in")



### FST

data <- fread("../novelannotations/analysis_tables/250219_novel_exon_fsts.tsv")

novelpval <- format(wilcox.test(data[novelty=="Novel", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)
novelfivethree <- format(wilcox.test(data[novelty=="Novel 5'/3'", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)
p <- ggplot(data, aes(mean_mean_fst, col=novelty)) + 
  stat_ecdf(geom = "step")+
  mytheme+
  scale_color_manual(values=c("darkgrey", "#AD5113","#519D8F"))+
  labs(x=expression("Mean " ~ F[ST] ~ "(CEU-All Populations)"), color="Exonic Region", y="Cumulative Proportion")+
  ggmagnify::geom_magnify(from = c(xmin = 0.025, xmax = 0.1, ymin = 0.75, ymax = 1), 
                          to = c(xmin = 0.15, xmax = 0.65, ymin = 0.1, ymax = 0.75))+
  annotate(geom="text", label=paste0("p=",novelpval), x=0.5, y=0.15, color="#AD5113")+
  annotate(geom="text", label=paste0("p=",novelfivethree), x=0.5, y=0.25, color="#519D8F")
ggsave(plot=p,"10_figures/01_plots/supp/10_fst/ECDF_fstCEUYRI.YRIdiscoveredExons.pdf", dpi=700, width = 6, height = 4,  units = "in")



# FST popsp exons
data <- fread("../novelannotations/analysis_tables/250221_pop_spec_exon_fsts.tsv")

novelpval <- format(wilcox.test(data[novelty=="Novel", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)
novelfivethree <- format(wilcox.test(data[novelty=="Novel 5'/3'", mean_mean_fst], data[novelty=="Known", mean_mean_fst])$p.value, scientific = TRUE, digits = 2)
p <- ggplot(data, aes(mean_mean_fst, col=pop_spec)) + 
  stat_ecdf(geom = "step")+
  mythemen+
  scale_color_manual(values=c("darkgrey", "#AD5113","#519D8F"))+
  labs(x=expression("Mean " ~ F[ST] ~ "(CEU-All Populations)"), color="Exonic Region", y="Cumulative Proportion")+
  # ggmagnify::geom_magnify(from = c(xmin = 0.025, xmax = 0.1, ymin = 0.75, ymax = 1), 
  #                         to = c(xmin = 0.15, xmax = 0.65, ymin = 0.1, ymax = 0.75))+
  # annotate(geom="text", label=paste0("p=",novelpval), x=0.5, y=0.15, color="#AD5113")+
  # annotate(geom="text", label=paste0("p=",novelfivethree), x=0.5, y=0.25, color="#519D8F")+
  facet_wrap(~pop_spec_pop)
ggsave(plot=p,"10_figures/01_plots/supp/10_fst/ECDF_fstCEUYRI.YRIdiscoveredExons.pdf", dpi=700, width = 6, height = 4,  units = "in")




##### SJ explainability

data <- fread("../novelannotations/analysis_tables/250221_personal_hg38_unique_novel_explainability.tsv")
data <- as.data.frame(data)
# Load packages
library(ComplexUpset)


median_fun <- function(x) {
  return(data.frame(label = round(median(x), 2)))  # Return median value for annotation
}

# Define the set columns
set_columns <- c("SS variant","Exonic variant\n+-10bp from SJ")

# Create an intersection identifier column
data$intersection <- apply(data[set_columns], 1, paste, collapse = "")
data$intersection <- factor(data$intersection, levels = rev(c("FALSEFALSE", "FALSETRUE", "TRUEFALSE", "TRUETRUE")))
# Create the UpSet plot with boxplots
upset(
  data, 
  sort_intersections_by='degree',
  set_sizes = FALSE,
  intersect = set_columns,
  base_annotations = list(),
  annotations = list(
    'Value Distribution' = (
      ggplot(data, aes(x = intersection, y = `% SJs`)) +
        geom_boxplot(fill = "#61814B", linewidth=0.25) +
        mytheme +
        labs(x = "", y = "% of Novel Splice Junctions\nExplained by Tool Filters") +
        theme(axis.text.x = element_blank()) +
        stat_summary(
          fun = median, 
          geom = "text", 
          aes( label = paste0(round( after_stat(y), 0), "%")),  # Round the median and use as label
          vjust = -2, size=6*0.35
        )+
        coord_cartesian(ylim=c(0, 74))
      
    )
  )
)+
  xlab("")+
  theme(text=element_text(family="Helvetica", size=7, face="bold"))
ggsave("10_figures/01_plots/main/personalizedhg38/upset_percentageExplainedNovelSJ.pdf", dpi=700, width = 3, height = 3,  units = "in")














# Define the plot with boxplot and median annotations
ggplot(data, aes(x = intersection, y = `% SJs`)) +
  geom_boxplot(fill = "#61814B") +
  mytheme +
  labs(x = "", y = "% of Novel Splice Junctions\nExplained by Tool Filters") +
  theme(axis.text.x = element_blank()) +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes( label = paste0(round( after_stat(y), 0), "%")),  # Round the median and use as label
    vjust = -2
  )+
  coord_cartesian(ylim=c(0, 70))

# Define the plot with boxplot and median annotations using `n_sj_w_explanation`
ggplot(data, aes(x = intersection, y = `% SJs`)) +
  geom_boxplot(fill = "#61814B") +
  mytheme +
  labs(x = "", y = "% of Novel Splice Junctions\nExplained by Tool Filters") +
  theme(axis.text.x = element_blank()) +
  stat_summary(
    fun = median, 
    geom = "point", 
    aes(y = after_stat(y)),  # Compute median for the `y` position
    size = 3, color = "red"  # Use a point for median
  ) +
  geom_text(
    aes(
      label = paste0(round(n_sj_w_explanation, 0), "%"),
      y = max(`% SJs`, na.rm = TRUE) + 5  # Position the label above the max value of `% SJs`
    ),
    vjust = 0  # Adjust vertical alignment of the text
  )

