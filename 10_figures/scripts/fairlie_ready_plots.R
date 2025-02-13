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


data <- fread("../novelannotations/analysis_tables/241031_perc_supported_ics_per_annot.tsv")
data[, dataset:=factor(dataset, levels=c("PODER", "ENCODE4", "GTEx", "CHESS3"))]
data$tech <- c("Long-Reads", "Short-Reads", "Long-Reads", "Long-Reads")
data$poder <- factor(c(0,0,0,1))
ggplot(data, aes(x=dataset, y=n_ic, fill=tech))+
  geom_col(aes(alpha=poder))+
  geom_text(aes(label=paste(round(perc, digits=1), "%", sep=" ")), vjust=-0.5)+
  mytheme+
  labs(x="", y="# Externally-supported\nNovel Transcripts", fill="")+
  scale_fill_manual(values=c("purple", "#4A90E2"))+
  scale_alpha_manual(values=c(0.6, 1))+
  guides(alpha="none")+
  theme(legend.position="top")+
  ylim(c(0,19000))
ggsave("10_figures/suppfig/barplot_PODER_ENCODE_GTEX_CHESS_validation.pdf", dpi=700, width = 11, height = 10,  units = "cm")



data <- fread("../novelannotations/analysis_tables/241031_n_supported_ics_by_novelty.tsv")
data[, structural_category:=factor(structural_category, levels=rev(c("FSM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")))]
data[, supported_by_external:=factor(supported_by_external, levels=c("TRUE", "FALSE"))]
ggplot(data, aes(x=structural_category, y=n_t,  fill=structural_category))+
  geom_col(aes(alpha=supported_by_external), position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label=n_t, group=supported_by_external), vjust=0.5,hjust=0,position = position_dodge2(width = 0.9, preserve = "single"), size=6*0.35)+
  mytheme+
  labs(x="", y="# PODER Transcripts", alpha="Externally\nValidated")+
  scale_fill_manual(values=colsqanti)+
  scale_alpha_manual(values=rev(c(0.5, 1)))+
  guides(fill="none",alpha = guide_legend(override.aes = list(size = 3)))+
  theme(legend.position=c(0.8,0.25))+
  ylim(c(0, 71000))+
  coord_flip()
ggsave("10_figures/01_plots/main/fig_02/barplot_PODER_validation_bysqanticat.pdf", dpi=700, width = 2.35, height = 3,  units = "in")



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
ggplot(data, aes(x=dataset, y=n_ic, fill=tech))+
  geom_col(aes(alpha=poder))+
  geom_text(aes(label=paste(round(perc, digits=1), "%", sep=" ")), vjust=-0.5, size=6*0.35)+
  geom_text(aes(label=n_ic), vjust=1.5, size=6*0.35)+
  mytheme+
  labs(x="", y="# Novel Transcripts\nnot Supported by External Effort", fill="")+
  scale_fill_manual(values=c("purple", "#4A90E2"))+
  scale_alpha_manual(values=c(0.6, 1))+
  guides(alpha="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        legend.position="top")+
  ylim(c(0,60000))
ggsave("10_figures/01_plots/supp/12_external_supp/barplot_PODER_ENCODE_GTEX_CHESS_nonvalidation.pdf", dpi=700, width = 3, height = 2.25,  units = "in")


chess <- fread("../novelannotations/analysis_tables/chess_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="CHESS3"]
encode <- fread("../novelannotations/analysis_tables/enc_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="ENCODE4"]
gtex <- fread("../novelannotations/analysis_tables/gtex_sqanti_classification.txt")[, .(isoform, structural_category)][, effort:="GTExv9"]

mydata <- rbind.data.frame(rbind.data.frame(chess, encode), gtex)
categories <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")
names(categories) <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "intergenic", "genic", "fusion", "antisense")
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")

mydata[, structural_category:=categories[structural_category]]
mydata[, structural_category:=factor(structural_category, levels=categories)]

ggplot(mydata[!is.na(structural_category)], aes(x=structural_category, fill=structural_category))+
  facet_wrap(~effort)+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  guides(fill="none")+
  labs(x="", y="# Transcripts")+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))+
  geom_text(aes(label = after_stat(count)), stat="count", angle=90, size=6*0.35, hjust=0)+
  ylim(c(0, 180000))
ggsave("10_figures/01_plots/supp/12_external_supp/barplot_ENCODE_GTEX_CHESS_sqanti.pdf", dpi=700, width = 6, height = 3,  units = "in")






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

data <- fread("../novelannotations/analysis_tables/250210_perc_novel_hg38_absent_sjs_w_variant_per_cell_line.tsv")
ggplot(data, aes(x="", y=perc))+
  geom_violin(fill="#61814B",alpha=0.5)+
  geom_boxplot(width=0.3,fill="#61814B")+
  ggbeeswarm::geom_quasirandom(size=0.5)+
  mytheme+
  labs(x="", y="% Novel Splice Junctions\nfound only in personalized-GRCh38\nexplained by tools filters")+
  theme(axis.ticks.x=element_blank())
ggsave("10_figures/01_plots/main/personalizedhg38/violin_percentage_NovelSJ_explained.pdf", dpi=700, width = 1.5, height = 2.25,  units = "in")
