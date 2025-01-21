## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/08_allele_specifics"
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
library(ggpmisc)
library(ggpubr)
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadataraw <- metadataraw[mixed_samples==FALSE]
metadata <- metadataraw[merged_run_mode==TRUE]
metadata[, samplecode:=paste(lab_number_sample, lab_sampleid, cell_line_id, sep="_")]

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
# load data
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
arrayenh <- fread("array_enhanced_gencode", header = F)
array <- rbind.data.frame(arraygen, arraypan)
array <- rbind.data.frame(array, arrayenh)

asts <- list()

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  temp <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_annotated_taggedresults.tsv"))
  temp[, `:=`(annot=TYPE, samplecode=SAMPLE)]
  asts <- append(asts, list(temp))}

astsraw <- rbindlist(asts, use.names=TRUE)
astsraw <- metadata[, .(samplecode, sample, population, map_reads_assemblymap)][astsraw, on="samplecode"]
astsraw[, annot := ifelse(annot=="gencode", "GENCODE", 
                          ifelse(annot=="enhanced_gencode", "Enhanced\nGENCODE", "PODER"))]
astsraw[, annot:=factor(annot, levels=c("GENCODE", "PODER", "Enhanced\nGENCODE"))]
# fwrite(astsraw, "data/ASTS_results_threeannots.tsv", quote=F, row.names = F, sep="\t")
astsraw <- fread("data/ASTS_results_threeannots.tsv")
asts <- astsraw[gene_testable==TRUE & FDR<0.05]
asts_test <- astsraw[gene_testable==TRUE]

# computed number of tested genes
asts[, tested_genes:=uniqueN(geneid.v), by=c("sample", "annot")]
asts_test[, tested_genes:=uniqueN(geneid.v), by=c("sample", "annot")]


asts[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "OOA")]
asts_test[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "OOA")]

# Significant genes per ancestry
ggplot(unique(asts[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=afr, y=tested_genes, fill=afr))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  scale_size_continuous(range = c(0.5, 1.5))+
  mytheme+
  labs(x="", y="# Significant ASTU Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=F),size=6*0.35)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  guides(fill="none")+
  facet_wrap(~annot)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("../10_figures/01_plots/supp/30_as_afr/violin_sigASTU_afrOOA.pdf", dpi=700, width = 4, height = 4,  units = "in")
ggplot(unique(asts_test[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=afr, y=tested_genes, fill=afr))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  scale_size_continuous(range = c(0.5, 1.5))+
  mytheme+
  labs(x="", y="# Tested ASTU Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE),size=6*0.35)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  guides(fill="none")+
  facet_wrap(~annot)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("../10_figures/01_plots/supp/30_as_afr/violin_ASTU_afrOOA.pdf", dpi=700, width = 4, height = 4,  units = "in")
ggplot(unique(asts[annot=="PODER", .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=afr, y=tested_genes, fill=afr))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.65)+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  scale_size_continuous(range = c(0.5, 1.5))+
  mytheme+
  labs(x="", y="# Significant ASTU Genes with PODER", size="Mapped\nReads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE),size=6*0.35)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  guides(fill="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  ylim(c(20, 85))
ggsave("../10_figures/01_plots/main/fig_04/violin_sigASTU_afrOOA_PODER.pdf", dpi=700, width = 3, height = 2.5,  units = "in")
# number of tested genes per annot?

ggplot(unique(asts[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=tested_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.5)+
  scale_size_continuous(range = c(0.5, 1.5))+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Significant ASTU Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  geom_pwc(method="t_test", label.size=7*0.35, p.adjust.by ="panel",label = "p.adj.format",step.increase = 0.16)+
  # stat_compare_means(comparisons = list(c("PODER", "GENCODEv47"), 
  #                                       c("Enhanced\nGENCODEv47", "PODER"),
  #                                       c("Enhanced\nGENCODEv47", "GENCODEv47")),method = "t.test",
  #                    method.args = list(alternative = "two.sided", paired=TRUE), size=6*0.35,p.adjust.method = "holm",label = "p.adj")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("GENCODEv47"="GENCODE", "Enhanced\nGENCODEv47"="Enhanced\nGENCODE"))+
  ylim(c(15, 125))
ggsave("../10_figures/01_plots/main/fig_04/violin_sigASTU_byAnnot.pdf", dpi=700, width = 2.75, height = 2.75,  units = "in")
ggsave("../10_figures/01_plots/supp/31_as_annots/violin_sigASTU_byAnnot.pdf", dpi=700, width = 4, height = 4,  units = "in")
ggplot(unique(asts_test[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=tested_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  scale_size_continuous(range = c(0.5, 1.5))+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Tested ASTU Genes", size="Mapped Reads (M)", col="Population", fill="")+
  geom_pwc(method="t_test", label.size=7*0.35, p.adjust.by ="panel",label = "p.adj.format",step.increase = 0.16)+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("GENCODEv47"="GENCODE", "Enhanced\nGENCODEv47"="Enhanced\nGENCODE"))
ggsave("../10_figures/01_plots/supp/31_as_annots/violin_testedASTU_byAnnot.pdf", dpi=700, width = 4, height = 4,  units = "in")




# Explore how many more genes are tested in each annot
astswide <- dcast(unique(asts_test[, .(sample, annot, tested_genes)]), sample ~ annot, value.var = "tested_genes") 
astswide[, diff_tested_genes:=PODER-GENCODE]
astswide[, ratio_tested_genes:=PODER/GENCODE]

astswidemeta <- metadata[, .(sample, population, map_reads_assemblymap)][astswide, on="sample"]
astswidemeta[, eur:=ifelse(population %in% c("LWK", "YRI"), "AFR",
                           ifelse(population %in% c("PEL", "HAC", "ITU"), "nonEUR-OOA", "EUR"))]
# ggplot(astswidemeta, aes(x=sample, y=diff_tested_genes, fill=population))+
#   geom_col()+
#   mytheme+
#   scale_fill_manual(values=popcols)


# check  normality and homocedasticity
car::leveneTest(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
bartlett.test(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
shapiro.test(astswidemeta$ratio_tested_genes)
qqnorm(astswidemeta$ratio_tested_genes)
qqline(astswidemeta$ratio_tested_genes, col = "red")
my_comparisons <- list( c("CEU", "HAC"), c("CEU", "ITU"), c("CEU", "LWK"),c("CEU", "PEL"), c("CEU", "YRI") )


anova_result <- aov(ratio_tested_genes ~ population, data = unique(astswidemeta[, .(sample, population, ratio_tested_genes)]))
summary(anova_result)
# Post-hoc Tukey Test
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
astswidemeta[, eur:=ifelse(population=="CEU", "European", "Non-European")]
astswidemeta[, mean:=mean(ratio_tested_genes), by=c("eur")]
astswidemeta[, sd:=sd(ratio_tested_genes), by=c("eur")]


# Create the ggplot barplot
ggplot(unique(astswidemeta[, .(eur, mean, sd)]), aes(x = eur, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "#545454") +
  geom_bar(stat = "identity", aes(fill = eur), alpha = 0.75) +
  labs(y = "Ratio of Tested ASTU Genes\nPODER vs GENCODE", x = "", color="Population") +
  mytheme +
  scale_fill_manual(values = c("#466995", "#A53860")) +
  ggbeeswarm::geom_quasirandom(data = astswidemeta, aes(x = eur, y = ratio_tested_genes, colour = population), alpha = 0.5, size = 1) +
  geom_segment(aes(x = 1, xend = 2, y = 1.38, yend = 1.38), color = "black", linewidth=1*0.35) + # Add bar manually
  annotate("text", x = 1.5, y = 1.4, label = paste0("p=", round(t.test(astswidemeta[eur=="European", ratio_tested_genes], astswidemeta[eur=="Non-European", ratio_tested_genes])$p.val, 4)),
           size = 6*0.35, family="Helvetica") + # Add p-value manually
  guides(fill = "none") +
  scale_color_manual(values=popcols)+
  coord_cartesian(y = c(1, 1.4)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("../10_figures/01_plots/main/fig_04/barplot_testedASTU_increase.pdf", dpi=700, width = 2.6, height = 2.25,  units = "in")


# if mean and sd is computed by pop
ggplot(unique(astswidemeta[, .(population, mean, sd)]), aes(x = population, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "#545454") + # Error bars
  geom_bar(stat = "identity", aes(fill = population), alpha=0.75) + # Bar plot
  labs(y = "Ratio Tested ASTU Genes with PODER  Transcripts", x="") +
  mytheme+
  scale_fill_manual(values=popcols)+
  ggbeeswarm::geom_quasirandom(data=astswidemeta,aes(x=population, y=ratio_tested_genes), alpha=0.5, size=0.5)+
  guides(fill="none")+
  coord_cartesian(y=c(1, 1.375))+
  geom_hline(yintercept=1, linetype="dashed", color="black")+
  geom_hline(yintercept=astswidemeta[population=="CEU", mean], linetype="dashed", color="darkgrey")



# ggplot(astswidemeta, aes(x=population, y=ratio_tested_genes))+
#   mytheme+
#   geom_boxplot(outliers=F)+
#   geom_jitter(alpha=0.8, aes(col=population, size=map_reads_assemblymap/10^6))+
#   scale_color_manual(values=popcols)+
#   labs(x="", y="Ratio of Tested Genes by PODER against GENCODE", size="Mapped Reads (M)", col="Population")+
#   stat_compare_means(comparisons = my_comparisons,method = "t.test",
#                      method.args = list(alternative = "less"))+
#   geom_hline(yintercept=1.1)
# 
# my_comparisons <- list( c("Eur", "non-Eur") )
# astswidemeta[, eur:=fifelse(population=="CEU", "Eur", "non-Eur")]
# ggplot(astswidemeta, aes(x=eur, y=ratio_tested_genes))+
#   mytheme+
#   geom_boxplot(outliers=F)+
#   geom_jitter(alpha=0.8, aes(col=population, size=map_reads_assemblymap/10^6))+
#   scale_color_manual(values=popcols)+
#   labs(x="", y="Ratio of significant genes by PODER against GENCODE", size="Mapped Reads (M)", col="Population")+
#   stat_compare_means(comparisons = my_comparisons,method = "t.test",
#                      method.args = list(alternative = "less"))+
#   geom_hline(yintercept=1.1)
# 
# ggplot(astswidemeta, aes(x=population, y=map_reads_assemblymap/10^6))+
#   mytheme+
#   geom_boxplot(outliers=F)+
#   geom_jitter(alpha=0.8, aes(col=population, size=ratio_tested_genes))+
#   scale_color_manual(values=popcols)+
#   labs(x="", y="Mapped reads (M)", size="Ratio tested genes\nPODER/GENCODE", col="Population")+
#   stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",
#                      method.args = list(alternative = "less")) 

# check biotypes of tested genes
astsgenes <- unique(asts[, .(annot, sample, population, map_reads_assemblymap, geneid.v, gene_biotype, gene_testable)])
astsgenes[gene_biotype=="novel/ambiguous gene", gene_biotype:="Novel/Ambiguous Gene"][gene_biotype=="protein_coding", gene_biotype:="Protein Coding"]
ggplot(unique(astsgenes[, .(annot, geneid.v, gene_biotype)]), aes(x=annot, fill=gene_biotype))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c("purple","#ee9b00","darkgrey","#4C8C36"))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)+
  scale_x_discrete(labels = c("gencode" = "GENCODEv47", "poder"="PODER"))+
  labs(x="", y="# ASTU significant Genes", fill="Gene Biotype")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 2))+
  scale_x_discrete(labels=c("GENCODEv47"="GENCODE", "Enhanced\nGENCODEv47"="Enhanced\nGENCODE"))+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))

ggsave("../10_figures/01_plots/supp/32_astu_biotype/barplot_sigASTU_geneBiotypes.pdf", dpi=700, width = 3, height = 3,  units = "in")
ggsave("../10_figures/02_panels/svg/supp/32_astu_biotype.svg", dpi=700, width = 3, height = 3,  units = "in")
ggsave("../10_figures/02_panels/png/supp/32_astu_biotype.png", dpi=700, width = 3, height = 3,  units = "in")

# which genes are tested with PODER and GENCODE
# astsgenesannot[, tested:="Tested"]

astsgenesannotwide <- dcast(unique(astsgenes[annot%in%c("GENCODEv47", "PODER"), .(annot, geneid.v, gene_biotype, gene_testable)]), ...~annot, fill="Not Tested")
astsgenesannotwide[, wheretested := case_when(
  GENCODEv47 == TRUE & PODER == TRUE ~ "Both",
  is.na(GENCODEv47) & PODER == TRUE ~ "Only PODER",
  GENCODEv47 == TRUE & is.na(PODER) ~ "Only GENCODE"
)]
astsgenesannotwide[, `:=`(wheretested=factor(wheretested, levels=c("Only GENCODE", "Both", "Only PODER")), 
                          gene_biotype=factor(gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene", "IG_C_gene")))]


astsgenesannotwideplussamples <- astsgenes[, .(sample, population, map_reads_assemblymap, geneid.v, annot)][astsgenesannotwide, on="geneid.v"]


novelcontaininggenes <- unique(asts[FDR<0.05 & gene_testable==TRUE][, novel_isoform_containing_gene:=any(grepl("transcript|novel", transcriptid.v)), by="geneid.v"][, .(geneid.v, novel_isoform_containing_gene)])
astsgenesannotwideplussamples <- novelcontaininggenes[astsgenesannotwideplussamples, on="geneid.v"]
astsgenesannotwideplussamples[, wheretested_novelty:=factor(fifelse(wheretested=="Only PODER",
                                                                    fifelse(novel_isoform_containing_gene==TRUE, "Only PODER", "Only PODER"), as.character(wheretested)),
                                                            levels=c("Only GENCODE", "Both", "Only PODER"))]
# astsgenesannotwideplussamples[, wheretested_novelty:=factor(fifelse(wheretested=="Only PODER",
#                                                              fifelse(novel_isoform_containing_gene==TRUE, "Only PODER\n(Contains Novel Transcripts)", "Only PODER\n(Only Annotated Transcripts)"), as.character(wheretested)),
#                                                             levels=c("Only GENCODEv47", "Both", "Only PODER\n(Contains Novel Transcripts)"))]

ggplot(unique(astsgenesannotwideplussamples[, .(wheretested_novelty, geneid.v, gene_biotype)]), 
       aes(x=gene_biotype, fill=wheretested_novelty))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)+
  theme(legend.position = c(0.7, 0.9))+
  labs(x="", y="# Significant ASTU Genes", fill="")+
  scale_alpha_manual(values=c(1, 0.7))+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("Novel/Ambiguous Gene"="Novel/Ambiguous\nGene"),
                   guide = guide_axis(n.dodge=2))
ggsave("../10_figures/01_plots/main/fig_04/barplot_sigASTU_geneBiotypes_intersection.pdf", dpi=700, width = 2.25, height = 2.5,  units = "in")
ggsave("../10_figures/01_plots/supp/33_astu_overlap/barplot_sigASTU_geneBiotypes_intersection.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")



# # are we testing (and finding sig) more novel genes than expected by whats in PODER?
# chiinput <-rbind.data.frame(summary(factor(unique(poder[, .(geneid.v, gene_biotype)])$gene_biotype)),
#                             summary(factor(unique(astsgenesannot[annot=="poder", .(geneid.v, gene_biotype)])$gene_biotype)))
# colnames(chiinput) <- c("lncRNA", "novel/ambiguous", "protein_coding")
# rownames(chiinput) <- c("background", "tested")
# chiinput$novellnc <- apply(chiinput[, c("lncRNA", "novel/ambiguous")], 1, sum)
# 
# fisher.test(chiinput[, c("novellnc", "protein_coding")])
# mosaicplot(chiinput[, c("novellnc", "protein_coding")])
# mosaicplot(chiinput[, c("protein_coding", "lncRNA","novel/ambiguous")])

# intersection of tested genes across samples 
astsgenesannotwideplussamples[, samplesharing := uniqueN(sample), by=c("geneid.v","wheretested")]
astsgenesannotwideplussamples[, samplesharingbin :=cut(samplesharing,
                                                       breaks = c(0,  1, 3, 10, Inf),
                                                       labels = c("1", "2-3", "4-10", ">10"),
                                                       right = TRUE), by="wheretested"]
ggplot(unique(astsgenesannotwideplussamples[, .(samplesharingbin, wheretested, geneid.v)]), aes(x=samplesharingbin, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  labs(x="# Samples", y="# Significant ASTU Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)
ggsave("../10_figures/01_plots/supp/33_astu_overlap/barplot_significantASTU_geneBiotypes_intersection.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")




##### NOw for tested only



astsgenes <- unique(asts_test[gene_testable=="TRUE", .(annot, sample, population, map_reads_assemblymap, geneid.v, gene_biotype, gene_testable)])

# which genes are tested with PODER and GENCODE
# astsgenesannot[, tested:="Tested"]

astsgenesannotwide <- dcast(unique(astsgenes[annot%in%c("GENCODEv47", "PODER"), .(annot, geneid.v, gene_biotype, gene_testable)]), ...~annot, fill="Not Tested")
astsgenesannotwide[, wheretested := case_when(
  GENCODEv47 == TRUE & PODER == TRUE ~ "Both",
  is.na(GENCODEv47) & PODER == TRUE ~ "Only PODER",
  GENCODEv47 == TRUE & is.na(PODER) ~ "Only GENCODE"
)]
astsgenesannotwide[, gene_biotype:=fifelse(gene_biotype=="protein_coding", "Protein Coding",
                                           fifelse(gene_biotype=="novel/ambiguous gene", "Novel/Ambiguous\nGene", gene_biotype))]

astsgenesannotwide[, `:=`(wheretested=factor(wheretested, levels=c("Only GENCODE", "Both", "Only PODER")), 
                          gene_biotype=factor(gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous\nGene", "IG_C_gene")))]
astsgenesannotwideplussamples <- astsgenes[, .(sample, population, map_reads_assemblymap, geneid.v, annot)][astsgenesannotwide, on="geneid.v"]


ggplot(unique(astsgenesannotwideplussamples[, .(wheretested, geneid.v, gene_biotype)]), 
       aes(x=gene_biotype, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)+
  theme(legend.position = c(0.7, 0.9))+
  labs(x="", y="# Tested ASTU Genes", fill="")+
  scale_alpha_manual(values=c(1, 0.7))+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("Novel/Ambiguous Gene"="Novel/Ambiguous\nGene"),
                   guide = guide_axis(n.dodge=1))
ggsave("../10_figures/01_plots/supp/33_astu_overlap/barplot_testedASTU_geneBiotypes_intersection.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")



# intersection of tested genes across samples 
astsgenesannotwideplussamples[, samplesharing := uniqueN(sample), by=c("geneid.v","wheretested")]
astsgenesannotwideplussamples[, samplesharingbin :=cut(samplesharing,
                                                       breaks = c(0,  1, 3, 10, Inf),
                                                       labels = c("1", "2-3", "4-10", ">10"),
                                                       right = TRUE), by="wheretested"]
ggplot(unique(astsgenesannotwideplussamples[, .(samplesharingbin, wheretested, geneid.v)]), aes(x=samplesharingbin, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  labs(x="# Samples", y="# Tested ASTU Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)
ggsave("../10_figures/01_plots/supp/33_astu_overlap/barplot_allASTU_geneBiotypes_intersection.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")


### check what % of genes replicate between annotaitons
astu <-unique(asts[FDR<0.05 & annot!="EnhancedGENCODE", .(sample, geneid.v, annot)])
astu <- dcast(astu, sample+geneid.v~annot)
astu[, gencodepoder:=fifelse(GENCODE=="GENCODE" & PODER=="PODER", "both", "one")]
astu[, gencodes:=uniqueN(geneid.v), by="sample"]
astu[gencodepoder=="both", boths:=uniqueN(geneid.v), by="sample"]
astu[, perc:=boths/gencodes*100]
