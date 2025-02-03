## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling/"
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
library(ggpubr)
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

# load metadata
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE]
metadata <- metadata[mixed_samples==FALSE]
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C")
names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC")

# load downsampling results
sqanti_refseq <- fread("../04_evaluation/02_sqanti/data/downsampling_espresso_refseq/downsamplig_espresso_refseq_classification.txt")
annot_refseq <- fread("../../../novelannotations/merge_espresso_downsample_refseq/merged.gtf") # merged annotation
sqanti_gencode <- fread("../04_evaluation/02_sqanti/data/downsampling_espresso_gencode/downsampling_espresso_classification.txt")
annot_gencode <- fread("../../../novelannotations/merge_espresso_downsample/merged_parsednames.gtf")


#### PARSE ANNOTATION
colnames(annot_gencode) <- c("contig", "platform", "feature", "start", "end", "unk1", "strand", "unk2", "info")
colnames(annot_refseq) <- c("contig", "platform", "feature", "start", "end", "unk1", "strand", "unk2", "info")

annot_gencode <- annot_gencode[feature=="transcript",]
annot_refseq <- annot_refseq[feature=="transcript",]

annot_gencode[, transcript := tstrsplit(info, "\"")[[2]]]
annot_gencode[, samples := tstrsplit(info, "\"")[[4]]]
annot_refseq[, transcript := tstrsplit(info, "\"")[[2]]]
annot_refseq[, samples := tstrsplit(info, "\"")[[4]]]

# Split the tool_sample_pairs into individual pairs
annot_gencode[, samples_split := strsplit(samples, ",")]
annot_refseq[, samples_split := strsplit(samples, ",")]


annot_gencode_long <- annot_gencode[, .(sample = unlist(strsplit(samples, ","))), by = transcript][, annot:="GENCODE"]
annot_refseq_long <- annot_refseq[, .(sample = unlist(strsplit(samples, ","))), by = transcript][, annot:="RefSeq"][, sample:=gsub("espresso_q_", "", sample)]

# add sqanti to annot_expanded to trace the origin of the transcripts
newannot_expanded_gencode <- sqanti_gencode[, .(isoform, structural_category)][annot_gencode_long, on=c("isoform"="transcript")]
newannot_expanded_refseq <- sqanti_refseq[, .(isoform, structural_category)][annot_refseq_long, on=c("isoform"="transcript")]


newannot <- rbind.data.frame(newannot_expanded_gencode, newannot_expanded_refseq)
newannot[, cell_line_id := tstrsplit(sample, "_")[[3]]]
newannotmeta <- metadata[, .(cell_line_id, population, map_reads_assemblymap, sample)][newannot, on="cell_line_id"]
newannotmeta[, trx_per_sample:=uniqueN(isoform), by=c("sample", "annot")]
newannotmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]
newannotmeta[, afr:=ifelse(population%in%c("MPC", "YRI", "LWK"), "African", "OOA")]
newannotmeta[, mean:=mean(trx_per_sample), by=c("annot", "population")]
newannotmeta[, sd:=sd(trx_per_sample), by=c("annot", "population")]

# Create the ggplot barplot
ggplot(unique(newannotmeta[, .(population, annot, mean, sd)]), aes(x = population, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, color = "#545454") + # Error bars
  geom_bar(stat = "identity", aes(fill = population), alpha=0.75) + # Bar plot
  labs(y = "# Discovered Transcripts", x="") +
  mytheme+
  scale_fill_manual(values=popcols)+
  ggbeeswarm::geom_quasirandom(data=unique(newannotmeta[, .(population, annot, sample, trx_per_sample)]),aes(x=population, y=trx_per_sample), alpha=0.5, size=0.5)+
  guides(fill="none")+
  facet_wrap(~annot)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave("../../10_figures/01_plots/main/fig_01/barplot_DiscoveredTranscripts_PerPop.pdf", dpi=700, width = 2.3, height = 2.4,  units = "in")

## test if there are differences-------------------------

# Test for Normality
shapiro.test(unique(newannotmeta[annot=="GENCODE", .(trx_per_sample, population, eur)])$trx_per_sample) # Apply separately for each group if needed
shapiro.test(unique(newannotmeta[annot=="RefSeq", .(trx_per_sample, population, eur)])$trx_per_sample) # Apply separately for each group if needed

# Test for Homogeneity of Variances
car::leveneTest(trx_per_sample ~ population, data = unique(newannotmeta[annot=="GENCODE", .(trx_per_sample, population, eur)]))
car::leveneTest(trx_per_sample ~ population, data = unique(newannotmeta[annot=="RefSeq", .(trx_per_sample, population, eur)]))

# Perform ANOVA
anova_result <- aov(trx_per_sample ~ eur, data = unique(newannotmeta[annot=="GENCODE", .(trx_per_sample, population, eur, annot)]))
anova_result <- aov(trx_per_sample ~ eur, data = unique(newannotmeta[annot=="RefSeq", .(trx_per_sample, population, eur, annot)]))
anova_result <- aov(trx_per_sample ~ population, data = unique(newannotmeta[annot=="GENCODE", .(trx_per_sample, population, eur, annot)]))
anova_result <- aov(trx_per_sample ~ population, data = unique(newannotmeta[annot=="RefSeq", .(trx_per_sample, population, eur, annot)]))
summary(anova_result)

# Post-hoc Tukey Test
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

#####-------------------------------------------------------


categories <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")
names(categories) <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "intergenic", "genic", "fusion", "antisense")
newannotmeta$structural_category <- categories[newannotmeta$structural_category]
newannotmeta[, structural_category :=factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense"))]
p <- ggplot(unique(newannotmeta[annot=="GENCODE" & structural_category%in%c("FSM", "ISM", "NIC", "NNC")][, trx_per_sample_cat := uniqueN(isoform), by=c("sample", "structural_category")][, .(trx_per_sample_cat, structural_category,population, eur)]), 
       aes(x = eur, y = trx_per_sample_cat, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.9)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.6, width = 0.2, dodge.width = 0.8, aes(col=population)) +
  mytheme +
  scale_fill_manual(values = c("#466995", "#A53860")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1000), 
               position = position_dodge(0.9), size=7*0.35) +
  labs(x = "", y = "# Discovered Transcripts", col="Population")+
  facet_wrap(~structural_category, nrow=1)+
  guides(fill="none")+
  scale_color_manual(values=popcols)+
  theme(legend.key.size = unit(0.2, "cm"), 
        axis.text.x=element_text(angle=45,vjust=0.75, hjust=0.5),
        legend.position = c(0.1, 0.4))+
  geom_pwc(ref.group="European",
           method="t_test", label.size=7*0.35
  )+ylim(c(-1000,31500))+
  scale_x_discrete(labels = c("European", "Non\nEuropean"))
ggadjust_pvalue(
  p=p,
  p.adjust.method = "BH",
  label = "p.adj.format",
  hide.ns = NULL,
  output = c("plot"))
ggsave("../../10_figures/01_plots/main/fig_01/violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC.pdf", dpi=700, width = 5, height = 2.7,  units = "in")


p <- ggplot(unique(newannotmeta[annot=="RefSeq" & structural_category%in%c("FSM", "ISM", "NIC", "NNC")][, trx_per_sample_cat := uniqueN(isoform), by=c("sample", "structural_category")][, .(trx_per_sample_cat, structural_category,population, eur)]), 
            aes(x = eur, y = trx_per_sample_cat, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.9)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.6, width = 0.2, dodge.width = 0.8, aes(col=population)) +
  mytheme +
  scale_fill_manual(values = c("#466995", "#A53860")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1000), 
               position = position_dodge(0.9), size=7*0.35) +
  labs(x = "", y = "# Discovered Transcripts", col="Population")+
  facet_wrap(~structural_category, nrow=1)+
  guides(fill="none")+
  scale_color_manual(values=popcols)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.position = c(0.1, 0.4))+
  geom_pwc(ref.group="European",
           method="t_test", label.size=7*0.35
  )+ylim(c(-1500,20000))+
  scale_x_discrete(labels = c("European", "Non\nEuropean"))
ggadjust_pvalue(
  p=p,
  p.adjust.method = "BH",
  label = "p.adj.format",
  hide.ns = NULL,
  output = c("plot"))
ggsave("../../10_figures/01_plots/main/fig_03/violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC_REFSEQ.pdf", dpi=700, width = 5, height = 2.5,  units = "in")




p<-ggplot(unique(newannotmeta[annot=="GENCODE" & structural_category%in%c("FSM", "ISM", "NIC", "NNC")][, trx_per_sample_cat := uniqueN(isoform), by=c("sample", "structural_category")][
    , .(trx_per_sample_cat, structural_category,population, afr)]), 
            aes(x = afr, y = trx_per_sample_cat, fill = afr)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.9)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.6, width = 0.2, dodge.width = 0.8, aes(col=population)) +
  mytheme +
  # ggpubr::stat_compare_means(aes(label=..p.adj..), comparisons=list(c("European", "non-European")),
  #                            method = "t.test", 
  #                            method.args = list(alternative = "two.sided"),
  #                            p.adjust.method ="fdr") +
  scale_fill_manual(values =  c("#F7D257", "#496F5D")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = 0), 
               position = position_dodge(0.9), size=7*0.35) +
  labs(x = "", y = "# Discovered Transcripts", col="Population")+
  facet_wrap(~structural_category, nrow=2)+
  guides(fill="none")+
  scale_color_manual(values=popcols)+
  theme(legend.key.size = unit(0.2, "cm"))+
  geom_pwc(ref.group="African",
           method="t_test", label.size=7*0.35
  )+ylim(c(0,31000))
ggadjust_pvalue(
  p=p,
  p.adjust.method = "BH",
  label = "p.adj.format",
  hide.ns = NULL,
  output = c("plot"))
ggsave("../../10_figures/01_plots/supp/03_afr_down_cat/01_violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC.pdf", dpi=700, width = 6, height = 6,  units = "in")
ggsave("../../10_figures/02_panels/png/supp/03_afr_down_cat.png", dpi=700, width = 6, height = 6,  units = "in")
ggsave("../../10_figures/02_panels/svg/supp/03_afr_down_cat.svg", dpi=700, width = 6, height = 6,  units = "in")





## NOW CHECK THE POPULATION SPECIFIC TRANSCRIPTS
newannotmeta <-newannotmeta[annot=="RefSeq" & !sample%in%c("ITU1", "CEU1", "LWK1", "YRI1", "YRI2", "AJI1", "AJI2", "PEL1", "PEL2", "HAC2", "HAC2")]
newannotmeta[, sample_detected:=1][, sample_per_pop_detected:=sum(sample_detected), by=c("isoform", "population")]
newannotmeta[, i.sample:=NULL][, cell_line_id:=NULL][, nrow:=NULL]
newannotmeta[, populations_detected := uniqueN(population), by = "isoform"]

colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353","#8EDE95")
names(colsqanti) <- unique(newannotmeta$structural_category)[c(1,4,3,2,7,5,8,6)]

newannotmeta[, novelty:=factor(fifelse(structural_category%in%c("FSM", "ISM"), "Annotated", "Novel"), levels=rev(c("Annotated", "Novel")))]
newannotmeta[, population:=factor(population, levels=c("MPC", "ITU", "CEU", "LWK", "YRI", "AJI", "PEL", "HAC"))]
# ggplot(newannotmeta[populations_detected==1 & sample_per_pop_detected>=2], aes(x=population, fill=structural_category))+
#   geom_bar(position="fill")+
#   mytheme+
#   scale_fill_manual(values=colsqanti)+
#   geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))
# ggplot(newannotmeta[populations_detected==1 & sample_per_pop_detected>=2 & structural_category%in%c("FSM", "ISM", "NIC", "NNC")], aes(x=population, fill=structural_category))+
#   geom_bar(position="fill")+
#   mytheme+
#   scale_fill_manual(values=colsqanti)+
#   geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))+
#   labs(x="", y="Proportion of Population Specific Transcripts", fill="")
# ggsave("figures/.pdf", dpi=700, width = 25, height = 14,  units = "cm")


popsp <-newannotmeta[populations_detected==1 & sample_per_pop_detected>=2 & structural_category%in%c("FSM", "ISM", "NIC", "NNC")]
isoform_counts <- popsp[novelty == "Annotated", uniqueN(isoform), by = population]

# Step 2: Reorder `population` based on the counts calculated in Step 1
# (Higher counts appear first)
popsp$population <- factor(popsp$population, 
                                  levels = isoform_counts[order(-V1)]$population)
popsp[, count:=uniqueN(isoform), by=c("population", "novelty")]
popsp[, knowncount := ifelse(any(novelty == "Annotated"), 
                             count[novelty == "Annotated"][1], 
                             NA), 
      by = population]
popsp <- unique(popsp[, .(novelty, population, count,knowncount)])
popsp[, perc:=paste0(round(count/sum(count)*100, digits=0), "%"), by="population"]
ggplot(unique(popsp[, .(perc,novelty, population, count,knowncount)]), 
       aes(x=reorder(population,-knowncount ),y=count,  fill=novelty))+
  geom_col()+
  mytheme+
  scale_fill_manual(values=rev(c("#5E8531","#A2331D")))+
  geom_text(aes(label=perc), position=position_stack(vjust=0.5), size=6*0.35)+
  labs(x="", y="# Population Specific Transcripts", fill="Intron Chain")+
  theme(legend.position="top")
ggsave("../../10_figures/01_plots/supp/18_popsp_val_refseq/barplotStack.downsampling_annotatedVSnovelIntronChain_POPspecificTrx.pdf", dpi=700, width = 3, height = 3,  units = "in")
# ggsave("../../10_figures/01_plots/supp/16_popsp_val_gencode/barplotStack.downsampling_annotatedVSnovelIntronChain_POPspecificTrx.pdf", dpi=700, width = 3, height = 3,  units = "in")



##### Check if any population has more transcripts than others
# select 4 random samples per pop
selected_samples <- unique(newannotmeta[, .(population, sample)])[, .SD[sample(.N, min(.N, 4))], by = population]$sample

foursample <- newannotmeta[annot=="GENCODE" & sample%in%selected_samples]
foursample <- unique(foursample[, .(population, isoform, structural_category, eur, afr)])[, nb_trx_pop:=uniqueN(isoform), by="population"]
foursamplep <- unique(foursample[, .(eur, afr, population, nb_trx_pop)])

ggplot(foursamplep, aes(x=population, y=nb_trx_pop))+
  geom_col()+
  mytheme
ggplot(foursamplep, aes(x=afr, y=nb_trx_pop))+
  ggbeeswarm::geom_quasirandom(aes(color=population))+
  mythemen+
  scale_color_manual(values=popcols)+
  ggpubr::stat_compare_means(comparisons = list(c("African", "OOA")), method="wilcox", method.args = list(alternative="two.sided"), label.y=73000)+
  labs(x="", y="# Unique Transcripts Detected in Population")
