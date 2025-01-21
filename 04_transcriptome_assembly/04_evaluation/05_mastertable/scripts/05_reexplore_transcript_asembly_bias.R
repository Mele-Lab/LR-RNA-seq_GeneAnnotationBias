## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation"
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



data <- fread("05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv")
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,2,3,5,4,8,6,7)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
data[, filter:=factor(filter, levels=c("pass", "fail"))]

# Identify columns that match the pattern (3 letters + 1 number)
cols_to_keep <- colnames(data)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(data))]
subdata <- data[, c(cols_to_keep, "filter", "sample_sharing", "structural_category", "isoform"), with = FALSE]
datalong <- melt(subdata, measure.vars= colnames(subdata)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(subdata))], variable.name = "sample", value.name = "detected")
datalong <- datalong[!is.na(detected)]
datalong <- datalong[detected==1]

datalongmeta <- metadata[, .(sample, map_reads_assemblymap,population )][datalong, on="sample"]

datalongmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "EUR", "nonEUR")]
datalongmeta[, trx_per_cat:=uniqueN(isoform), by=c("sample", "structural_category")]
datalongmeta[, trx_per_cat_norm :=trx_per_cat/map_reads_assemblymap*10^6]
datalongmeta[filter=="pass", trx_per_cat_poder:=uniqueN(isoform), by=c("sample", "structural_category")]
datalongmeta[filter=="pass", trx_per_cat_norm_poder :=trx_per_cat_poder/map_reads_assemblymap*10^6]

# PLOT
library(ggbeeswarm)
library(ggpubr)
datalongmeta[, eur:=fifelse(eur=="EUR", "European", "Non-European")]

### UMA
p<-ggplot(unique(datalongmeta[structural_category %in% c("FSM", "ISM","NIC", "NNC"), 
                              .(structural_category, trx_per_cat_norm, eur, sample,population, filter)]), 
          aes(x = eur, y = trx_per_cat_norm, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.8) + # Adjust dodge width if necessary
  geom_boxplot(outliers=F, width=0.05,position = position_dodge(0.9))+
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.1,  dodge.width = 0.6,aes(color=population))+
  mytheme+
  facet_wrap(~structural_category)+
  labs(x="", y="# Discovered Transcripts in UMA\nper Million Reads per Sample", color="Population")+
  geom_pwc(ref.group = "European", method = "t_test") +
  scale_fill_manual(values=c("#466995", "#A53860"))+
  scale_color_manual(values=popcols)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)+
  guides(fill="none")+
  ylim(c(0, 4250))
ggadjust_pvalue(p)
ggsave("../../10_figures/01_plots/supp/13_uma_biasval_cat/violin_UMA_discovered_transcripts.EURnonEUR_PerSqantiCategory.pdf", dpi=700, width = 6, height = 6,  units = "in")
ggsave("../../10_figures/02_panels/svg/supp/13_uma_biasval_cat.svg", dpi=700, width = 6, height = 6,  units = "in")
ggsave("../../10_figures/02_panels/png/supp/13_uma_biasval_cat.png", dpi=700, width = 6, height = 6,  units = "in")


## PODER
p<-ggplot(unique(datalongmeta[structural_category %in% c("FSM","NIC", "NNC") & filter=="pass", 
                              .(structural_category, trx_per_cat_norm_poder, eur, sample,population, filter)]), 
          aes(x = eur, y = trx_per_cat_norm_poder, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + # Adjust dodge width if necessary
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.1,  dodge.width = 0.6,aes(color=population), size=1)+
  geom_boxplot(outliers=F, width=0.3,position = position_dodge(0.9),lwd=0.2)+
  mytheme+
  facet_wrap(~structural_category)+
  labs(x="", y="# Discovered Transcripts in PODER\nper Million Reads per Sample", color="Population")+
  geom_pwc(ref.group = "European", method = "t_test", label.size=6*0.35) +
  scale_fill_manual(values=c("#466995", "#A53860"))+
  scale_color_manual(values=popcols)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5, size=6*0.35)+
  guides(fill="none")+
  ylim(c(0, 4000))+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1.1))
ggadjust_pvalue(p)
ggsave("../../10_figures/01_plots/main/fig_03//violin_PODER_discovered_transcripts.EURnonEUR_PerSqantiCategory.pdf", dpi=700, width = 4.25, height = 2.5,  units = "in")




ggplot(unique(datalongmeta[structural_category %in% c("FSM", "ISM"), 
                           .(structural_category, trx_per_cat_norm_poder, eur, sample,population)]), 
       aes(x = eur, y = trx_per_cat_norm_poder, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.8) + # Adjust dodge width if necessary
  geom_boxplot(outliers=F, width=0.05,position = position_dodge(0.9))+
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.2,  dodge.width = 0.9,)+
  mytheme+
  facet_wrap(~structural_category)+
  labs(x="", y="# Assembled Transcripts/Million Reads per Sample", fill="")+
  ggpubr::stat_compare_means( comparisons=list(c("EUR", "nonEUR")),method = "t.test", method.args = list(alternative="greater"))+
  scale_fill_manual(values=c("#466995", "#A53860"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)
datalongmeta[, eur := as.factor(eur)]
ggplot(unique(datalongmeta[structural_category %in% c("NIC", "NNC"), 
                           .(structural_category, trx_per_cat_norm_poder, eur, sample,population)]), 
       aes(x = eur, y = trx_per_cat_norm_poder, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.8) + # Adjust dodge width if necessary
  geom_boxplot(outliers=F, width=0.05,position = position_dodge(0.9))+
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.2,  dodge.width = 0.9,)+
  mytheme+
  facet_wrap(~structural_category)+
  labs(x="", y="# Assembled Transcripts/Million Reads per Sample", fill="")+
  ggpubr::stat_compare_means( comparisons=list(c("EUR", "nonEUR")),method = "t.test", method.args = list(alternative="less"))+
  scale_fill_manual(values=c("#466995", "#A53860"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=200), vjust=0.5)




popsp <-data[sample_sharing>=2 & population_sharing==1, ]
cols_to_keep <- colnames(popsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(popsp))]
subpopsp <- popsp[, c(cols_to_keep, "filter", "sample_sharing", "structural_category", "isoform", "trx_per_asstrx_count_incomplete-splice_match"), with = FALSE]
subpopsplong <- melt(subpopsp, measure.vars= colnames(subpopsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(subpopsp))], variable.name = "sample", value.name = "detected")
subpopsplong <- subpopsplong[!is.na(detected)]
subpopsplong <- subpopsplong[detected==1]

subpopsplongmeta <- metadata[, .(sample, map_reads_assemblymap,population )][, total_throughput:=sum(map_reads_assemblymap), by="population"][subpopsplong, on="sample"]
subpopsplongmeta[, trx_per_cat_per_pop := uniqueN(isoform), by=c("population", "structural_category")]
subpopsplongmeta[, trx_per_cat_per_pop_filter := uniqueN(isoform), by=c("population", "structural_category", "filter")]
subpopsplongmeta[, trx_per_cat_per_pop_norm :=trx_per_cat_per_pop/total_throughput*10^6]
subpopsplongmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]

# LINE PLOT FOR 2 SAMPLES SHARING
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific UMA Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  geom_segment(
    aes(xend = total_throughput/10^6, # Adjust for arrow length
        yend = ifelse(eur == "European", trx_per_cat_per_pop + 50, trx_per_cat_per_pop - 50)),
    arrow = arrow(length = unit(0.2, "cm")), # Arrow size
    color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(
    aes(x = total_throughput/10^6,  # Position text away from point
        y = ifelse(eur == "European", trx_per_cat_per_pop + 75, trx_per_cat_per_pop - 75),
        label = population,
      color=population),
    fontface="bold",
    family="Helvetica",
    size = 4.5,
    alpha=0.6)+
  scale_color_manual(values=popcols)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../10_figures/suppfig/line_UMA_popSpecific_transcripts.2samplesSharing.faceted.pdf", dpi=700, width = 24, height = 20,  units = "cm")

# LINE PLOT FOR 3! SAMPLES SHARING
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific UMA Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  geom_segment(
    aes(xend = total_throughput/10^6, # Adjust for arrow length
        yend = ifelse(eur == "European", trx_per_cat_per_pop + 5, trx_per_cat_per_pop - 5)),
    arrow = arrow(length = unit(0.2, "cm")), # Arrow size
    color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(
    aes(x = total_throughput/10^6,  # Position text away from point
        y = ifelse(eur == "European", trx_per_cat_per_pop + 7.5, trx_per_cat_per_pop - 7.5),
        label = population,
        color=population),
    fontface="bold",
    family="Helvetica",
    size = 4.5,
    alpha=0.6)+
  scale_color_manual(values=popcols)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../10_figures/suppfig/line_UMA_popSpecific_transcripts.3samplesSharing.faceted.pdf", dpi=700, width = 24, height = 20,  units = "cm")
