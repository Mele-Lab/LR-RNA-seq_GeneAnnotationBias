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
                              .(structural_category, trx_per_cat_norm, eur, sample,population)]), 
          aes(x = eur, y = trx_per_cat_norm, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.8) + # Adjust dodge width if necessary
  geom_boxplot(outliers=F, width=0.05,position = position_dodge(0.9))+
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.1,  dodge.width = 0.6,aes(color=population))+
  mytheme+
  facet_wrap(~structural_category, nrow=2)+
  labs(x="", y="# Discovered Transcripts in UMA\nper Million Reads per Sample", color="Population")+
  geom_pwc(ref.group = "European", method = "t_test", label.size = 7*0.35) +
  scale_fill_manual(values=c("#466995", "#A53860"))+
  scale_color_manual(values=popcols)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5, size=6*0.35)+
  guides(fill="none", color="none")+
  ylim(c(0, 4250))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggadjust_pvalue(p)
ggsave("../../10_figures/01_plots/main/fig_03/violin_UMA_discovered_transcripts.EURnonEUR_PerSqantiCategory.pdf", dpi=700, width = 2.5, height = 5,  units = "in")



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
ggsave("../../10_figures/01_plots/main/fig_03/violin_PODER_discovered_transcripts.EURnonEUR_PerSqantiCategory.pdf", dpi=700, width = 4.25, height = 2.5,  units = "in")

datalongmeta[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "Out of African")]
p<-ggplot(unique(datalongmeta[structural_category %in% c("FSM","NIC", "NNC") & filter=="pass", 
                              .(structural_category, trx_per_cat_norm_poder, afr, sample,population, filter)]), 
          aes(x = afr, y = trx_per_cat_norm_poder, fill = afr)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + # Adjust dodge width if necessary
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.1,  dodge.width = 0.6,aes(color=population), size=1)+
  geom_boxplot(outliers=F, width=0.3,position = position_dodge(0.9),lwd=0.2)+
  mytheme+
  facet_wrap(~structural_category)+
  labs(x="", y="# Discovered Transcripts in PODER\nper Million Reads per Sample", color="Population")+
  geom_pwc(ref.group = "African", method = "t_test", label.size=6*0.35) +
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  scale_color_manual(values=popcols)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5, size=6*0.35)+
  guides(fill="none")+
  ylim(c(0, 4000))+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1.1))
ggadjust_pvalue(p)
ggsave("../../10_figures/01_plots/main/fig_03/violin_PODER_discovered_transcripts.AFRooA_PerSqantiCategory.pdf", dpi=700, width = 4.25, height = 2.5,  units = "in")




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




popsp <-data[, popsp :=fifelse(population_sharing==1 & sample_sharing>=2, "popsp", "nonPopsp") ]
cols_to_keep <- colnames(popsp)[grepl("^[A-Z]{3}$", colnames(popsp))]
subpopsp <- popsp[, c(cols_to_keep,  "structural_category", "isoform", "popsp"), with = FALSE]
subpopsplong <- melt(subpopsp, measure.vars= colnames(subpopsp)[grepl("^[a-zA-Z]{3}$", colnames(subpopsp))], variable.name = "population", value.name = "detected")
subpopsplong <- subpopsplong[detected!=0]


subpopsplong[, trx_category:=factor(fifelse(structural_category%in%c("FSM", "ISM"), "known", "novel"), levels=c("known", "novel"))]



poptrxlong <-subpopsplong


### COMPUTE ODDS RATIO
counts <- poptrxlong[, .N, by=.(trx_category, population,popsp)]

# Run Fisher's exact test by population
results <- lapply(unlist(unique(counts[,population])), function(pop) {
  # Filter data for the population
  pop_data <- counts[population == pop]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[popsp == "nonPopsp" & trx_category == "known", N],
      pop_data[popsp == "nonPopsp" & trx_category == "novel", N],
      pop_data[popsp == "popsp" & trx_category == "known", N],
      pop_data[popsp == "popsp" & trx_category == "novel", N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("known", "novel"),
      rev(c("popsp", "nonPopsp"))
    )
  )
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, conf.int = T)
  
  list(population = pop, p.value = test_result$p.value, odds.ratio = test_result$estimate, ci=list(test_result$conf.int))
})

resultss <- rbindlist(results)
# Extract lower and upper CI values into new columns
resultss[, `:=`(ci_lower = sapply(ci, `[`, 1),  # Extract first value of CI
                ci_upper = sapply(ci, `[`, 2))] 

resultss[, eur:=fifelse(population%in%c("AJI", "CEU"), "EUR", "nonEUR")]
resultss[, fdr:=p.adjust(p.value, method = "BH")]
resultss[, FDR:=factor(fifelse(fdr<0.05, "FDR<0.05", "FDR>=0.05"))]
resultss[, population:=factor(population, levels=levels(poptrxlong$population))]

order_samples <-counts[popsp=="popsp", total:=sum(N), by=.(population)][, known_per:=N/total][trx_category=="known" & popsp=="popsp"][order(known_per), population]

p3 <- ggplot(resultss, aes(x=odds.ratio, y=population, color=eur))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5, show.legend=F)+
  geom_vline(xintercept=1, linetype="dashed")+
  annotate(geom="rect", 
           xmin=min(resultss[eur=="nonEUR", ci_lower]),
           xmax=max(resultss[eur=="nonEUR", ci_upper]),
           ymin=0.5,
           ymax=8.5,
           fill="#A53860",
           alpha=0.2)+
  annotate(geom="rect", 
           xmin=min(resultss[eur=="EUR", ci_lower]),
           xmax=max(resultss[eur=="EUR", ci_upper]),
           ymin=0.5,
           ymax=8.5,
           fill="#466995",
           alpha=0.2)+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_alpha_manual(values=rev(c(0.45, 1)))+
  mytheme+
  labs(x="Odds Ratio\n(log10)", y="", color="", alpha="")+
  theme(legend.position = "top")+
  scale_x_continuous(trans="log10", breaks=c(0.5, 1, 2, 3))+
  scale_y_discrete(limits=order_samples)+
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha=guide_legend(override.aes = list(size = 3), nrow=2))+
  annotation_logticks(sides="b")+
  guides(color="none")

dt1 <- unique(counts[popsp=="popsp", total:=sum(N), by=.(population)][, known_per:=N/total][, eur:=ifelse(population%in%c("CEU", "AJI"), "European", "non-European")][popsp=="popsp"])

p1 <- ggplot(dt1, 
             aes(x=population, fill=eur, y=known_per))+
  geom_col(position="fill", aes(alpha=trx_category))+
  geom_text(data=dt1,
            aes(label = paste0(round(known_per*100, 1), "%"), group=trx_category),
            position = position_fill(vjust = 0.5),
            size=6*0.35 )+
  scale_x_discrete(limits=order_samples)+
  labs(x="", y="Proportion of\nPopulation-Specific\nTranscripts", fill="")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  scale_alpha_manual(values=c( 0.25,0.75),name="",labels = c("known"="Known", "novel"="Novel"))+
  scale_color_manual(values=rep("black", 2))+
  guides(color="none")+
  coord_flip()+
  mytheme+
  theme(legend.position = "top",
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(margin = margin(r = -1, unit = "pt")),
        legend.box.margin=margin(-10,-10,-10,-45),
        legend.box = "horizontal",               # Ensure horizontal alignment
        legend.box.just = "center",              # Center the legend box
        legend.box.spacing = unit(0.5, "cm"),    # Adjust spacing between rows
        legend.direction = "horizontal",         # Horizontal direction for multiple rows
        legend.justification = "left",
        legend.spacing.y = unit(0, "cm"),
        plot.margin = margin(0, 3, 5, -5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 3),order = 2,nrow = 2),
    alpha=guide_legend(order = 1, nrow = 2, override.aes = list(size = 3))
  )+
  scale_y_continuous(expand =  c(0, 0, 0, 0))

p2 <- ggplot(unique(poptrxlong[popsp=="popsp", .(isoform,trx_category, population)])[, eur:=fifelse(population%in%c("AJI", "CEU"), "European", "Non-European")], aes(x=population, fill=eur))+
  geom_bar()+
  mytheme+
  coord_flip()+
  xlab("")+
  guides(fill="none")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  ylab("# Transcripts")+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),  plot.margin = margin(0, 0, 5, -10),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  scale_y_continuous(expand =  c(0, 0, 0.1, 0), n.breaks=3)+
  scale_x_discrete(limits=order_samples)

library(patchwork)
# Ensure both p1 and p2 have the same theme
p1 <- p1 + theme(legend.position = "top")

# Combine the plots using patchwork
p3+p1+p2 +
  plot_layout(guides="collect",  widths=c(3,5.5, 0.75))&
  theme(legend.position='top')
ggsave("../../10_figures/01_plots/supp/16_popsp_validations/barplot_PODER_popSpecificTrx_trx.bycategory_UMA.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")

