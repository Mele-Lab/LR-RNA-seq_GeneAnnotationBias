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
ggsave("../../10_figures/01_plots/supp/13_bias_downsamplings/violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC.pdf", dpi=700, width = 5, height = 2.7,  units = "in")


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
ggsave("../../10_figures/01_plots/supp/13_bias_downsamplings/violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC_REFSEQ.pdf", dpi=700, width = 5, height = 2.5,  units = "in")


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



## NOW CHECK THE POPULATION SPECIFIC TRANSCRIPTS REFSEQ
newannotmeta <-newannotmeta[annot=="RefSeq" & !sample%in%c("ITU1", "CEU1", "LWK1", "YRI1", "YRI2", "AJI1", "AJI2", "PEL1", "PEL2", "HAC2", "HAC2")]
newannotmeta[, sample_detected:=1][, sample_per_pop_detected:=sum(sample_detected), by=c("isoform", "population")]
newannotmeta[, i.sample:=NULL][, cell_line_id:=NULL][, nrow:=NULL]
newannotmeta[, populations_detected := uniqueN(population), by = "isoform"]

colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353","#8EDE95")
names(colsqanti) <- unique(newannotmeta$structural_category)[c(1,4,3,2,7,5,8,6)]

newannotmeta[, trx_category:=factor(fifelse(structural_category%in%c("FSM", "ISM"), "known", "novel"), levels=c("known", "novel"))]
newannotmeta[, population:=factor(population, levels=c("MPC", "ITU", "CEU", "LWK", "YRI", "AJI", "PEL", "HAC"))]

newannotmeta[, popsp :=fifelse(populations_detected==1 & sample_per_pop_detected>=2, "popsp", "nonPopsp")]



poptrxlong <-newannotmeta


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
ggsave("../../10_figures/01_plots/supp/16_popsp_validations/barplot_PODER_popSpecificTrx_trx.bycategory_DOWNSAMPLING_REFSEQ.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")





## NOW CHECK THE POPULATION SPECIFIC TRANSCRIPTS GENCODE
newannotmeta <-newannotmeta[annot=="GENCODE" & !sample%in%c("ITU1", "CEU1", "LWK1", "YRI1", "YRI2", "AJI1", "AJI2", "PEL1", "PEL2", "HAC2", "HAC2")]
newannotmeta[, sample_detected:=1][, sample_per_pop_detected:=sum(sample_detected), by=c("isoform", "population")]
newannotmeta[, i.sample:=NULL][, cell_line_id:=NULL][, nrow:=NULL]
newannotmeta[, populations_detected := uniqueN(population), by = "isoform"]

colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353","#8EDE95")
names(colsqanti) <- unique(newannotmeta$structural_category)[c(1,4,3,2,7,5,8,6)]

newannotmeta[, trx_category:=factor(fifelse(structural_category%in%c("FSM", "ISM"), "known", "novel"), levels=c("known", "novel"))]
newannotmeta[, population:=factor(population, levels=c("MPC", "ITU", "CEU", "LWK", "YRI", "AJI", "PEL", "HAC"))]

newannotmeta[, popsp :=fifelse(populations_detected==1 & sample_per_pop_detected>=2, "popsp", "nonPopsp")]



poptrxlong <-newannotmeta


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
ggsave("../../10_figures/01_plots/supp/16_popsp_validations/barplot_PODER_popSpecificTrx_trx.bycategory_DOWNSAMPLING_GENCODE.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")













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




