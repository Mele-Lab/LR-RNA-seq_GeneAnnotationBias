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
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- metadata$color_pop
names(popcol) <- metadata$population





data <- fread("../novelannotations/analysis_tables/250206_personalized_hg38s_n_ic.tsv")


ggplot(data, aes(x=map_genome, y=n_ic))+
  geom_point()+
  mytheme+
  geom_line(aes(color=cell_line_id, group=cell_line_id))+
  labs(x="Genome Assembly", y="# Intron Chains", color="Cell Line")+
  scale_x_discrete(labels=c("hap1"="Haplotype 1\n(personalized\nGRCh38 1)","hap2"="Haplotype 2\n(personalized\nGRCh38 2)", "hg38"="GRCh38" ))
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/lineplot_numberIntronChains_perHaplotype.pdf", dpi=700, width = 3.75, height = 2.75,  units = "in")


data <- fread("../novelannotations/analysis_tables/250206_personalized_hg38s_n_sj.tsv")
ggplot(data, aes(x=map_genome, y=n_sj))+
  geom_point()+
  mytheme+
  geom_line(aes(color=cell_line_id, group=cell_line_id))+
  labs(x="Genome Assembly", y="# Splice Junctions", color="Cell Line")+
  scale_x_discrete(labels=c("hap1"="Haplotype 1\n(personalized\nGRCh38 1)","hap2"="Haplotype 2\n(personalized\nGRCh38 2)", "hg38"="GRCh38" ))
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/lineplot_numberSpliceJunctions_perHaplotype.pdf", dpi=700, width = 3.75, height = 2.75,  units = "in")


# load data
# data <- fread("../novelannotations/analysis_tables/250206_perc_novel_ics_uniq_to_pers_hg38.tsv")

data <- fread("../novelannotations/analysis_tables/250226_ratio_novel_ics_in_haps_vs_hg38.tsv")
data <- data[, .(cell_line_id, n_novel_ic, n_non_hg38_novel_ic, n_novel_ic_hg38, ratio_novel, diff, population, map_reads_assemblymap)]

data[, increase:=(ratio_novel-1)*100]
data[, mean:=mean(increase), by="population"]
data[, sd:=sd(increase), by="population"]

ggplot(unique(data[, .(mean, sd, population)]), aes(x=reorder(population, mean), fill=population))+
  geom_col(aes(y=mean), alpha=0.8)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, linewidth=0.2)+
  ggbeeswarm::geom_quasirandom(data=data, mapping=aes(y=increase),size=0.3)+
  mytheme+
  labs(x="", y="% increase in Novel Transcripts\nUnique to personalized-GRCh38s")+
  scale_fill_manual(values=popcol)+
  guides(fill="none")+
  annotate(geom="text", x=1.5, y=4.5, label="p=0.0004", size=6*0.35)
ggsave("10_figures/01_plots/main/personalizedhg38/barplot_IntronChainsUniquetopersonalizedGRCh38s.pdf", dpi=700, width = 1.8, height = 2.25,  units = "in")

## test if there are differences-------------------------

# Test for Normality
shapiro.test(unique(data[, .(increase, population)])$increase) # Apply separately for each group if needed

# Test for Homogeneity of Variances
car::leveneTest(increase ~ population, data = unique(data[, .(increase, population)]))

# Perform ANOVA
anova_result <- aov(increase ~ population, data = unique(data[, .(increase, population)]))
anovares <-summary(anova_result)
anovares[[1]]$`Pr(>F)`[1]
# Post-hoc Tukey Test
tukey_result <- TukeyHSD(anova_result)
tukey_result$population
# Split the comparison into two groups
mydata <- as.data.frame(tukey_result$population)
mydata <- rownames_to_column(mydata, var = "comparison")
tukey_result <- separate(mydata, comparison, into = c("Group1", "Group2"), sep = "-")
colnames(tukey_result)[6] <- "padj"

setDT(tukey_result)
tukey_result[, sig:=ifelse(padj<0.05, "*", "")]
tukey_result <- tukey_result[, .(Group1, Group2, diff, padj, sig)]
tukey_result <- rbind.data.frame(tukey_result, c(rep("CEU", 2), NA,NA, NA))
tukey_result <- rbind.data.frame(tukey_result, c(rep("HAC", 2), NA,NA, NA))
tukey_result <- rbind.data.frame(tukey_result, c(rep("ITU", 2), NA, NA,NA))
tukey_result <- rbind.data.frame(tukey_result, c(rep("YRI", 2), NA, NA, NA))
tukey_result <- rbind.data.frame(tukey_result, c(rep("LWK", 2), NA, NA, NA))
tukey_result <- rbind.data.frame(tukey_result, c(rep("PEL", 2), NA, NA, NA))
tukey_result[Group1=="ITU" & Group2=="HAC", `:=`(Group1="HAC", Group2="ITU")]
tukey_result[Group1=="PEL" & Group2=="LWK", `:=`(Group1="LWK", Group2="PEL")]

#61814B
#a33c34
p <-ggplot(tukey_result, aes(x=Group1, y=Group2, fill=abs(as.numeric(diff))))+
  geom_tile()+
  geom_text(aes(label=sig))+
  labs(x="", y="", fill="Mean Increase in\nNovel Transcripts\nDifference (%)")+
  scale_fill_continuous(low="white", high="#61814B", na.value="white")+
  scale_x_discrete(limits=unique(data$population)[order(unique(data$mean))])+
  scale_y_discrete(limits=unique(data$population)[order(unique(data$mean))])+
  geom_vline(xintercept=4.5, linetype="dashed")+
  annotate(geom="text", x=4.5,y=6, label="African", hjust=-0.2, size=6*0.35)+
  annotate(geom="text", x=4.5,y=6, label="OOA", hjust=1.2,  size=6*0.35)+
  mytheme+
  theme(legend.position = c(0.4,0.7))

# Add colored squares for populations (horizontal for Group1)
p + 
  geom_tile(data = tukey_result, aes(x = Group1, y = 0.375),  # Adjust `y` to place above the plot
            fill = popcol[tukey_result$Group1], width = 1, height = 0.2, alpha=0.8) +
  
  # Add colored squares for Group2 (vertical annotations)
  geom_tile(data = tukey_result, aes(x = 0.5, y = Group2),  # Adjust `x` to place left of plot
            fill = popcol[tukey_result$Group2], width = 0.2, height = 1, alpha=0.8)
ggsave("10_figures/01_plots/main/personalizedhg38/heatmap_AnovaSignificantDifferencesbetweenPOPS_perHaplotype.pdf", dpi=700, width = 2.35, height = 2.25,  units = "in")


# data <- fread("../novelannotations/analysis_tables/250206_perc_novel_sjs_uniq_to_pers_hg38.tsv")
# ggplot(data, aes(x=population, y=perc_uniq_to_haplotypes, fill=population))+
#   geom_boxplot(outliers = F, show.legend = T)+
#   ggbeeswarm::geom_quasirandom()+
#   mytheme+
#   labs(x="", y="% of Splice Junctions Unique to personalized-GRCh38s")+
#   scale_fill_manual(values=popcol)
