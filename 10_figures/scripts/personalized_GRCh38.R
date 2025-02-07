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



data <- fread("../novelannotations/analysis_tables/250206_perc_novel_ics_uniq_to_pers_hg38.tsv")
ggplot(data, aes(x=population, y=perc_uniq_to_haplotypes, fill=population))+
  geom_boxplot(outliers = F, show.legend = T, size=0.25)+
  ggbeeswarm::geom_quasirandom(size=1)+
  mytheme+
  labs(x="", y="% of Intron Chains Unique to personalized-GRCh38s")+
  scale_fill_manual(values=popcol)+
  guides(fill="none")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/boxplot_IntronChainsUniquetopersonalizedGRCh38s.pdf", dpi=700, width = 3, height = 2.75,  units = "in")

## test if there are differences-------------------------

# Test for Normality
shapiro.test(unique(data[, .(perc_uniq_to_haplotypes, population)])$perc_uniq_to_haplotypes) # Apply separately for each group if needed

# Test for Homogeneity of Variances
car::leveneTest(perc_uniq_to_haplotypes ~ population, data = unique(data[, .(perc_uniq_to_haplotypes, population)]))

# Perform ANOVA
anova_result <- aov(perc_uniq_to_haplotypes ~ population, data = unique(data[, .(perc_uniq_to_haplotypes, population)]))
summary(anova_result)

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
ggplot(tukey_result, aes(x=Group1, y=Group2, fill=padj))+
  geom_tile()+
  geom_text(aes(label=sig))+
  mytheme+
  labs(x="", y="", fill="FDR")+
  scale_fill_continuous(low="darkred", high="white")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/heatmap_AnovaSignificantDifferencesbetweenPOPS_perHaplotype.pdf", dpi=700, width = 3, height = 2.75,  units = "in")


data <- fread("../novelannotations/analysis_tables/250206_perc_novel_sjs_uniq_to_pers_hg38.tsv")
ggplot(data, aes(x=population, y=perc_uniq_to_haplotypes, fill=population))+
  geom_boxplot(outliers = F, show.legend = T)+
  ggbeeswarm::geom_quasirandom()+
  mytheme+
  labs(x="", y="% of Splice Junctions Unique to personalized-GRCh38s")+
  scale_fill_manual(values=popcol)
