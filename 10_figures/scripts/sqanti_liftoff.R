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

data <- fread("10_figures/data/sqanti_blast/collapsed_sqanti_classification.tsv")
data <- data[structural_category!="structural_category"]
data[, sample:=gsub(".*\\.", "", gsub("#.*", "", isoform))]
data[structural_category=="genic_intron", structural_category:="genic"]
plotme <- unique(data[, count:=uniqueN(isoform), by=.(sample, structural_category)][, .(sample, structural_category, count)])

colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
cats <- c("full-splice_match","incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "intergenic", "genic", "fusion", "antisense")
newcats <- c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense")
names(newcats) <- cats
names(colsqanti) <- newcats
plotme[, structural_category:=newcats[structural_category]]
plotme[, structural_category:=factor(structural_category, levels=newcats)]


median_fun <- function(x,y) {
  return(data.frame(y = y, label =  paste0(round(median(x), 0))))
}

ggplot(plotme, aes(x=structural_category, y=as.integer(count), fill=structural_category))+
  geom_boxplot(size=0.25)+
  ggbeeswarm::geom_quasirandom(show.legend = F, size=0.5)+
  scale_fill_manual(values=colsqanti)+mytheme+
  guides(fill="none")+
  labs(y="# Transcripts", x="Structural Category of Personal Assembly-built Annotation\nLifted off to GRCh38 compared to enhanced GENCODE")+
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 5000), size=6*0.35)
ggsave("10_figures/01_plots/supp/38_map_pg/boxplot.sqanti_of_personal_assembliesGTF_liftedofftoGRCh38_comparedToEnhancedGENCODE.pdf", dpi=700, width = 3.5, height = 2.25,  units = "in")


# load expression data
expr <- fread("10_figures/data/expression/collapsed_expression_in_personalassemblies.tsv")
expr <- expr[TPM!="TPM" & `#feature_id`!="__unassigned"]
colnames(expr) <- c("isoform", "TPM")
mix <- expr[data[, .(isoform, structural_category)], on=.(isoform)]

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(10^median(x), 2)))
}

mix[, structural_category:=newcats[structural_category]]
mix[, structural_category:=factor(structural_category, levels=newcats)]
ggplot(mix, aes(x=structural_category, y=as.numeric(TPM), fill=structural_category))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1, outliers = F, show.legend = F)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1), size=5*0.35) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = -1.3), size=6*0.35)+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=colsqanti)+
  guides(fill="none")+
  labs(y="# TPM", x="Structural Category of Personal Assembly-built Annotation\nLifted off to GRCh38 compared to enhanced GENCODE")+
  annotation_logticks(sides="l")
ggsave("10_figures/01_plots/supp/38_map_pg/boxplot.sqanti_of_personal_assembliesGTF_liftedofftoGRCh38_comparedToEnhancedGENCODE_expression.pdf", dpi=700, width = 3.5, height = 2.25,  units = "in")

  