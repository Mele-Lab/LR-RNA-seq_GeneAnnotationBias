## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
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
library(edgeR)
# load data
poder <- fread("../../novelannotations/kallisto_quant/matrix.abundance.genelevel.tsv")
gencode <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.genelevel.tsv")
metadata<- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==FALSE]
metadata <- metadata[merged_run_mode==TRUE]

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
# keep only genes identified in both
poder <- poder[geneid.v%in%gencode$geneid.v]
gencode <- gencode[geneid.v%in%poder$geneid.v]

# compute cpm
podercpm <-cpm(DGEList(poder))
podercpm<-data.table(data.frame(podercpm))
podercpm$geneid.v <- poder$geneid.v
gencodecpm <-cpm(DGEList(gencode))
gencodecpm<-data.table(data.frame(gencodecpm))
gencodecpm$geneid.v <- gencode$geneid.v
# format long and merge
poderlong <- melt(podercpm, id.vars = "geneid.v", value.name = "poder_expression", variable.name = "sample")
gencodelong <- melt(gencodecpm, id.vars = "geneid.v", value.name = "gencode_expression", variable.name = "sample")

data <- poderlong[gencodelong, on=c("geneid.v", "sample")]
data[, population := gsub(".$", "", sample)]


# plot
ggplot(data, aes(x=gencode_expression+1, y=poder_expression+1, col=population))+
  geom_point(alpha=0.5)+
  mytheme+
  scale_color_manual(values=popcols)+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(trans="log10")+
  facet_wrap(~population, nrow = 2)+
  annotation_logticks(sides="bl")+
  ggpmisc::stat_poly_line(color="black") +
  ggpmisc::stat_poly_eq(use_label(c("eq", "adj.R2")),label.y = 0.9, color="black")+
  ggpmisc::stat_poly_eq(use_label(c( "p", "n")),label.y = 0.8, color="black")+
  guides(color="none")+
  labs(x="GENCODEv47 quantification", y="PODER quantification")


data[, ratiopodergencode := poder_expression/gencode_expression]
data[, ratiodisruption:=ifelse(ratiopodergencode>1.1, ">10% PODER overexpression", 
                               ifelse(ratiopodergencode<0.9, ">10% GENCODE overexpression",
                                      "No differences"))]
data <- data[poder_expression!=0 | gencode_expression!=0,]

ggplot(data, aes(x=population, y=ratiopodergencode, fill=population))+
  geom_violin(scale="width", alpha=0.7)+
  geom_boxplot(outliers = F, width=0.1)+
  mytheme+
  scale_fill_manual(values=popcols)+
  scale_y_continuous(trans="log10")+
  annotation_logticks(sides = "l")

ggplot(data, aes(x=ratiodisruption, col=sample, fill=population))+
  geom_bar(position="dodge")+
  mytheme+
  scale_fill_manual(values=popcols)+
  scale_color_manual(values=rep("black",43))+
  guides(col="none")+
  labs(y="# Genes")+
  