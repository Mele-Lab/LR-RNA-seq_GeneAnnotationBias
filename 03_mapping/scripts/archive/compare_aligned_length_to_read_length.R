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
sample_vector <- list.files("ONT_preprocessing/scripts/data/","2024*")

# 1- if a samples has finished retrieve its % of original reads that remain after the filter
names_percent_original <-c()
percent_original <- list()
for(i in sample_vector){
  thispath <- paste0(getwd(),"ONT_preprocessing/scripts/data/", i,"/qc/",i,"_qualities_lengths.tsv")
  if(file.exists(thispath)){
    names_percent_original <- c(names_percent_original, i)
    percent_original <- append(percent_original, list(fread(paste0("ONT_preprocessing/scripts/data/", i,"/qc/",i,"_qualities_lengths.tsv"))))}
}

names(percent_original) <- names_percent_original
stats <- rbindlist(percent_original, idcol="sample")
stats <- stats[status =="all",]


# 2- load alignment info
sample_vector <- list.files("03_mapping/data_pseudogene_masked/genomic_nanoplot/","*_nanoplotNanoStats.txt")

names_align <-c()
align <- list()
for(i in sample_vector){
    names_align <- c(names_align, i)
    align <- append(align, list(fread(paste0(getwd(),"/03_mapping/data_pseudogene_masked/genomic_nanoplot/", i))))
}

names(align) <- names_align
align <- rbindlist(align, idcol="sample")
align <- align[Metrics =="median_read_length",]
align[, jointsamples := gsub("2024...._HS_","",gsub("_nanoplotNanoStats.txt", "", sample))]

# compute ponderated raw lengths
stats[, jointsamples := gsub("2024...._HS_", "", sample)]
stats[, totalnum := sum(number), by="jointsamples"][, propnum := number/totalnum][, lengths_pon := median_lengths*propnum][,rawlength := sum(lengths_pon), by="jointsamples" ]

# merge data
tog <- align[unique(stats[,.(rawlength, jointsamples)]), on="jointsamples"]

library(ggpmisc)
ggplot(tog, aes(x=rawlength, y=as.numeric(dataset)))+
  labs(x="median raw read length",
       y="median aligned(MAPQ>10) read length")+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey")+
  xlim(c(465, 620))+
  ylim(c(465, 620))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.x=0.3)+
  geom_point()+
  coord_flip()+
  mytheme
