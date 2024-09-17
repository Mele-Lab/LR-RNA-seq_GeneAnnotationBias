## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/03_mapping"
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

# load data
data <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/03_mapping/data/count_mapped_reads.tsv")
colnames(data) <- c("sample", "ratio", "genome_mappedQ10", "sirv_mappedQ10", "processed_reads")


# Plot genome mapped reads MAPQ>10 ~ reads after preprocessing pipeline
ggplot(data, aes(x=processed_reads/10^6, y=genome_mappedQ10/10^6))+
  geom_abline(intercept=0, slope=1, linetype="dashed", col="grey")+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.9)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  geom_point()+
  guides(line=F)+
  labs(x="Processed reads (M)", y="Mapped reads with MAPQ>=10 (M)")+
  theme

# Plot all mapped reads MAPQ>10 ~ reads after preprocessing pipeline
ggplot(data, aes(y=genome_mappedQ10+sirv_mappedQ10, x=processed_reads))+
  geom_abline(intercept=0, slope=1, linetype="dashed", col="grey")+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.9)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  geom_point()+
  guides(line=F)+
  theme

# Plot ratio of all mapped reads MAPQ>10 / reads after preprocessing pipeline
ggplot(data, aes(x="",y=(genome_mappedQ10+sirv_mappedQ10)/processed_reads))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads MAPQ>=10 (%)", x="")+
  theme
