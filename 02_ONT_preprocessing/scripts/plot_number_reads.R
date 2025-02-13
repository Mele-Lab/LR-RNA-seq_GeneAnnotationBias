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
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(median(x), 2)))
}
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE & mixed_samples==FALSE]

metadatalong <- melt(metadata[, .(sample, reads_fastq10, reads_fastq7, map_reads_generalmap, map_reads_assemblymap)], id.vars = "sample", variable.name = "filter", value.name = "reads")
metadatalong[, type := ifelse(filter%in%c("reads_fastq10", "reads_fastq7"), "Reads", "Mapped Reads")]
ggplot(metadatalong, aes(x=filter, y=reads/1e6, fill=type))+
  geom_violin(alpha=0.5)+
  geom_boxplot(outliers =F, width=0.05, show.legend = F)+
  mytheme+
  labs(x="", y="# Final Reads (Million)", fill="")+
  scale_x_discrete(labels=c("reads_fastq10"="Min\nQ10",
                            "reads_fastq7"="Min\nQ7",
                            "map_reads_generalmap"="Genome Mapping\nQ7",
                            "map_reads_assemblymap"="Genome Mapping\nQ10"),
                   limits=c("reads_fastq7","map_reads_generalmap", "reads_fastq10","map_reads_assemblymap"))+
  scale_fill_manual(values=rev(c("darkgrey", "#4c9141")))+
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 7.5), size=6*0.35) # Adjust y for median labels
ggsave("10_figures/01_plots/supp/01_preproc_res1/violin_reads_per_sample.pdf", dpi=700, width = 3.5, height = 2.25,  units = "in")
