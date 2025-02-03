## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/02_ONT_preprocessing"
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
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE & mixed_samples==FALSE]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
q7 <- fread("data/q7/all_counts")[, q:="Q7"]
q10 <- fread("data/q10/all_counts")[, q:="Q10"]
finalreads <- rbind.data.frame(q7, q10)
finalreads[, lab_number_sample:=gsub("_.*", "", V1)][, Reads:=V2][, `:=`(V1=NULL, V2=NULL)]


data <- metadata[, .(lab_number_sample, sample, population)][finalreads, on="lab_number_sample"]
data[ , Reads:=Reads/4]
# Plot final number of reads per quality
ggplot(data, aes(x=q, y=Reads/10^6))+
  geom_violin()+
  geom_boxplot(width=0.15, outliers=F)+
  ggbeeswarm::geom_quasirandom(aes(color=population))+
  mytheme+
  scale_color_manual(values=popcol)+
  labs(x="", y="Reads (M)", color="Population")

# Plot general filterings-----------------
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  paste0(round(median(x), 1), "%")))
}

# get samples names that have started the preprocessing
sample_vector <- list.files("data/qc_data_multiplerun_samples","2024*", full.names = T)
sample_vector <- c(sample_vector, list.files("data/qc_data_singlerun_samples","2024*", full.names = T))
# 1- if a samples has finished retrieve its % of original reads that remain after the filter
names_percent_original <-c()
percent_original <- list()
for(i in sample_vector[grepl("reads_assessed_filtered", sample_vector)]){
  names_percent_original <- c(names_percent_original, i)
  percent_original <- append(percent_original, list(fread(i)))}
names(percent_original) <- gsub( "_reads_assessed_filtered.tsv","" , gsub(".*/", "", names_percent_original))



percentoriginal <- rbindlist(percent_original, idcol="sample")
ggplot(percentoriginal, aes(x=phenomena2, y=results2, fill=phenomena2))+
  geom_violin(alpha=0.5, scale = "width")+
  geom_boxplot(outliers = F, width=0.05)+
  mytheme+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_fill_manual(values=c("#4c9141", "#4c9141", "darkred", "#4c9141","darkgrey", "#4c9141"))+
  guides(fill="none")+
  labs(x="", y="% of Original Reads")+
  scale_x_discrete(labels=c("final_duplex_%"="Duplex Reads", 
                            "duplication_rate"="PCR\nDuplicates", 
                            "reads_dup-assessed"="Assessed for\nPCR duplication",
                            ">Q7"="Q>=7",
                            ">Q10"="Q>=10"), 
                   limits=c("reads_dup-assessed", "duplication_rate", "final_duplex_%", ">Q7", ">Q10"))+
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = -5), size=6*0.35) # Adjust y for median labels
ggsave("../10_figures/01_plots/supp/01_preproc_res1/violin_percentages of reads.pdf", dpi=700, width = 7, height = 3,  units = "in")

