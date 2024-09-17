## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/02_isoquant"
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

# prepare list of files
sqantistats_filelist <-list.files("stats", pattern=".novel_vs_known.SQANTI-like.tsv_sqantistats")

# load files
sqantistats_names <- c()
sqantistats_list <- list()
for(sample in sqantistats_filelist){
  sqantistats_names <- c(sqantistats_names, gsub("\\..*", "", sample))
  sqantistats_list <- append(sqantistats_list, list(fread(paste0("stats/",sample))))
  }
names(sqantistats_list) <- sqantistats_names
# parse files
sqantistats <- rbindlist(sqantistats_list, idcol="sample")
colnames(sqantistats)<-c("sample", "novel_isoforms", "category")
# add metadata
sqantistats[, total_novel:=sum(novel_isoforms), by=sample]
sqantistats[, ancestry:=gsub(".$","",lapply(strsplit(sqantistats$sample,"_"), \(x) x[[4]]))]
sqantistats$continent <- ifelse(sqantistats$ancestry %in% c("NI", "PY", "KE"), "AFR", "OOA")

# add info about the number of reads used
# 2- if a samples has finished retrieve its read number track
readnum <- list()
for(i in sqantistats_names){
  thispath <- paste0("../../ONT_preprocessing/scripts/data/", i,"/qc/count/",i,"_readnum_track.txt")
  if(file.exists(thispath)){
    readnum <- append(readnum, list(fread(thispath)))}
}

names(readnum) <- sqantistats_names

# parse data and add number of reads
throughput <- rbindlist(readnum, idcol="sample")
colnames(throughput)[3] <- "reads"
sqantistats <- throughput[V1=="12.final_count_afterfiltering_trimmedQ7",.(sample, reads)][sqantistats, on="sample"]
# add info about detected isoforms in general, not only novel
detectediso <- fread("stats/count_discovered_transcripts.tsv")
colnames(detectediso) <- c("sample", "detected")
detectediso$sample <- gsub(".*/", "", gsub("\\..*", "", detectediso$sample))
sqantistats <-detectediso[sqantistats, on="sample"]

sqantistats[, known_isoforms :=detected-total_novel]

sqantistats_wide <-dcast(sqantistats, sample+total_novel+detected+ancestry+continent+reads+known_isoforms~category, value.var="novel_isoforms" )
sqantistats_long <- melt(sqantistats_wide, 
                         id.vars = c("sample", "reads", "continent", "ancestry", "detected", "total_novel"),
                         variable.name="category",
                         value.name="isoforms",
                         value.factor = )
# plots

ggplot(sqantistats)+
  geom_col(aes(x=sample, y=novel_isoforms, fill=category))
ggplot(sqantistats_long, aes(x=sample, y=isoforms, fill=category))+
  geom_col(pos="fill")+
  xlab("")+
  facet_wrap(~continent, scales="free_x")+
  mytheme
  
library(ggpmisc)
ggplot(sqantistats_long, aes(x=reads/10^6, y=isoforms, col=category))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2","p")))  +
  geom_point()+
  xlab("reads(M)")+
  ylab("Detected isoforms")+
  scale_color_manual(values=c("black", "green", "blue", "red" , "yellow", "orange", "darkgrey", "grey" ))+
  facet_wrap(~continent)+
  mytheme

metadata <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/00_metadata/pantranscriptome_samples_metadata.tsv")

sqantistats <- metadata[,.(two_letter_pop, color, population)][sqantistats, on=c(two_letter_pop="ancestry"),allow.cartesian=TRUE]

mycolors <- unique(metadata$color)
names(mycolors) <- unique(metadata$population)

ggplot(unique(sqantistats[category=="intergenic", .(population, novel_isoforms, reads, color)]),
       aes(x=population, y=novel_isoforms, size=reads/10^6, color=population))+
  scale_color_manual(values=mycolors)+
  geom_point(alpha=0.5)+
  ylab("intergenic loci")+
  xlab("")+
  mytheme
