## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling"
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
metadata <- fread(paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/00_metadata/data/pantranscriptome_samples_metadata.tsv"))
myfilelist <- list()
for(file in list.files("data", pattern="*count*")){
  thisfile <- fread(paste0("data/",file), header=F)
  myfilelist <- append(myfilelist, list(thisfile))
}

data <- rbindlist(myfilelist)
colnames(data) <- c("count", "category", "sampleiter")

# prepare names and colors
mycols <- c( "#0D53A3","#448AFF", "#009688", "#8BC34A", "#FFC107", "#FF9800", "#F44336", "#AD1457", "purple")
sqcategory <-c("FSM", "ISM","NIC","NNC", 
               "Genic","Intronic","Antisense","Fusion","Intergenic")
category_long <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", 
                   "genic", "genic_intron","antisense", "fusion", "intergenic")
names(sqcategory) <- category_long
names(mycols) <- sqcategory


# parse 
data[, sample := tstrsplit(sampleiter, "_")[[2]]][, iter := gsub("iter","",tstrsplit(sampleiter, "_")[[5]])][, sampleiter:=NULL]
data <- metadata[, .(lab_sampleid, population, color_pop, ooa, color_ooa, sex, color_sex)][data, on=c("lab_sampleid"="sample")]
data[, category:= sqcategory[category]]
# plot
ggplot(data, aes(x=as.integer(iter), y=count, col=category))+
  geom_jitter(alpha=0.8)+
  mytheme+
  xlab("iteration")+
  scale_color_manual(values=mycols)+
  scale_x_continuous(breaks = seq(1, 10, by = 1))
ggplot(data, aes(x=population, y=count, color=category))+
  geom_jitter(alpha=0.5)+
  mytheme+
  scale_color_manual(values=mycols)+
  xlab("")

# compute means and plot
meandata <-data[, .(mean_count = mean(count),sd = sd(count)), by=.(lab_sampleid, category)]

meandata <- metadata[, .(lab_sampleid, population, color_pop, color_ooa, color_sex, ooa, sex)][meandata, on="lab_sampleid"]

ggplot(meandata, aes(x=population, y=mean_count, color=category))+
  geom_point(position = position_jitter(seed = 1))+
  mytheme+
  scale_color_manual(values=mycols)+
  xlab("")+
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, col=category), width=.3,
                position=position_jitter(seed=1), alpha=0.5)

ggplot(meandata, aes(x=ooa, y=mean_count, fill=color_ooa))+
  geom_violin(alpha=0.5)+
  geom_boxplot(outliers=F, width=0.2)+
  geom_jitter()+
  mytheme+
  scale_fill_identity()+
  facet_wrap(~category, scales = "free")+
  ggpubr::stat_compare_means(method="wilcox.test", label="p.format", label.x.npc=0.45, label.y.npc=0.95)



# sum all sqanti categories and plot
totaldata <- metadata[, .(lab_sampleid, population, color_pop, color_ooa, color_sex, ooa, sex)][meandata[, .(count= sum(mean_count)), by="lab_sampleid"], on="lab_sampleid"]

ggplot(totaldata, aes(x=sex, y=count, fill=color_sex))+
  geom_violin(alpha=0.5)+
  geom_boxplot(outliers=F, width=0.2)+
  geom_jitter()+
  mytheme+
  scale_fill_identity()+
  ggpubr::stat_compare_means(method="wilcox.test", label="p.format", label.x.npc=0.45, label.y.npc=0.95)
ggplot(totaldata, aes(x=ooa, y=count, fill=color_ooa))+
  geom_violin(alpha=0.5)+
  geom_boxplot(outliers=F, width=0.2)+
  geom_jitter()+
  mytheme+
  scale_fill_identity()+
  ggpubr::stat_compare_means(method="wilcox.test", label="p.format", label.x.npc=0.45, label.y.npc=0.95)

# there are no differences at the transcript detection level, are there differences at the quantification level?


####### 2nd part COMPARE QUANTIFICATIONS

# load data
metadata <- fread(paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/00_metadata/pantranscriptome_samples_metadata.tsv"))
myfilelist <- list()
filenames <- c()
for(iter in list.files("data/flair/02_collapse/")){
  for(file in list.files(paste0("data/flair/02_collapse/", iter), pattern="*counts*")){
  thisfile <- fread(paste0("data/flair/02_collapse/", iter, "/", file), header=F)
  myfilelist <- append(myfilelist, list(thisfile))
  filenames <- c(filenames, file)
}}
names(myfilelist) <- filenames

data <- rbindlist(myfilelist, idcol = "sample_iter")
colnames(data)<- c("sample_iter","transcript", "counts")

# parse samples
data[, lab_sampleid:=tstrsplit(sample_iter, "_")[[2]]]
data <- metadata[, .(lab_sampleid, population, color_pop, color_ooa, color_sex, ooa, sex, sample)][data, on="lab_sampleid"]
data <- metadata[, .(lab_sampleid, sample)][data, on="lab_sampleid"]


#### Check if the discovered transcripts are always the same across iterations
myfilelist <- list()
filenames <- c()
for(file in list.files("data/sqanti/", pattern="*gtf")[71:80]){
    thisfile <- fread(paste0("data/sqanti/", file, "/", file, "_corrected.gtf"), header=F)
    myfilelist <- append(myfilelist, list(thisfile))
    filenames <- c(filenames, file)
  }
names(myfilelist) <- filenames

# parse
data <- rbindlist(myfilelist, idcol = "sample_iter")
colnames(data)<- c("sample_iter","contig", "tool", "feature", "start", "end", "unk1", "strand", "unk2", "info")

# keep only transcripts and extract ids
data <- data[feature=="transcript",][, .(sample_iter, info)]
data[, "transcriptid":= tstrsplit(info, "\"")[[2]]][, "geneid":= tstrsplit(info, "\"")[[4]]]
data <- data[, .(sample_iter, transcriptid, geneid)]
data[, sample_iter:=tstrsplit(sample_iter, "_")[[5]]]
data <- data[grepl("ENST", transcriptid)]
library(UpSetR)
binary_matrix <- data %>%
  distinct(geneid, sample_iter) %>%
  mutate(value = 1) %>%
  spread(key = sample_iter, value = value, fill = 0)

# Set the row names to geneid and remove the geneid column
rownames(binary_matrix) <- binary_matrix$geneid
binary_matrix <- binary_matrix[ , -1]

# Plot the UpSet plot
upset(as.data.frame(binary_matrix), 
      sets = colnames(binary_matrix),
      order.by = "freq")
