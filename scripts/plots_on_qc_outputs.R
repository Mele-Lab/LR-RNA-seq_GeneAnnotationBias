library(tidyverse)
library(data.table)
library(ggpmisc)

# theme
theme <- theme_minimal() + theme(axis.ticks = element_line(linewidth = 0.2, color = "black"), axis.text = element_text(size = 11, color="black"),
                                 axis.title = element_text(size=12, vjust = -0.5, color = "black"),
                                 legend.text = element_text(size=12), legend.title = element_text(size = 12, face = "bold"),
                                 panel.border = element_rect(linewidth = 0.4, fill = FALSE), panel.background = element_blank(),   panel.grid = element_line(linewidth =0.2),
                                 strip.text = element_text(size=11),   strip.background = element_blank(),
                                 legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
                                 legend.key.size = unit(15, "pt"))


# get samples names
sample_vector <- list.files("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/data/", "2024*")

# if a samples has finished retrieve its data qc output
namesvec <-c()
samplelist <- list()
for(i in sample_vector){
  thispath <- paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/data/", i,"/qc/",i,"_reads_assessed_filtered.tsv")
  if(file.exists(thispath)){
    namesvec <- c(namesvec, i)
    samplelist <- append(samplelist, list(fread(paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/data/", i,"/qc/",i,"_reads_assessed_filtered.tsv"))))}
}

names(samplelist) <- namesvec

# load tsv with throughput
throughput <- fread("Downloads/Final Cell lines Coriell.xlsx - ONT sample ID.tsv")

# parse data and add number of reads
data <- rbindlist(samplelist, idcol="sample")
dt_wide <- dcast(data, sample ~ phenomena2, value.var = "results2")
dt_wide <- throughput[,c("Reads", "Sequenced sample"), with=F][dt_wide, on=.(`Sequenced sample`=sample)]
dt_wide$Reads <- as.numeric(gsub(",", ".", dt_wide$Reads))
setnames(dt_wide, grep("reads_dup-assessed", colnames(dt_wide)), "dup_assessed")
setnames(dt_wide, grep("Sequenced sample", colnames(dt_wide)), "sample")
setnames(dt_wide, grep("Reads", colnames(dt_wide)), "reads")

# Plot
ggplot(dt_wide, aes(x="",y=dup_assessed))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Assessed reads (%)", x="")+
  theme
ggplot(dt_wide, aes(x="",y=duplication_rate))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Duplication rate (%)", x="")+
  theme

# duplication rate ~ reads
cor.test(dt_wide$duplication_rate, dt_wide$reads, use="complete.obs")
ggplot(dt_wide, aes(x=reads, y=duplication_rate, size=dup_assessed, color=dup_assessed))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.9)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Duplication(%)")+
  scale_color_gradient(low="grey", high="red")+
  theme



ggplot(dt_wide, aes(x="", y=`final_duplex_%`))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Duplex reads (% of processed reads)", x="")+
  theme

ggplot(dt_wide, aes(x=reads, y=`final_duplex_%`))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Duplication(%)")+
  theme

# 2nd part

# if a samples has finished retrieve its data qc output
namesvec <-c()
samplelist2 <- list()
for(i in sample_vector){
  thispath <- paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/data/", i,"/qc/",i,"_reads_percentage_bystep.tsv")
  if(file.exists(thispath)){
    namesvec <- c(namesvec, i)
    samplelist2 <- append(samplelist2, list(fread(thispath)))}
}

names(samplelist2) <- namesvec
# parse data and add number of reads
data2 <- rbindlist(samplelist2, idcol="sample")
dt2_wide <- dcast(data2, sample ~ phenomena, value.var = "results")
dt2_wide <- throughput[,c("Reads", "Sequenced sample"), with=F][dt2_wide, on=.(`Sequenced sample`=sample)]
dt2_wide$Reads <- as.numeric(gsub(",", ".", dt2_wide$Reads))
setnames(dt2_wide, grep("Reads", colnames(dt2_wide)), "reads")



# splitting ~ thorughput
ggplot(dt2_wide, aes(x="", y=new_splitted_reads))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads from concatamers\n(% of processed reads)", x="")+
  theme
ggplot(dt2_wide, aes(x=reads, y=new_splitted_reads))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Concatamers(%)")+
  theme

# long format for qualitites
dt2_long <- melt(dt2_wide, measure.vars = c("Q<7", "Q<10"), variable.name = "Qfilter", value.name = "reads_filtered")

ggplot(dt2_long, aes(x=Qfilter, y=-reads_filtered, fill=Qfilter))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads removed (%)", x="")+
  guides(fill="none")+
  scale_fill_manual(values=c("#98B473", "#D9D76E"))+
  theme


# MOVE TO LENGTHS AND QUALITIES
# if a samples has finished retrieve its data qc output
namesvec <-c()
samplelist <- list()
for(i in sample_vector){
  thispath <- paste0("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/data/", i,"/qc/",i,"_qualities_lengths.tsv")
  if(file.exists(thispath)){
    namesvec <- c(namesvec, i)
    samplelist <- append(samplelist, list(fread(thispath)))}
}

names(samplelist) <- namesvec

# parse data and add number of reads
data <- rbindlist(samplelist, idcol="sample")
dt_wide3 <- throughput[,c("Reads", "Sequenced sample"), with=F][data, on=.(`Sequenced sample`=sample)]
dt_wide3$Reads <- as.numeric(gsub(",", ".", dt_wide3$Reads))
setnames(dt_wide3, grep("Reads", colnames(dt_wide3)), "reads")

# check lengths

ggplot(dt_wide3, aes(x=status, y=median_lengths, fill=status))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Median length (bp)", x="")+
  guides(fill="none")+
  scale_x_discrete(limits=c("all", "Simplex", "Duplex", "Not_Splitted", "Splitted"))+
  scale_fill_manual(values=c("grey", "#52822A", "#7B8CDE", "#90BB6C", "#3249BB"))+
  theme

ggplot(dt_wide3, aes(x=status, y=median_quals, fill=status))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Median Quality (Phred Score)", x="")+
  guides(fill="none")+
  scale_x_discrete(limits=c("all", "Simplex", "Duplex", "Not_Splitted", "Splitted"))+
  scale_fill_manual(values=c("grey", "#52822A", "#7B8CDE", "#90BB6C", "#3249BB"))+
  theme

ggplot(dt_wide3, aes(x=reads, y=median_lengths))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Median Lengths")+
  facet_wrap(~status)+
  theme

ggplot(dt_wide3, aes(x=reads, y=median_quals))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.65)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.45)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Median Qualities (Phred Score)")+
  facet_wrap(~status)+
  theme


# final amount of reads
ggplot(dt_wide3[status=="all"], aes(x="", y=number/10^6))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Clean #reads (M)", x="")+
  theme
