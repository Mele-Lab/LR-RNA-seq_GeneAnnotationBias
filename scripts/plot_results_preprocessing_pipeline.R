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
setup_script(relative_path, 3, 48) # sets wd, load tidyverse, data.table, sets DT threads
##
##    catch arguments from standard input 
##      catch_args(number_of_args, "obj_name1", "obj_name2")
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------0

library(ggpmisc)

# get samples names that have started the preprocessing
sample_vector <- list.files("data/qc_data_multiplerun_samples","2024*", full.names = T)
sample_vector <- c(sample_vector, list.files("data/qc_data_singlerun_samples/","2024*", full.names = T))
# 1- if a samples has finished retrieve its % of original reads that remain after the filter
names_percent_original <-c()
percent_original <- list()
for(i in sample_vector[grepl("reads_percentage_bystep", sample_vector)]){
    names_percent_original <- c(names_percent_original, i)
    percent_original <- append(percent_original, list(fread(i)))}
}

names(percent_original) <- names_percent_original

# 2- if a samples has finished retrieve its read number track
names_readnum <-c()
readnum <- list()
for(i in sample_vector){
  thispath <- paste0(getwd(),"/scripts/data/", i,"/qc/count/",i,"_readnum_track.txt")
  if(file.exists(thispath)){
    names_readnum <- c(names_readnum, i)
    readnum <- append(readnum, list(fread(thispath)))}
}

names(readnum) <- names_readnum

# parse data and add number of reads
track <- rbindlist(readnum, idcol="sample")
colnames(track)[3] <- "reads"
track_wide <- dcast(track, sample ~ V1, value.var = "reads")
track_wide[, reads_kept:=`12.final_count_afterfiltering_trimmedQ7`/`1.bam_with_duplex_status`*100]


# parse data and add number of reads
percent_original <- rbindlist(percent_original, idcol="sample")
percent_original_wide <- dcast(percent_original, sample ~ phenomena2, value.var = "results2")
percent_original_wide <- track[V1=="1.bam_with_duplex_status", .(sample, reads)][percent_original_wide, on="sample"]
setnames(percent_original_wide, grep("reads_dup-assessed", colnames(percent_original_wide)), "dup_assessed")

# Plot
a<-ggplot(percent_original_wide, aes(x="",y=dup_assessed))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  ggiraph::geom_jitter_interactive(tooltip=sample)+
  labs(y="Assessed reads (%)", x="")+
  mytheme
ggplot(percent_original_wide, aes(x="",y=duplication_rate))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Duplication rate (%)", x="")+
  mytheme

# duplication rate ~ reads
cor.test(percent_original_wide$duplication_rate, percent_original_wide$reads, use="complete.obs")
ggplot(percent_original_wide, aes(x=reads/10^6, y=duplication_rate))+
  geom_point(aes(size=dup_assessed, color=dup_assessed))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.9)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Duplication(%)")+
  scale_color_gradient(low="grey", high="red")+
  mytheme


# % of duplex reads in the end
ggplot(percent_original_wide, aes(x="", y=`final_duplex_%`))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Duplex reads (% of processed reads)", x="")+
  mytheme

ggplot(percent_original_wide, aes(x=reads/10^6, y=`final_duplex_%`))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Duplex reads(%)")+
  mytheme

# 3- if a samples has finished retrieve its percentage of reads change by step
names_bystep <-c()
bystep <- list()
for(i in sample_vector){
  thispath <- paste0(getwd(), "/scripts/data/", i,"/qc/",i,"_reads_percentage_bystep.tsv")
  if(file.exists(thispath)){
    names_bystep <- c(names_bystep, i)
    bystep <- append(bystep, list(fread(thispath)))}
}

names(bystep) <- names_bystep
# parse data and add number of reads
bystep <- rbindlist(bystep, idcol="sample")
bystep_wide <- dcast(bystep, sample ~ phenomena, value.var = "results")
bystep_wide <- track[V1=="1.bam_with_duplex_status", .(sample, reads)][bystep_wide, on="sample"]



# splitting ~ thorughput
ggplot(bystep_wide, aes(x="", y=new_splitted_reads))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads from concatamers\n(% of processed reads)", x="")+
  mytheme
ggplot(bystep_wide, aes(x=reads, y=new_splitted_reads))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Concatamers(%)")+
  mytheme

# long format for qualitites
bystep_long <- melt(bystep_wide, measure.vars = c("Q<7", "Q<10"), variable.name = "Qfilter", value.name = "reads_filtered")

ggplot(bystep_long, aes(x=Qfilter, y=-reads_filtered, fill=Qfilter))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads removed (%)", x="")+
  guides(fill="none")+
  scale_fill_manual(values=c("#98B473", "#D9D76E"))+
  mytheme


# MOVE TO LENGTHS AND QUALITIES
# if a samples has finished retrieve its data qc output
names_ql <-c()
qualen <- list()
for(i in sample_vector){
  thispath <- paste0(getwd(), "/scripts/data/", i,"/qc/",i,"_qualities_lengths.tsv")
  if(file.exists(thispath)){
    names_ql <- c(names_ql, i)
    qualen <- append(qualen, list(fread(thispath)))}
}

names(qualen) <- names_ql

# parse data and add number of reads
qualen <- rbindlist(qualen, idcol="sample")
qualen_wide <- track[V1=="1.bam_with_duplex_status", .(sample, reads)][qualen, on="sample"]


# check lengths

ggplot(qualen_wide, aes(x=status, y=median_lengths, fill=status))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter(aes(size=number/10^6),alpha=0.4)+
  labs(y="Median length (bp)", x="", size="Reads (M)")+
  guides(fill="none")+
  scale_x_discrete(limits=c("all", "Simplex", "Duplex", "Not_Splitted", "Splitted"))+
  scale_fill_manual(values=c("grey", "#52822A", "#7B8CDE", "#90BB6C", "#3249BB"))+
  mytheme

ggplot(qualen_wide, aes(x=status, y=median_quals, fill=status))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Median Quality (Phred Score)", x="")+
  guides(fill="none")+
  scale_x_discrete(limits=c("all", "Simplex", "Duplex", "Not_Splitted", "Splitted"))+
  scale_fill_manual(values=c("grey", "#52822A", "#7B8CDE", "#90BB6C", "#3249BB"))+
  mytheme

ggplot(qualen_wide, aes(x=reads, y=median_lengths))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Median Lengths")+
  facet_wrap(~status)+
  mytheme

ggplot(qualen_wide, aes(x=reads, y=median_quals))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.65)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.45)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Median Qualities (Phred Score)")+
  facet_wrap(~status)+
  mytheme


# final amount of reads
ggplot(qualen_wide[status=="all"], aes(x="", y=number/10^6))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Clean #reads (M)", x="")+
  mytheme

# Q ~ length
ggplot(qualen[status=="all"], aes(x=median_lengths, y=median_quals))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  guides(line=F)+
  labs(x="reads(M)", y="Median Lengths")+
  facet_wrap(~status)+
  mytheme

lm(median_quals~median_lengths, data[status=="all"])
cor.test(qualen[status=="all",median_quals], qualen[status=="all",median_lengths])



# MOVE TO % of reads kept overall
ggplot(track_wide, aes(x="", y=reads_kept))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads kept (%)", x="")+
  mytheme

mapped <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/03_mapping/data/number_of_mapped_reads.txt")
colnames(mapped) <- c("sample", "mapped_reads")

newdata <- mapped[track[V1=="12.final_count_afterfiltering_trimmedQ7", .(sample, reads)], on="sample"]
# Mapped reads ~ reads
newdata[,mapping_ratio := mapped_reads/reads*100]

ggplot(newdata, aes(x=reads/10^6, y=mapped_reads/10^6))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.9)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  geom_point()+
  guides(line=F)+
  labs(x="Raw reads(M)", y="Mapped reads (M)")+
  geom_abline(slope=1, linetype="dashed", col="grey")+
  mytheme
# add ancestry info from the names
newdata[, ancestry:=sub(".*_.*_.*_([A-Z]{2})\\d+_.*", "\\1", sample)]
newdata$continent <- ifelse(newdata$ancestry%in%c("PY", "NI", "KE"), "AFR", "OOA")

ggplot(newdata, aes(x=continent, y=mapping_ratio))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter(aes(col=ancestry, size=reads/10^6))+
  labs(y="Mapped reads (%)", x="", size="reads (M)")+
  mytheme
wilcox.test(newdata[continent=="AFR", mapped_reads], newdata[ continent=="OOA", mapped_reads])



# Check N50s

names_nanostats <-c()
nanostats <- list()
paths <-list.files("nanoplotstats", pattern="_corr", full.names = T)
for(i in paths){
  names_nanostats <- c(names_nanostats, gsub("_preprocessed.fastq.gz_corrQlength_n50.tsv", "", gsub("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/nanoplotstats/", "", i)))
  nanostats <- append(nanostats, list(fread(i)))
}

names(nanostats) <- names_nanostats
nanostats <- rbindlist(nanostats, idcol="sample")
nanostats$sample <- gsub(".*/", "" , nanostats$sample)
# add number of reads info
nanostats <- track[V1=="12.final_count_afterfiltering_trimmedQ7", .(sample, reads)][nanostats, on=.(sample)]

# correlation of Q-length
ggplot(nanostats, aes(x="", y=corr_Qlength))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter(aes(size=reads/10^6), alpha=0.5)+
  labs(y="Pearson corr coeficient", x="", size="reads (M)")+
  mytheme

# n50

ggplot(nanostats, aes(x="", y=n50))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter(aes(size=reads/10^6), alpha=0.5)+
  labs(y="N50", x="", size="reads (M)")+
  mytheme

# reads over 1k bp

ggplot(nanostats, aes(x="", y=reads_longer1k))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter(aes(size=reads/10^6), alpha=0.5)+
  labs(y="Reads longer than 1kb (%)", x="", size="reads (M)")+
  mytheme


ggplot(nanostats, aes(x=n50, y=reads_longer1k))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.9)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  geom_point(aes(size=reads), alpha=0.5)+
  guides(line=F)+
  labs(x="N50", y="% reads longer 1kb")+
  geom_abline(slope=1, linetype="dashed", col="grey")+
  mytheme


dat <- track[V1=="12.final_count_afterfiltering_trimmedQ7"][, sampleid:=gsub(".*_","", sample)]
dat[, totalreads:=sum(reads), by=sampleid]
ggplot(dat, aes(x="", y=totalreads/10^6))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.25)+
  geom_jitter()+
  labs(y="Reads per sample (M)", x="")+
  mytheme
dat[totalreads<10000000]
dat[reads<10000000]

dat[sampleid%in% dat[reads<10000000, sampleid]][order(sampleid),]
