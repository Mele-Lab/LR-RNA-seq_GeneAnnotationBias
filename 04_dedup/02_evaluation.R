library(tidyverse)
library(data.table)

setwd("/home/pclavell/mounts/projects/Projects/gencode_diversity/deduplication/")

# read data
data <- fread("04_dedup/01_dedup_results/FAX00275_juncbed_transcriptome.sorted_statistics.tsv", header=F)
colnames(data) <- c("Count", "Repeats")

# add edit distance info
data <- data[Repeats!="final_umi_count"][, edit_distance:=c(rep(0,6),rep(1,6), rep(2,6), rep(3,7), rep(4,8))]
data$Repeats <- factor(data$Repeats, labels=c("Unique", 2:8))

data$edit_distance <- factor(data$edit_distance)

# compute the percentage
data <-data[, .(total=sum(Count)), by=edit_distance][data, on="edit_distance"]
data[, Percentage:=Count/total*100]

ggplot(data)+
  geom_col(aes(x=Repeats, y=Percentage, fill=edit_distance), position="dodge")+
  theme_bw()+
  ylab("% reads")+
  xlab("Times repeated")+
  labs(fill="UMI Edit distance")+
  theme(legend.position = "top",
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,4,1,1), "lines"))+
  scale_fill_manual(values=c("grey", "#76a651", "#cae34b","#eda832","#ed5732" ))+
  scale_y_continuous(expand=c(0,0))
