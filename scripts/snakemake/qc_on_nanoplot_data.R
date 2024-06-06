## 0---------HEADER----------0
##
## AUTHOR: Pau Clavell-Revelles
## EMAIL: pauclavellrevelles@gmail.com
## CREATION DATE: 2024-05-29
##
## NOTES:
## 
##   
##
## ---------------------------


## DEFINE ARGUMENTS
args <- commandArgs(trailingOnly=TRUE)

if(length(args)>0){
  cat("CAPTURING ARGUMENTS...\n\n", sep= "")
  SAMPLE <- args[1]
  NANOUTPUT <- args[2]
  READCOUNT <- args[3]
  OUTPDF <- args[4]
  OUTSTATS <- args[5]
  OUTREADS <- args[6]
  PATHWD <- args[7]
  OUTREADS_2 <- args[8]
}
## DEFINE RELATIVE PATH (to projects/bsc83)
relative_path <- PATHWD


## set working directory for local and mn5
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")

machine <- ifelse(Sys.info()[7]=="pclavell", "local",
                  ifelse(grepl("bsc", Sys.info()[7]), "mn5",
                         "ERROR"))
if(machine=="ERROR"){
  stop("ERROR: User not recognized, this is not your laptop nor your MN5 login")
}
cat(paste0("YOU ARE WORKING IN ", machine, "... \n\n"))

if(machine=="local"){
  mn5projects <- "/home/pclavell/mounts/mn5/"
}else if(machine=="mn5"){
  mn5projects <- "/gpfs/projects/bsc83/"
}

wd <- paste0(mn5projects, relative_path, "/")
if(!dir.exists(wd)){dir.create(wd, recursive = T)}
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, "... \n\n",sep = "")
## ---------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 6, digits = 4) # non-scientific notation

## ---------------------------

## load up the packages we will need:
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
library(tidyverse)
library(data.table)
if(machine=="mn5"){setDTthreads(threads=48)}
if(machine=="local"){setDTthreads(threads=4)}
## 0---------END OF HEADER-------------------------------------------------------------------0 ##
library(patchwork)

library(ggside, lib.loc = "/gpfs/projects/bsc83/utils/Rpackages")
library(hexbin, lib.loc = "/gpfs/projects/bsc83/utils/Rpackages")

# load theme
theme <- theme_minimal() + theme(axis.ticks = element_line(linewidth = 0.2, color = "black"), axis.text = element_text(size = 11, color="black"),
                                 axis.title = element_text(size=12, vjust = -0.5, color = "black"),
                                 legend.text = element_text(size=12), legend.title = element_text(size = 12, face = "bold"),
                                 panel.border = element_rect(linewidth = 0.4, fill = FALSE), panel.background = element_blank(),   panel.grid = element_line(linewidth =0.2),
                                 strip.text = element_text(size=11),   strip.background = element_blank(),
                                 legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
                                 legend.key.size = unit(15, "pt"))


# Load data
data <- fread(NANOUTPUT)
# setwd("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/")
# data <- fread("data/4kbam/qc/4kbam_nanostats.tsv.gz")
# newreads <- fread("data/4kbam/qc/4kbam_readnum_track.txt", header=F)


# Add metadata from readID
data$duplex <- ifelse(grepl(":1", data$id), "Duplex", "Simplex")
data$splitted <- ifelse(grepl("_", data$id), "Splitted", "Not_Splitted")
data$duplex <- ordered(data$duplex, c("Simplex", "Duplex"))
data$splitted <- ordered(data$splitted, c("Not_Splitted", "Splitted"))

# Compute medians by group
data[,median_quality_duplex := median(quals), by=duplex]
data[,median_lengths_duplex := median(lengths), by=duplex]
data[,median_quality_splitted := median(quals), by=splitted]
data[,median_lengths_splitted := median(lengths), by=splitted]


duplex_median_quals <- data[, .(median_quals = median(quals)), by=duplex]
duplex_median_length <- data[, .(median_lengths = median(lengths)), by=duplex]
splitted_median_quals <- data[, .(median_quals = median(quals)), by=splitted]
splitted_median_length <- data[, .(median_lengths = median(lengths)), by=splitted]

# Compute bottom 97.5% lengths to exclude outliers in the plots
ninetysevenpointfivebottomlength <- quantile(data$lengths, 0.975)


# Plot by duplex status
a <- ggplot(data, aes(x=lengths, y=quals) ) +
  geom_hex(bins = 150) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~duplex)+
  xlim(0, ninetysevenpointfivebottomlength)+
  ggtitle(SAMPLE, "Quality~Length by Duplex Status")+
  geom_hline(data = duplex_median_quals, aes(yintercept = median_quals), linewidth=1.5, linetype=2, col="brown")+
  geom_vline(data = duplex_median_length, aes(xintercept = median_lengths), linewidth=1.5, linetype=2, col="brown")+
  geom_xsidedensity(fill="grey")+
  geom_ysidedensity(fill="grey")+
  xlab("lengths (top 2.5% skipped)")+
  ylab("read quality")+
  guides(fill = "none")+
  geom_text(data = duplex_median_quals, aes(x = ninetysevenpointfivebottomlength, y = median_quals, label = round(median_quals, 2)), vjust = -0.5, hjust=1,color = "brown") +
  geom_text(data = duplex_median_length, aes(x = median_lengths, y = 0, label = round(median_lengths, 2)), hjust = -0.5, vjust=0, color = "brown") +
  theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))

# Plot by splitting status
b <- ggplot(data, aes(x=lengths, y=quals) ) +
  geom_hex(bins = 150) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~splitted)+
  xlim(0, ninetysevenpointfivebottomlength)+
  ggtitle(SAMPLE, "Quality~Length by Split Status")+
  geom_hline(data = splitted_median_quals, aes(yintercept = median_quals), linewidth=1.5, linetype=2, col="brown")+
  geom_vline(data = splitted_median_length, aes(xintercept = median_lengths), linewidth=1.5, linetype=2, col="brown")+
  geom_xsidedensity(fill="grey")+
  geom_ysidedensity(fill="grey")+
  ylab("read quality")+
  xlab("lengths (top 2.5% skipped)")+
  guides(fill = "none")+
  geom_text(data = splitted_median_quals, aes(x = ninetysevenpointfivebottomlength, y = median_quals, label = round(median_quals, 2)), vjust = -0.5, hjust=1,color = "brown") +
  geom_text(data = splitted_median_length, aes(x = median_lengths, y = 0, label = round(median_lengths, 2)), hjust = -0.5, vjust=0, color = "brown") +
  theme+
  theme(axis.text.x = element_text(angle =45, hjust = 1))



# Reads that would be lost with Q10 filter
reads_Q10kept <- nrow(data[quals>=10,])
reads_Q10filtered <- nrow(data[quals<10,])

# Count Reads in each step
newreads <- fread(READCOUNT, header=F)
# Format data
colnames(newreads) <- c("step", "number")
newreads$number <- as.integer(newreads$number)

newreads[, original := number[1]]
newreads[, percentage := (number/original*100)]
newreads[, percentage_diff := percentage-100]

# Assign names to data
original <- newreads$number[1] # reads that have been basecalled
after_split <- newreads$number[5] # reads that we have after splitting
parent_simplex <- newreads$number[11] # amount of parent simplex
final <- nrow(data) # final amount
after_dedup <- newreads$number[13] # reads after deduplication
reads_notUMI <- newreads$number[10]
reads_UMI_unmapped <- newreads$number[12]
reads_UMI <- newreads$number[9]
preQ7filter <- newreads$number[14]
# calculations
reads_UMI_mapped <- reads_UMI-reads_UMI_unmapped
removed_duplicate_reads <- reads_UMI_mapped-after_dedup
duplicationper <- removed_duplicate_reads/reads_UMI_mapped*100
split_won <- after_split-original

## Make checks to see if numbers fit
# The reads after the splitting should be equal to the reads with UMI either mapped or unmapped and the reads without UMI
after_split==reads_UMI_mapped+reads_UMI_unmapped+reads_notUMI



# Compute interesing stats
results <- c(-parent_simplex/(after_split-removed_duplicate_reads), 
             split_won/original, 
             -removed_duplicate_reads/after_split,
             (final-preQ7filter)/preQ7filter,
             -reads_Q10filtered/final)*100
phenomena <-c("duplex", 
              "new_splitted_reads", 
              "duplicates_removed" ,
              "Q<7",
              "Q<10")
finalstats <- data.frame(results, phenomena)


results2 <- c(removed_duplicate_reads/reads_UMI_mapped,
              reads_UMI_mapped/after_split,
              nrow(data)/original,
              reads_Q10kept/original,
              sum(grepl(":1", data$id))/nrow(data))*100
phenomena2 <- c("duplication_rate", 
                "reads_dup-assessed", 
                ">Q7",
                ">Q10",
                "final_duplex_%")
finalstats2 <- data.frame(results2, phenomena2)



c <-ggplot(finalstats,aes(x=phenomena, y=results, fill = results > 0))+
  geom_col()+
  scale_fill_manual(values = c("#8B232D", "#618E3B"))+
  guides(fill = "none")+
  ylab("% of reads change")+
  xlab("")+
  scale_x_discrete(limits = c("Q<10",
                              "Q<7",
                              "duplicates_removed", 
                              "new_splitted_reads", 
                              "duplex"
                              ))+
  scale_y_continuous(expand = c(0.3, 0))+
  geom_text(aes(label=round(results, digits=1), hjust=0.5-sign(results)/2))+
  coord_flip()+
  theme

d <-ggplot(finalstats2,aes(x=phenomena2, y=results2))+
  geom_col()+
  guides(fill = "none")+
  ylab("%")+
  xlab("")+
  scale_y_continuous(expand = expansion(mult = c(0, .3)))+
  scale_x_discrete(limits = c("duplication_rate",
                              "reads_dup-assessed",
                              "final_duplex_%", 
                              ">Q7", 
                              ">Q10"
  ))+
  geom_text(aes(label=round(results2, digits=1), vjust=-0.5))+
  theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# save plots
multiplot1 <-a+b+plot_layout(nrow=1, heights = c(0.5, 0.5), guides="collect", axis_titles = "collect")
multiplot2 <-c+d+plot_layout(nrow=1, widths = c(0.6, 0.4), guides="collect")

multiplot <-multiplot1 / multiplot2+plot_layout(ncol=1, heights = c(0.8, 0.2),tag_level="new")


# Save medians in a tsv
cbind(duplex_median_length, duplex_median_quals)

duplexmerge <- duplex_median_length[duplex_median_quals, on="duplex"]
colnames(duplexmerge)[1] <- "status"
splittedmerge <- splitted_median_length[splitted_median_quals, on="splitted"]
colnames(splittedmerge)[1] <- "status"

allmerge <- rbind.data.frame(duplexmerge,splittedmerge)
# compute global stats
globalsstats<-data.frame(matrix(c("all", median(data$lengths), median(data$quals)), ncol=3))
colnames(globalsstats) <- colnames(allmerge)
allmerge <- rbind.data.frame(allmerge, globalsstats)

allmerge$number <- c(sum(data[,duplex=="Duplex"]),
                     sum(data[,duplex=="Simplex"]),
                     sum(data[,splitted=="Not_Splitted"]),
                     sum(data[,splitted=="Splitted"]),
                     nrow(data))

#### OUTPUTS
ggsave(OUTPDF, multiplot, scale=4.5)
fwrite(allmerge, OUTSTATS, quote = F, row.names = F, col.names = T, sep = "\t")
fwrite(finalstats, OUTREADS, quote = F, row.names = F, col.names = T, sep = "\t")
fwrite(finalstats2, OUTREADS_2, quote = F, row.names = F, col.names = T, sep = "\t")

