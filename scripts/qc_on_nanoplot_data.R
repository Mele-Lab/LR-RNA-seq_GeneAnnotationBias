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

## DEFINE RELATIVE PATH (to projects/bsc83)
relative_path <- "Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts"

## DEFINE ARGUMENTS
args <- commandArgs(trailingOnly=TRUE)

if(length(args)>0){
  cat("CAPTURING ARGUMENTS...\n\n", sep= "")
  SAMPLE <- args[1]
  NANOUTPUT <- args[2]
  READCOUNT <- args[3]
  OUTPDF <- args[4]
  OUTTSV <- args[5]
}


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

library(ggside)
data <- fread("data/fourkbam/qc/fourkbam_nanostats.tsv.gz")
reads <- fread("data/fourkbam/qc/fourkbam_readnum_track.txt", header=F)
SAMPLE <- "mocksample"
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

ninetysevenpointfivebottomlength <- quantile(data$lengths, 0.975)


# Plot by duplex status
a <- ggplot(data, aes(x=lengths, y=quals) ) +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~duplex)+
  xlim(0, ninetysevenpointfivebottomlength)+
  ggtitle(SAMPLE, "Quality~Length by Duplex Status")+
  geom_hline(data = duplex_median_quals, aes(yintercept = median_quals), size=1.5, linetype=2, col="brown")+
  geom_vline(data = duplex_median_length, aes(xintercept = median_lengths), size=1.5, linetype=2, col="brown")+
  geom_xsidedensity(fill="grey")+
  geom_ysidedensity(fill="grey")+
  geom_text(data = duplex_median_quals, aes(x = 0, y = median_quals, label = round(median_quals, 2)), vjust = -0.5, color = "brown") +
  geom_text(data = duplex_median_length, aes(x = median_lengths, y = 0, label = round(median_lengths, 2)), hjust = -0.5, color = "brown") +
  theme

# Plot by splitting status
b <- ggplot(data, aes(x=lengths, y=quals) ) +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~splitted)+
  xlim(0, ninetysevenpointfivebottomlength)+
  ggtitle(SAMPLE, "Quality~Length by Split Status")+
  geom_hline(data = splitted_median_quals, aes(yintercept = median_quals), size=1.5, linetype=2, col="brown")+
  geom_vline(data = splitted_median_length, aes(xintercept = median_lengths), size=1.5, linetype=2, col="brown")+
  geom_xsidedensity(fill="grey")+
  geom_ysidedensity(fill="grey")+
  geom_text(data = splitted_median_quals, aes(x = 0, y = median_quals, label = round(median_quals, 2)), vjust = -0.5, color = "brown") +
  geom_text(data = splitted_median_length, aes(x = median_lengths, y = 0, label = round(median_lengths, 2)), hjust = -0.5, color = "brown") +
  theme

multiplot <-a+b+plot_layout(ncol=1, heights = c(0.5, 0.5), guides="collect")
ggsave(OUTPDF, multiplot, scale=4.5)


# Reads that would be lost with Q10 filter
percentage_reads_notfiltered <- nrow(data[quals>=10,]) / nrow(data) *100
reads_notfiltered <- nrow(data[quals>=10,])
reads_notfiltered <- data.table(matrix(c("Q10 filter", reads_notfiltered), byrow = TRUE, nrow = 1))
# Medians
# Reads in each step

reads <- fread(READCOUNT, header=F)
reads <- fread("data/fourkbam/qc/fourkbam_readnum_track.txt", header=F)

newreads <- rbind(reads, reads_notfiltered)
colnames(newreads) <- c("step", "number")
newreads$number <- as.integer(newreads$number)

newreads[, original := number[1]]
newreads[, percentage := (number/original*100)]
newreads[, percentage_diff := percentage-100]

original <- newreads$number[1]
after_split <- newreads$number[5]
after_split <- newreads$number[5]
reads_UMI_mapped <- newreads$number[9]-newreads$number[11]


c("% duplex", "% concatamers", "% duplicates", )

ggplot(newreads)+
  geom_col(aes(x=step, y=percentage))
