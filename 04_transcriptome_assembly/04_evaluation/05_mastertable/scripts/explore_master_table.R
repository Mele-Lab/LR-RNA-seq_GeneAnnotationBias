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
library(ggpmisc)

data <- fread("04_transcriptome_assembly/04_evaluation/02_sqanti/data/240909merge/240909merge_reslongmeta_annot_sqanti_sj.tsv")
metadata <- fread("00_metadata/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
# EXPLORE sample SHARING
sample_sharing <- as.data.frame(data[,.SD, .SDcols = c(grep("isoform", colnames(data)),grep("^[A-Z]{3}[0-9]{1}$",colnames(data)))])
sample_sharing <- column_to_rownames(sample_sharing, var="isoform")
UpSetR::upset(sample_sharing,
              sets = colnames(sample_sharing),
              order.by = "freq")

# this plots make me wonder about the correlation between mapped reads and total discovered isoforms
sample_discovery <-as.data.frame(colSums(sample_sharing))
colnames(sample_discovery) <- "detected_isoforms"
sample_discovery <- rownames_to_column(sample_discovery, var="sample")
sample_discovery <- metadata[, .(sample, map_reads_assemblymap, population)][sample_discovery, on="sample"]

ggplot(sample_discovery, aes(x=map_reads_assemblymap/10^6, y=detected_isoforms/10^3, color=population))+
  stat_poly_line(color="darkgrey") +
  geom_point()+
  scale_color_manual(values=popcol)+
  xlab("Mapped Reads (M)")+
  ylab("# Detected Transcripts (x1000)")+
  mytheme+
  stat_poly_eq(use_label(c("eq", "R2")), color="black")
# what about transcript sharing across samples
transcript_share <- as.data.frame(table(rowSums(sample_sharing)))
colnames(transcript_share) <- c("samples", "transcripts")
transcript_share$samples <-as.numeric(transcript_share$samples)
ggplot(transcript_share, aes(x=samples, y=transcripts))+geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme

# Define custom breaks
breaks <- c(0, 1,3, 43)
# Add the 'bins' column
transcript_share$bins <- cut(transcript_share$samples, breaks = breaks, include.lowest = F, right = TRUE, labels=c("1", "2-3", ">3"))
setDT(transcript_share)
transcript_share[, transcripts_bin := sum(transcripts), by="bins"]
ggplot(unique(transcript_share[, .(bins, transcripts_bin)]), aes(x=bins, y=transcripts_bin))+geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  geom_text(aes(label=transcripts_bin), vjust=-0.5)+
  ylim(0, max(transcript_share$transcripts_bin) * 1.1)  # Increase the y-axis limit by 10%

# now the same by sqanti category
transcript_sharing_sqanti <- as.data.frame(rowSums(sample_sharing))
colnames(transcript_sharing_sqanti) <- c("samples")
transcript_sharing_sqanti$isoform <- rownames(transcript_sharing_sqanti)

transcript_sharing_sqanti <- data[, .(isoform, structural_category)][transcript_sharing_sqanti, on="isoform"]
transcript_sharing_sqanti[, transcripts := .N, by="structural_category"]

colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C", "#B5B5B5", "#4F4F4F", "#6E5353", "#000000")
names(colsqanti) <- unique(data$structural_category)[c(1,2,3,5,8,6,7,4)]

transcript_sharing_sqanti$structural_category <- factor(transcript_sharing_sqanti$structural_category, levels=names(colsqanti))

# plot
ggplot(unique(transcript_sharing_sqanti[, .(samples,isoform,structural_category)]), aes(x=samples, fill=structural_category))+
  geom_bar(position="stack")+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  scale_fill_manual(values=colsqanti)
ggplot(unique(transcript_sharing_sqanti[, .(samples,isoform,structural_category)]), aes(x=samples, fill=structural_category))+
  geom_bar(position="fill")+
  labs(x="# Samples where a transcript is found", y="Proportion of  Transcripts")+
  mytheme+
  scale_fill_manual(values=colsqanti)

# Define custom breaks
breaks <- c(0, 1,3, 43)
# Add the 'bins' column
transcript_sharing_sqanti$bins <- cut(transcript_sharing_sqanti$samples, breaks = breaks, include.lowest = F, right = TRUE, labels=c("1", "2-3", ">3"))
setDT(transcript_sharing_sqanti)
transcript_sharing_sqanti[, transcripts_bin := .N, by=c("bins", "structural_category")]
toplot<- unique(transcript_sharing_sqanti[, .(bins, transcripts_bin, structural_category)])
ggplot(toplot, aes(x=bins, y=transcripts_bin, fill=structural_category))+
  geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  geom_text(aes(label=transcripts_bin),position = position_stack(vjust = 0.5))+
  ylim(0, 190000)+
  scale_fill_manual(values=colsqanti)

# ABOUT SJ
# is recount support correlated with sample sharing?
ggplot(data, aes(x=sample_sharing, y=log10(sj_less_recountsupported_counts+1), color=structural_category))+
  stat_poly_line(color="darkred") +
  geom_point(alpha=0.2)+
  mytheme+
  ylab("Recount3 counts of\nleast supported sj")+
  scale_color_manual(values=colsqanti)+
  facet_wrap(~structural_category)+
  stat_poly_eq(use_label(c("eq", "R2")), color="black")
library(ggridges)

# how is the distribution of recount support by structural category?
data$structural_category <- factor(data$structural_category, levels=rev(names(colsqanti)))
ggplot(data, aes( x=log10(sj_less_recountsupported_counts+1), y=structural_category, fill=structural_category))+
  geom_density_ridges() +  
  mytheme+
  xlab("Recount3 log10(counts) of\nleast supported SJ")+
  ylab("")+
  scale_fill_manual(values=colsqanti)

ggplot(data, aes( x=log10(sj_less_recountsupported_counts+1), y=structural_category, fill=structural_category))+
  geom_density_ridges() +  
  mytheme+
  xlab("Recount3 log10(counts) of\nleast supported SJ")+
  ylab("")+
  scale_fill_manual(values=colsqanti)

breaks <- c(0, 1,3, 43)
