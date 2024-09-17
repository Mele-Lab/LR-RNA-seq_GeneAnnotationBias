## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/06_quantification"
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

# load gene quantifications
gene <- fread("data/unmerged/gene_tpm_matrix.tsv")
colnames(gene) <- c("geneid.v","CH1_1","NI2_1","KE2_1","KE2_2","CH1_2","NI2_2")
gene <-gene[, sumrow:=CH1_1+CH1_2+NI2_1+NI2_2+KE2_1+KE2_2][sumrow>0,]

# correlations and plot
lm(gene$NI2_1~ gene$NI2_2)
cor(gene$NI2_1,gene$NI2_2)
ggplot(gene, aes(x=NI2_1, y=NI2_2))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.85)+
  mytheme
sum(gene$NI2_1>0)
sum(gene$NI2_2>0)

# load transcript quantifications
transcript <- fread("data/unmerged/transcript_tpm_matrix.tsv")
colnames(transcript) <- c("geneid.v","CH1_1","NI2_1","KE2_1","KE2_2","CH1_2","NI2_2")
transcript <-transcript[, sumrow:=CH1_1+CH1_2+NI2_1+NI2_2+KE2_1+KE2_2][sumrow>0,]

# correlations and plot

lm(transcript$NI2_1~ transcript$NI2_2)
cor(transcript$NI2_1,transcript$NI2_2)
cor.test(transcript$NI2_1,transcript$NI2_2)
sum(transcript$NI2_1>0)
sum(transcript$NI2_2>0)

ggplot(transcript, aes(x=NI2_1, y=NI2_2))+
  geom_point()+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")))  +
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.85)+
  mytheme

# merge and plot together
merging <- data.frame("feature"=c("gene", "gene","transcript", "transcript"),
                      "expressed"=c(sum(gene$NI2_1>0),sum(gene$NI2_2>0), sum(transcript$NI2_1>0),sum(transcript$NI2_2>0)),
                      "sample"=c("NI2_1", "NI2_2", "NI2_1", "NI2_2"))
ggplot(merging, aes(y=expressed, x=feature, fill=sample))+
  geom_col(position="dodge")+
  mytheme+
  labs(x="", y="expressed features", fill="")+
  scale_fill_manual(values=c("#533745","#9D9171"))

### NOW check if reruns express the same transcripts
# GENE
gene[, (2:ncol(gene)) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = 2:ncol(gene)]
newgene <-column_to_rownames(gene, var="geneid.v")
# Plot the UpSet plot
upset(as.data.frame(newgene)[,c(2,6)], 
      sets = colnames(newgene)[c(2,6)],
      order.by = "freq")

# TRANSCRIPT
oldtranscript <- copy(transcript)
transcript[, (2:ncol(transcript)) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = 2:ncol(transcript)]
newtranscript<-column_to_rownames(transcript, var="geneid.v")
# Plot the UpSet plot
upset(as.data.frame(newtranscript)[,c(2,6)], 
      sets = colnames(newtranscript)[c(2,6)],
      order.by = "freq")

transcript[, countrow:=CH1_1+CH1_2+NI2_1+NI2_2+KE2_1+KE2_2]
sharingquant <- oldtranscript[transcript[, .(geneid.v, countrow)], on="geneid.v"]
sharingquant[, meancount:=(CH1_1+CH1_2+NI2_1+NI2_2+KE2_1+KE2_2)/6]
# are the shared transcripts more expressed than the non-shared?

ggplot(sharingquant, aes(x=countrow, y=log10(meancount)))+
  geom_violin(aes(fill=as.factor(countrow)), alpha=0.5)+
  geom_boxplot(width=0.1, aes(fill=as.factor(countrow)))+
  mytheme+
  scale_x_continuous(breaks = 0:6)+
  guides(fill="none")+
  labs(x="Sample Sharing", y="Log10(Mean counts)")+
  geom_hline(yintercept=0, linetype="dashed", col="grey")
# samples are repeated (reruns) but considering the results of the following upset plot i think that it is comparable
upset(as.data.frame(newtranscript)[,1:6], 
      sets = colnames(newtranscript[1:6]),
      order.by = "freq")
