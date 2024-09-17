## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/03_mapping"
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

data <- fread("data_unmasked_unfiltered/1_PY1_GM10492_softclipping_filtered.tsv", header = F, sep = "\t")
colnames(data) <- c("softclipping", "readid", "alignment","geneid.v", "gene_biotype")
data[, gene_biotype := gsub(";", "", gsub("\"", "", gene_biotype))]
data[, gene_biotype := ifelse(grepl("pseudogene", gene_biotype), "pseudogene", "protein_coding")]

# Plot
ggplot(data, aes(fill= type, x=V1))+
  geom_histogram(binwidth = 40, position = "dodge")+
  xlim(c(0,200))+
  facet_wrap(~V3)+mytheme+xlab("Softclipped bases")
ggplot(data, aes(fill= type, x=type, y=V1))+
  geom_boxplot(outliers = F)+
  facet_wrap(~V3)+
  mytheme+
  ylab("Softclipped bases")+xlab("")


# widen
# datawide <- dcast(data, readid+geneid.v+gene_biotype~alignment, value.var = "softclipping")
# datawide <- dcast(data, readid~alignment, value.var = "softclipping")
# 
# datcomp<- datawide[!is.na(datawide$`tp:A:P`) & !is.na(datawide$`tp:A:P`),][, ratio := `tp:A:S`-`tp:A:P`]
# 


# Separate the tp:A:P and tp:A:S rows
tpA_P <- data[alignment == "tp:A:P"][, alignment:=NULL]
colnames(tpA_P)[1] <- "softclipping_p" 
tpA_S <- data[alignment == "tp:A:S"][, gene_biotype:=NULL][, alignment:=NULL]
colnames(tpA_S)[1] <- "softclipping_s" 


mymerge <- tpA_S[tpA_P, on=c("readid", "geneid.v")]
mydif <- mymerge[complete.cases(mymerge),][, difference := softclipping_s-softclipping_p]

ggplot(mydif, aes(y=difference, x=gene_biotype, fill=gene_biotype))+
  geom_violin(  bounds = c(-Inf, 2000))+
  geom_boxplot(outliers = F, width=0.2)+
  mytheme+
  labs(x="", y="Softclipping difference\n(2ary - 1ary)")
