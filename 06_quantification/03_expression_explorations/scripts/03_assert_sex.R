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

library(data.table)
library(tidyverse)
expr <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/06_quantification/data/masked_genomic_goodindex/gencode_v47/gene_tpm_matrix.tsv")
metadata <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/00_metadata/pantranscriptome_samples_metadata.tsv")
throughput <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/03_mapping/data_pseudogene_masked/count_mapped_reads.tsv")
colnames(throughput) <- c("seq_sample_id", "fastq", "mapped", "genomic_mapped", "sirv_mapped")
#ENSG00000229807
xist <- expr[grepl("ENSG00000229807",expr$geneid.v),]
xist_long <- melt(xist, id.vars="geneid.v", variable.name = "seq_sample", value.name = "tpm")
xist_long[, seq_sample_id := gsub("gencode_", "", seq_sample)]
xist_long[, sample:= gsub(".*_", "", seq_sample)]
xist_meta <- metadata[xist_long, on=c(cell_line_id="sample")]

xist_meta_thr <- throughput[xist_meta, on="seq_sample_id"]

ggplot(xist_meta_thr, aes(x=sex, y=tpm, color=population, size=genomic_mapped/10^6))+
  geom_jitter(position = position_jitter(seed = 1,height=0))+
  mytheme+
  labs(x="", size="Mapped reads (M)")+
  geom_text(aes(label = sample), position = position_jitter(seed = 1, height=0), vjust=-1)+
  scale_color_manual(values=unique(xist_meta_thr$color_pop)[order(unique(xist_meta_thr$population))], labels=unique(xist_meta_thr$population)[order(unique(xist_meta_thr$population))])

# i need to add the througput to have a correcpopulation# i need to add the througput to have a correct sense of what is happeing


# load table of LCL sex-biased genes from The impact of sex on gene expression across human tissues
library(openxlsx)
sexbias <- read.xlsx("/home/pclavell/Downloads/aba3066-table-s2.xlsx", sheet = "LCL")
setDT(sexbias)
sexbias[, geneid := gsub("\\..*", "", ENSEMBL_gene_id)][, sexbiased :="yes"]
expr_long <- melt(expr, id.vars="geneid.v", variable.name = "seq_sample", value.name = "tpm")
expr_long[, geneid := gsub("\\..*", "", geneid.v)]

expr_long <- metadata[expr_long[, cell_line_id := gsub(".*_", "", seq_sample)], on="cell_line_id"]
exprsexbias <- unique(sexbias[, .(geneid, sexbiased)])[unique(expr_long[, .(geneid, tpm, seq_sample)]), on="geneid"]
exprsexbias <- exprsexbias[sexbiased=="yes",]
exprsexbias[, sample := tstrsplit(seq_sample, "_")[[2]]]
# remove genes with expression 0
exprsexbias <-exprsexbias[, mean := mean(tpm), by="geneid"][mean!=0,]
exprsexbias_wide <- dcast(unique(exprsexbias[, .(geneid, seq_sample,tpm)]), geneid~seq_sample, value.var = "tpm")
library(factoextra)
res.pca <- prcomp(exprsexbias_wide[, 2:ncol(exprsexbias_wide)], scale = TRUE)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
)
