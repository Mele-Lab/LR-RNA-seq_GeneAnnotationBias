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

# read data
masked <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/06_quantification/data/masked_genomic_goodindex/gencode_v47/gene_tpm_matrix.tsv")
colnames(masked) <- gsub(".*_", "", colnames(masked))
unmasked <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/06_quantification/data/unmasked_genomic_wrong_index/gencode_v47/gene_tpm_matrix.tsv")
colnames(unmasked) <- gsub(".*_", "", colnames(unmasked))

removethissamples <-colnames(unmasked)[duplicated(colnames(unmasked))]
unmasked <- unmasked[,!(colnames(unmasked)%in%removethissamples),  with = FALSE]
masked <- masked[,!(colnames(masked)%in%removethissamples) & colnames(masked)%in%colnames(unmasked),  with = FALSE]


# melt

unmasked_l <- melt(unmasked, variable.name = "sample", value.name = "tpm_unmasked")
masked_l <- melt(masked, variable.name = "sample", value.name = "tpm_masked")

# merge

merged <- unmasked_l[masked_l, on=c("sample", "geneid.v")]
# load annotation info
annot <- fread("/home/pclavell/mounts/mn5/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
annot <- unique(annot[, .(geneid.v, gene_biotype)])

merged <- annot[merged, on="geneid.v"][geneid.v!="__unassigned" & c(tpm_masked>0 | tpm_unmasked>0),]

ggplot(merged[sample=="GM19117" & c(gene_biotype%in%c("protein_coding", "lncRNA") | grepl("pseudogene", gene_biotype)),], aes(x=tpm_masked, y=tpm_unmasked, colour = gene_biotype))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()


means <- merged[, mean_masked := mean(tpm_masked), by="geneid.v"]
means[, mean_unmasked := mean(tpm_unmasked), by="geneid.v"]  
means <- unique(means[, .(geneid.v, gene_biotype, mean_masked, mean_unmasked)])

biotypesvec <- unique(means[c(gene_biotype%in%c("protein_coding", "lncRNA") | grepl("pseudogene", gene_biotype)), gene_biotype])
colorvec <-c("#5f8e0a", "red","red", "#0a7c8e",  "red","red","red","red","red","red","red","red","red","red")
names(colorvec) <- biotypesvec


ggplot(means[c(gene_biotype%in%c("protein_coding", "lncRNA") | grepl("pseudogene", gene_biotype)),], aes(x=mean_masked, y=mean_unmasked, colour = gene_biotype))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  mytheme+
  scale_color_manual(values=colorvec)+
  labs(x="Log10 Mean gene expression (MASKED)", y="Log10 Mean gene expression (UNMASKED)")

