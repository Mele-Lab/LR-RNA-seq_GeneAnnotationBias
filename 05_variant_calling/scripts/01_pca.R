## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/05_variant_calling"
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

# pcs <- fread("data/pca_mysamplesonly.eigenvec")
# pcs[, "#IID" := gsub("2024...._HS_", "", `#IID`)]
# ggplot(pcs, aes(x=PC1, y=PC2))+
#   geom_point()+
#   geom_point(data=pcs[grepl("_", `#IID`),], aes(x=PC1, y=PC2), col="red")
# 
# 
# 
# 
metadata <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv")
# metadata[, processid := paste(lab_number_sample, lab_sampleid, cell_line_id, sep="_")]
# 
# mypcs <- metadata[pcs, on=c(processid="#IID")]
# ggplot(mypcs[grepl("_", processid),], aes(x=PC1, y=PC2))+
#   geom_point(size=8, aes(col=population))+
#   scale_color_manual(values=unique(mypcs[grepl("_", processid),color_pop])[order(unique(mypcs$population))])+
#   mytheme+
#   geom_text(data=mypcs[grepl("_", processid) & grepl("PY5|NI4|CH2|IS2", processid),], aes(x=PC1, y=PC2,label=sample))
# 

# load vcf
data <- fread("finalparsed.vcf", skip=340)

# remove info columns
data <- data[,c(1:9) :=NULL]
data <- data[rowSums(data==".|.")==0,]
#dataf <- data[, c(grepl("_", colnames(data))), with=FALSE]
# 
recod <- lapply(data, function(x) ifelse(x == "0|0", 0,
                                          ifelse(x %in% c("0|1", "1|0"), 1,
                                                 ifelse(x == "1|1", 2, NA))))

# recod <- lapply(data, \(x) sub("./.:.:.:.:.", "nd", x))
recod <- as.data.frame(recod)
recod <- recod[complete.cases(recod),]
# recod <- lapply(recod, \(x) gsub(":.*", "", x))
# recod <- as.data.frame(recod)
# recod <- lapply(recod, \(x) ifelse(grepl("1/1", x), "c",
#                                    ifelse(grepl("0/1", x), "b",
#                                           ifelse(grepl("1/0", x), "b",
#                                                  ifelse(grepl("0/0", x), "a", "nd")))))

recod <- as.data.frame(recod)
# recod <- lapply(recod, \(x) as.factor(x))
# recod <- as.data.frame(recod)
# 
# recod[, countnd := grepl("nd", recod)]
recodt <- t(recod)
# recodfiltered <- recod[rowSums(recod=="nd")<30,]
pc <- prcomp(recodt)

res.pca <-rownames_to_column(as.data.frame(pc[["x"]]), var="sample_m")
setDT(res.pca)
res.pca[, cell_line_id := gsub(".*_", "", res.pca$sample_m)]
# change the cell id of the MPC5 because it had a typo
res.pca[grepl("GM10497",res.pca$cell_line_id), "cell_line_id"] <- "GM10496"
merged <- metadata[, .(cell_line_id, population, color_pop, sample)][res.pca, on="cell_line_id"]

ggplot(merged, aes(x=PC1, y=PC2))+
  geom_point(aes(col=population), alpha=0.8, size=3)+
  scale_color_manual(values=unique(merged$color_pop)[order(unique(merged$population))])+
  mytheme+
  geom_repel_text(data=merged[grepl("IS2|CH2|NI4|PY5|KE1|PE5|CH6", sample_m)], aes(x=PC1, y=PC2, label=sample_m))

library(factoextra)
fviz_pca_ind(pc, geom="text")


library("FactoMineR")
pc <- PCA(recodfiltered)
