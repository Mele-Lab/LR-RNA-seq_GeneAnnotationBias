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
# metadata <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv")
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
data <- fread("data/fortyfiveincompletesamplestestmerge.vcf", skip=340)

# remove info columns
data <- data[,c(1:9) :=NULL]

# 
# recod <- lapply(data, function(x) ifelse(x == "0|0", 0, 
#                                           ifelse(x %in% c("0|1", "1|0"), 1, 
#                                                  ifelse(x == "1|1", 2, NA))))

recod <- lapply(data, \(x) sub("./.:.:.:.:.", "nd", x))
recod <- as.data.frame(recod)
recod <- lapply(recod, \(x) gsub(":.*", "", x))
recod <- as.data.frame(recod)
recod <- lapply(recod, \(x) ifelse(grepl("1/1", x), "c",
                                   ifelse(grepl("0/1", x), "b",
                                          ifelse(grepl("1/0", x), "b",
                                                 ifelse(grepl("0/0", x), "a", "nd")))))

recod <- as.data.frame(recod)
recod <- lapply(recod, \(x) as.factor(x))
recod <- as.data.frame(recod)

pc <- prcomp(as.data.frame(recod))

library(factoextra)
fviz_pca_var(pc, geom="text")


library("FactoMineR")
pc <- PCA(t(recod))
