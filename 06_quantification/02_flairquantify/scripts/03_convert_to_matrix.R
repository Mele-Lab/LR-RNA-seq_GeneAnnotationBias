## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/06_quantification/02_flairquantify/"
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

data <- fread("data/240917_merge_geneEntry_correctedScaffolds_nochrEBV/concat_flair_stringent_counts.tsv")
colnames(data) <- c("isoform", "counts", "sample")

# convert to wide
datawide <- dcast(data, isoform~sample, value.var="counts", fill = 0)
datawidena <- dcast(data, isoform~sample, value.var="counts")

# total counts
stats <- data.frame("isoform"=datawide$isoform)
stats$total_counts <- rowSums(datawide[, .SD, .SDcols = 2:ncol(datawide)])
stats$mean_counts <- rowMeans(datawidena[, .SD, .SDcols = 2:ncol(datawide)], na.rm=T)
stats$expressed_samples <- rowSums(datawide[, .SD, .SDcols = 2:ncol(datawide)]>0)
stats$min_counts <- apply(datawidena[, .SD, .SDcols = 2:ncol(datawide)], 1, min, na.rm=T)
stats$max_counts <- apply(datawidena[, .SD, .SDcols = 2:ncol(datawide)], 1, max, na.rm=T)

colnames(stats)[!grepl("isoform", colnames(stats))] <- paste0("flair_", colnames(stats)[!grepl("isoform", colnames(stats))] )

fwrite(stats, "data/240917_merge_geneEntry_correctedScaffolds_nochrEBV/quantification_stats_stringent_for_master_table.tsv", quote = F, row.names = F, sep="\t")
