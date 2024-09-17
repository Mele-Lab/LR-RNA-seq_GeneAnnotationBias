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
un <- fread("data_pseudogene_masked/number_of_mapped_reads.txt")
ub <- fread("data_pseudogene_masked_ub/number_of_mapped_reads.txt")
or <- fread("data/count_mapped_reads.tsv")

colnames(or) <- c("lab_id", "ratio", "genomic", "sirv", "fastq")
colnames(un) <- c("sample", "mapped_un")
colnames(ub) <- c("sample", "mapped_ub")

together <- un[ub, on="sample"]

# modify or
or[, sample:= gsub("2024...._HS_", "", lab_id)]
or[, mapped_unmasked := genomic+sirv]
or[, mapped_unmasked_sample := sum(mapped_unmasked), by="sample"]
or[, fastq_reads := sum(fastq), by="sample"]

# merge
fmerge <- unique(or[, .(sample, mapped_unmasked_sample, fastq_reads)])[together, on="sample"]


# plot to compare -un and -ub parameters
ggplot(together, aes(x= mapped_un, y=mapped_ub))+
  geom_point()
cor(together$mapped_un, together$mapped_ub)

ggplot(fmerge, aes(x=mapped_unmasked_sample/10^6, y=mapped_un/10^6))+
  geom_point()+
  geom_abline(slope=1, color="grey", linetype="dashed")+
  coord_fixed(ratio = 1, xlim = c(8, 20), ylim = c(8, 20), expand = TRUE, clip = "on")+
  mytheme+
  labs(x="Reads mapped to genome (M)", y="Reads mapped to\npseudogene-masked genome (M)")
