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
catch_args(3, "DIR_REL_TO_06", "FEATURE", "MAGNITUDE")
##
## 0----------------------------END OF HEADER----------------------------------0
# DIR_REL_TO_06 <- "01_isoquantify/data/gencode_v47"
# FEATURE <- "gene"
# MAGNITUDE <- "tpm"

# load data
samples_names <- c()
quantifications <- list()
for(sample in list.files(DIR_REL_TO_06)){
  loaded_sample <- fread(paste0(DIR_REL_TO_06, "/", sample, "/", sample, ".", FEATURE, "_", MAGNITUDE, ".tsv"))
  quantifications <- append(quantifications, list(loaded_sample))
  samples_names <- c(samples_names, sample)
}
names(quantifications) <- samples_names
quantifications <- rbindlist(quantifications, idcol="sample")
featureidv <-paste0(FEATURE,"id.v")
colnames(quantifications)[2] <- featureidv

# prepare formula and dcast
myformula <- as.formula(substitute(feature~sample, list(feature=as.name(featureidv))))
quanti <- dcast(quantifications, myformula,  fill=NA )


# save data
fwrite(quanti, paste0(DIR_REL_TO_06, "/", FEATURE, "_", MAGNITUDE, "_matrix.tsv"), sep="\t")

