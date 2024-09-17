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

# load library
library("variancePartition")
library("limma")
library("edgeR")

# load data
counts <- fread("01_isoquantify/data/gencodev47/gene_counts_matrix.tsv")
colnames(counts) <- gsub(".*_", "",colnames(counts))

metadata <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv")

counts <- column_to_rownames(counts, var="geneid.v")
metadata <- metadata[cell_line_id %in% colnames(counts),]

# choose technical variables
technicalvars <- c("population", "map_reads_generalmap",   "captrap_batch", "libprep_batch", "sex")


# all batches technicalvars <- c("sex", "Date_Trizol_1", "Date_Trizol_2", "Date_RNA_extraction_1", "Date_RNA_extraction_2",  "captrap_batch", "Date_ONT_seq", "population", "generalmap")
# technicalvars <- c("sex", "Date_Trizol_1", "Date_Trizol_2", "Date_RNA_extraction_1", "Date_RNA_extraction_2",  "captrap_batch", "Date_ONT_seq", "population", "generalmap")
# technicalvars <- c( "captrap_batch","Date_Trizol_1", "Date_Trizol_2", "Date_RNA_extraction_1", "Date_RNA_extraction_2", "Date_ONT_seq")
info <- metadata[,c("cell_line_id",technicalvars), with=FALSE]
# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
form <- as.formula(paste0("~",paste(technicalvars, collapse="+")))
#"Date_Trizol_1", "Date_Trizol_2", "Date_RNA_extraction_1", "Date_RNA_extraction_2", 

# identify genes that pass expression cutoff
isexpr <- rowSums(cpm(counts) > 1) >= 0.5 * ncol(counts)

# create data structure with only expressed genes
gExpr <- DGEList(counts = counts[isexpr, ])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
info[, names(info) := lapply(.SD, as.factor)]
info[, (c(3)) := lapply(.SD, as.numeric), .SDcols = c(3)]
info <- column_to_rownames(info, var="cell_line_id")
info <- info[colnames(counts),]
design <- model.matrix(form, info)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel(vobjGenes, form, info)

# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)

# melt and plot
vplong <- melt(rownames_to_column(as.data.frame(vp), var = "geneid.v"), value.name="variance_explained")

ggplot(vplong, aes(y=variance_explained, x=variable, fill=variable))+
  geom_violin(alpha=0.5,scale="width")+
  geom_boxplot(width=0.15, outliers=F)+
  mytheme+
  guides(fill="none")+
  xlab("")
