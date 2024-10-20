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
catch_args(2, "TYPE", "LEVEL")
##
## 0----------------------------END OF HEADER----------------------------------0

TYPE<-"pantrx"
LEVEL <-"transcriptlvl"

library(edgeR)
library(FactoMineR)
library(factoextra)
# load data
if(TYPE=="gencode"){
  counts <- fread("../../novelannotations/quantifications/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  counts <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene.tsv", header=F)
  colnames(annot) <- c("transcriptid.v","geneid.v")
}


metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE]

if(LEVEL=="genelvl"){
  # Sum transcript counts per gene
  counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
  counts <- counts[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
  
  counts <- column_to_rownames(counts, var="geneid.v")
}else if(LEVEL=="transcriptlvl"){
  counts <- column_to_rownames(counts, var=colnames(counts)[grep("transcript", colnames(counts))])
}
colnames(counts) <- gsub("_.*", "",colnames(counts))
metadata <- metadata[cell_line_id %in% colnames(counts),]

# choose technical variables
factorvars <- c("sex", "trizol_batchA", "trizol_batchB", "rna_extraction_batchA", "rna_extraction_batchB",  "captrap_batch", "libprep_batch",  "population")
covariates <- c(factorvars, "cell_line_id")
# keep only expressed genes
isexpr <- rowSums(cpm(counts) > 1) >= round(0.2 * ncol(counts))
counts <- counts[isexpr,] 

# normalise counts to tmm
# make the DGEList
y <- DGEList(counts)
# calculate TMM normalization factors
y <- calcNormFactors(y)
# get the normalized counts
cpms <- cpm(y, log=FALSE)
cpms<- rownames_to_column(as.data.frame(cpms), var="geneid")
setDT(cpms)

# transposing the data.table
melted_dt <- melt(cpms, id.vars = "geneid")
transposed_dt <- dcast(melted_dt, variable ~ geneid, value.var = "value")
# add variables
metadata[, (factorvars) := lapply(.SD, as.factor), .SDcols = factorvars]
topca <- metadata[, covariates, with=F][transposed_dt[, cell_line_id := gsub(".*_", "", variable)], on="cell_line_id"]
topca <- column_to_rownames(topca, var="cell_line_id")
topca$variable <- NULL
# PCA
mypca <- FactoMineR::PCA(topca,ncp = 10, quali.sup=c(1:(length(covariates)-1)))

# # variance explained by each PC
# fviz_screeplot(mypca, ncp=10, addlabels=T)
# 
# Captrap batch (PC1 and PC2 depend on it)
# fviz_pca_ind(mypca, habillage = 6, geom = "point",mean.point = FALSE, pointsize=5,
#              repel = TRUE, axes=c(1,2))

# # ONT_seq
# fviz_pca_ind(mypca, habillage = 7, geom = "point",mean.point = FALSE, pointsize=5,
#              axes = c(1,2), addEllipses = F, repel=T)
# 
# # population
# fviz_pca_ind(mypca, habillage = 8, geom = "point",mean.point = FALSE, pointsize=5,
#              axes = c(1, 2), addEllipses = F)

# # Compute the correlation matrix
# library("corrplot")
# corrplot(mypca[["quali.sup"]][["eta2"]], is.corr=FALSE)


mypcares <- as.data.frame(mypca$ind$coord)
colnames(mypcares) <- paste0("pc", 1:10)
mypcares <- rownames_to_column(mypcares, var="cell_line_id")
fwrite(mypcares, paste0("../07_differential_expressions/data/01_PCA_", TYPE,"_",LEVEL,".tsv"), quote = F, sep="\t", row.names = F)

