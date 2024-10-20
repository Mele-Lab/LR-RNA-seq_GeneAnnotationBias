## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
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
catch_args(1, "TYPE")
##
## 0----------------------------END OF HEADER----------------------------------0
library("variancePartition")
library("edgeR")
library("BiocParallel")

# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,"_genelvl.tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  countsprefilt <- fread("../../novelannotations/quantifications/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  countsprefilt <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene.tsv", header=F)
  colnames(annot) <- c("transcriptid.v","geneid.v")
}

# PARSE DATA
metadataraw <- metadataraw[mixed_samples==FALSE]
metadataraw <- metadataraw[merged_run_mode==TRUE]
metadataraw <-metadataraw[order(metadataraw$cell_line_id),]
metadataraw <- mypca[metadataraw, on="cell_line_id"]
metadataraw[, population:=factor(population, levels=c("CEU", "AJI", "ITU", "HAC", "PEL", "LWK", "YRI", "MPC"))]

# prepare new names vector
samplesnames <- metadataraw$sample
names(samplesnames) <- metadataraw$quantification_id
samplesnames <- c(samplesnames, "transcript_id"="transcriptid.v", "geneid.v"="geneid.v")
metadataraw <- column_to_rownames(metadataraw, var="sample")
metadataunfiltered <- metadataraw[order(rownames(metadataraw)), c(covariates)]

# Sum transcript counts per gene
countsprefilt <- annot[, .(transcriptid.v, geneid.v)][countsprefilt, on=c("transcriptid.v"= "transcript_id")]
countsprefilt <- countsprefilt[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
colnames(countsprefilt) <- gsub("_1$", "", colnames(countsprefilt))
setnames(countsprefilt, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
countsprefilt <- column_to_rownames(countsprefilt, var="geneid.v")
colnames(countsprefilt)  <- gsub(".*_", "", colnames(countsprefilt))
countsprefilt <- countsprefilt[, order(colnames(countsprefilt))]

# drop one of the CEU

degreslist <- list()
droppedsample<-c()

for(i in 1:5){
droppedceu <- paste0("CEU", i)
print(droppedceu)
metadata <- metadataunfiltered[rownames(metadataunfiltered)!=droppedceu,]
counts <- countsprefilt[!grepl(droppedceu, colnames(countsprefilt))]
# if(TYPE=="pantrx"){
#   fwrite(rownames_to_column(counts, var="geneid.v"), 
#          "../../novelannotations/kallisto_quant/matrix.abundance.genelevel.tsv", row.names = F, quote = F, sep="\t")
# }else if(TYPE=="gencode"){
#   fwrite(rownames_to_column(counts, var="geneid.v"), 
#          "../../novelannotations/v47_kallisto_quant/matrix.abundance.genelevel.tsv", row.names = F, quote = F, sep="\t")
# }
# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

form <- ~population+sex+pc1+pc2

design <- model.matrix(form, metadata)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

# PREPARE CONTRASTS #########################################################################........................................................

combinations<- combn(levels(metadata$population)[!grepl("CEU", levels(metadata$population))], m=2, simplify=F)

# Generate the contrasts
contrasts <- lapply(combinations, function(pair) {
  paste0(pair[1], "vs", pair[2], " = \"population", pair[1], " - population", pair[2], "\"")
})
mycon <-paste(contrasts, collapse=", ")

# Convert the contrasts into a named vector
contrast_list <- eval(parse(text = paste0("list(", mycon, ")")))

# fit dream model with contrasts
# Define contrast matrix to test ITU vs YRI
contrast_list <- append(contrast_list,
                        list("YRIvsCEU" = "populationYRI",
                             "MPCvsCEU" = "populationMPC",
                             "LWKvsCEU" = "populationLWK",
                             "AJIvsCEU" = "populationAJI",
                             "ITUvsCEU" = "populationITU",
                             "HACvsCEU" = "populationHAC",
                             "PELvsCEU" = "populationPEL"))
contrast.matrix <- do.call(makeContrasts, c(contrast_list, list(levels = design)))


# Apply contrasts to the model fit
fit <- dream(vobjDream, form, metadata)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- variancePartition::eBayes(fit)

# extract results from first contrast
results_list <- lapply(colnames(fit), function(contrast) {
  variancePartition::topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
})
names(results_list) <- colnames(fit)

resultsdt <- rbindlist(lapply(results_list, rownames_to_column, var="geneid.v"), idcol="contrast")
resultsdt$comparison <- tstrsplit(resultsdt$contrast, "vs")[[1]]
resultsdt$reference <- tstrsplit(resultsdt$contrast, "vs")[[2]]

resultsdt[, upDEG := sum(logFC >= 0.5 & adj.P.Val < 0.05), by = "contrast"]
resultsdt[, downDEG := sum(logFC <= -0.5 & adj.P.Val < 0.05), by = "contrast"]


colnames(resultsdt)[grep("adj.P.Val", colnames(resultsdt))] <- "fdr"
colnames(resultsdt)[grep("contrast", colnames(resultsdt))] <- "contrast_name"
colnames(resultsdt)[grep("comparison", colnames(resultsdt))] <- "contrast"
colnames(resultsdt)[grep("P.Value", colnames(resultsdt))] <- "pval"

degreslist <- append(degreslist, list(resultsdt))
droppedsample <- c(droppedsample,droppedceu)}


names(degreslist) <- droppedsample
data <- rbindlist(degreslist, idcol="droppedCEU")
data[, totaldeg := upDEG+downDEG]
data[, popincontrast := factor(ifelse((reference == "CEU" & contrast == "MPC") | (reference == "MPC" & contrast == "CEU"), 
                                      "CEU & MPC in Contrast",
                                      ifelse(reference == "CEU" | contrast == "CEU", "CEU in Contrast",
                                             ifelse(reference == "MPC" | contrast == "MPC", "MPC in Contrast",
                                                    "Other Populations"))), levels=c("CEU in Contrast", "MPC in Contrast","CEU & MPC in Contrast", "Other Populations" ))]

# plot the number of degs per samples dropped
ggplot(unique(data[, .(droppedCEU, contrast_name, totaldeg, popincontrast)]), aes(x=droppedCEU, y=totaldeg, fill=popincontrast))+
  geom_col()+
  mytheme+
  facet_wrap(~contrast_name)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c(unique(metadataraw[metadataraw$population=="CEU", c("color_pop")]),
                             unique(metadataraw[metadataraw$population=="MPC", c("color_pop")]),
                             "darkgreen",
                             "darkgrey"))+
  labs(x="Dropped CEU", y="# DEGs", fill="")
