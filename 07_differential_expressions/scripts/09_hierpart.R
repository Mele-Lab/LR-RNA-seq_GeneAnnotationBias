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

# Load libraries ####
suppressMessages(library(edgeR))
library(hier.part)
library(ggplot2)
library(ggpubr)


#Functions
current.model.mod <- function (y, current.comb, xcan, SStot=0,family = c("gaussian","quasibinomial"), 
                               link = c("logit"), gof = c("Rsqu","RMSPE"), ...)  {
  comb.data <- data.frame(xcan[, current.comb])
  colnames(comb.data) <- colnames(xcan)[current.comb]
  data <- data.frame(y, comb.data)
  depv <- names(data)[1]
  n.comb <- dim(comb.data)[2]
  xs <- vector("character", n.comb)
  for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
                                           "+", sep = "")
  xs[n.comb] <- names(comb.data)[n.comb]
  xss <- paste(xs, collapse = " ", sep = "")
  formu <- stats::formula(paste(depv, "~", xss, sep = ""))
  if (gof == "RMSPE") gf <- sqrt(sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu") {
    if (family == "quasibinomial") 
      gf <- (SStot-sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))/SStot
    if (family == "gaussian") 
      gf <- summary(stats::lm(formu, data = data))$r.squared
  }
  gf
}
all.regs.mod <- function (y, xcan, family = c("gaussian", "quasibinomial"), link = c("logit"), gof = c("Rsqu","RMSPE"),...) { 
  if (sum(is.na(xcan)) > 0) {
    missing <- is.na(apply(xcan, 1, FUN = sum))
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
            call. = FALSE)
  }
  if (sum(is.na(y)) > 0) {
    missing <- is.na(y)
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
            call. = FALSE)
  }
  pcan <- dim(xcan)[2]
  n <- (2^pcan) - 1
  combs <- combos1(pcan)$ragged
  SStot <- sum((y-mean(y))^2)
  
  if (gof == "RMSPE")  gfs <- sqrt(sum((stats::glm(y ~ 1, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu")   gfs <- 0
  
  for (i in 1:n) {
    if (i%%500 == 0) 
      cat(i, "regressions calculated:", n - i, "to go...\n")
    current.comb <- as.vector(combs[i, ][combs[i, ] > 0]) 
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse = "")
    # if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    # if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}
hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2] #5
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
}
#####---------------------------------------------------------------------------------------------------------------------------------------

## OBTAIN RESIDUALS
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

TYPE <- "pantrx"

# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,"_genelvl.tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  counts <- fread("../../novelannotations/quantifications/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  counts <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
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
metadata <- metadataraw[order(rownames(metadataraw)), c(covariates)]
# Sum transcript counts per gene
counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
counts <- counts[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
colnames(counts) <- gsub("_1$", "", colnames(counts))
setnames(counts, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
counts <- column_to_rownames(counts, var="geneid.v")
colnames(counts)  <- gsub(".*_", "", colnames(counts))
counts <- counts[, order(colnames(counts))]

# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

form <- ~pc1+pc2

design <- model.matrix(form, metadata)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)


# Apply contrasts to the model fit
fit <- dream(vobjDream, form, metadata)
residuals <- residuals(fit)

variables_OI <- c("population","sex")
# hier.part
print("Hier Part")
hier.part <- lapply(rownames(residuals), function(gene)
  hier.part.mod(y=residuals[gene,], x=metadata[variables_OI], fam = "gaussian", gof = "Rsqu"))
names(hier.part) <- rownames(residuals)

# Parse information
rsq <- sapply(rownames(residuals), function(gene) 
  sum(hier.part[[gene]]$IJ[,1])) 
names(rsq) <- rownames(residuals)

rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part), function(gene) 
                      as.numeric(unlist(hier.part[[gene]]$I.perc))))
rownames(rel_perc) <- rownames(residuals)
colnames(rel_perc) <- variables_OI

abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part), function(gene) 
                       hier.part[[gene]]$IJ[,1]*100)
)
rownames(abs_perc) <- rownames(residuals)
colnames(abs_perc) <- variables_OI

hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(variables_OI,"_rel"),paste0(variables_OI,"_abs"))

mean(hier_data$population_abs)
ggplot(hier_data, aes(x=population_abs))+
  geom_density()

fwrite(hier_data, "data/05_hierpart_geneexpression_minuspc1pc2_pluspopsex.tsv", sep="\t", quote = F, row.names = T)



######### HIER PART TRANSCRIPT EXPRESSION

# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,"_transcriptlvl.tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  counts <- fread("../../novelannotations/quantifications/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  counts <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
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
metadata <- metadataraw[order(rownames(metadataraw)), c(covariates)]
# # Sum transcript counts per gene
# counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
# counts <- counts[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
colnames(counts) <- gsub("_1$", "", colnames(counts))
setnames(counts, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
counts <- column_to_rownames(counts, var="transcriptid.v")
# colnames(counts)  <- gsub(".*_", "", colnames(counts))
counts <- counts[, order(colnames(counts))]

# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

form <- ~pc1+pc2

design <- model.matrix(form, metadata)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)


# Apply contrasts to the model fit
fit <- dream(vobjDream, form, metadata)
residuals <- residuals(fit)

variables_OI <- c("population","sex")
# hier.part
print("Hier Part")
hier.part <- lapply(rownames(residuals), function(gene)
  hier.part.mod(y=residuals[gene,], x=metadata[variables_OI], fam = "gaussian", gof = "Rsqu"))
names(hier.part) <- rownames(residuals)

# Parse information
rsq <- sapply(rownames(residuals), function(gene) 
  sum(hier.part[[gene]]$IJ[,1])) 
names(rsq) <- rownames(residuals)

rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part), function(gene) 
                      as.numeric(unlist(hier.part[[gene]]$I.perc))))
rownames(rel_perc) <- rownames(residuals)
colnames(rel_perc) <- variables_OI

abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part), function(gene) 
                       hier.part[[gene]]$IJ[,1]*100)
)
rownames(abs_perc) <- rownames(residuals)
colnames(abs_perc) <- variables_OI

hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(variables_OI,"_rel"),paste0(variables_OI,"_abs"))

mean(hier_data$population_abs)
ggplot(hier_data, aes(x=population_abs))+
  geom_density()

fwrite(hier_data, "data/05_hierpart_transcriptexpression_minuspc1pc2_pluspopsex.tsv", sep="\t", quote = F, row.names = T)


######### SPLICING RATIOS 

