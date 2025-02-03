## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions/"
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
# Loading libraries
suppressMessages(library(caret))
suppressMessages(library(parallel))
suppressMessages(library(hier.part))
suppressMessages(library(sandwich))
suppressMessages(library(lmtest))
suppressMessages(library(multcomp))
suppressMessages(library(optparse))
options(warn=-1)

#### FUNCTIONS-------------------------------------------------------------------------------------------------------------
# Create functions to be used:
is.event.exprs <- function(event.id, event_annot,tpm.threshold = 0.5){
  if(which(rownames(psi) == event.id) %% 1000 == 0){ system(paste("echo 'Processed: ", which(rownames(psi) == event.id)," out of ", nrow(psi), "events'"))}
  # The two most abundant isforms in numeratos and denominator[!numerator] median TPM >= 1
  # Isoforms that include the event
  isoforms.in <- unlist(strsplit(event_annot[event_id==event.id, alternative_transcripts], ","))
  isoforms.all <- unlist(strsplit(event_annot[event_id==event.id, total_transcripts], ","))
  # Isoforms that excluded the event
  isoforms.out <- isoforms.all[!isoforms.all%in%isoforms.in]
  
  # Isoforms exprs value for each tissue sample
  isoforms.tpm <- lapply(metadata$sample, function(sample)
    sapply(c(isoforms.in, isoforms.out), function(isoform)
      transcript.tpm[isoform,sample]
    ))
  names(isoforms.tpm) <- metadata$sample
  
  # At least 20% of the samples express the most abundant isoform above threshold TPM ----
  # Most abundant isoform expression (median per isoform across samples)
  iso.in.most_abundant <- names(which.max(sapply(isoforms.in, function(isoform)
    median(sapply(metadata$sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))
  iso.out.most_abundant <- names(which.max(sapply(isoforms.out, function(isoform)
    median(sapply(metadata$sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))
  
  # 20% of the samples express the most abundant above "threshold" TPM
  condition1 <- sum(sapply(metadata$sample, function(sample)
    isoforms.tpm[[sample]][iso.in.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))
  condition2 <- sum(sapply(metadata$sample, function(sample)
    isoforms.tpm[[sample]][iso.out.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))
  
  return(condition1 & condition2)
  
}

mountBeta <- function(glmBatch,batchReg,mydata) {
  # mount beta coefficients
  
  nbeta <- 1
  ncoef <- 2
  beta  <- numeric(1)
  
  for (i in batchReg) {
    if (!is.factor(mydata[,i])) {
      beta[nbeta] <- glmBatch$coefficients[ncoef]
      nbeta <- nbeta + 1
      ncoef <- ncoef + 1
    }
    if (is.factor(mydata[,i])) {
      nlev        <- nlevels(mydata[,i])
      b1          <- -sum(glmBatch$coefficients[ncoef:(ncoef+nlev-2)])/nlev
      beta[nbeta] <- b1
      beta[(nbeta+1):(nbeta+nlev-1)] <- b1+glmBatch$coefficients[ncoef:(ncoef+nlev-2)] 
      nbeta <- nbeta + nlev
      ncoef <- ncoef + (nlev-1)
    }
  }
  return(beta)
}

computeResiduals <- function(mydata, batchReg, traitReg) {
  
  regressors <- c(batchReg,traitReg)
  myformeq <- paste("y", paste(colnames(mydata)[colnames(mydata)%in%regressors], collapse=" + "), sep=" ~ ")
  glmFull  <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  ylogit    <- log(mydata$y/(1-mydata$y))
  myformeqB <- paste("~ ", paste(colnames(mydata)[colnames(mydata)%in%batchReg], collapse=" + ")," - 1")
  if (length(batchReg) > 1) index   <- sapply(mydata[,batchReg],is.factor)
  if (length(batchReg) == 1) index  <- is.factor(mydata[,batchReg])
  fact      <- batchReg[index]
  if (length(fact) > 1) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors],
                                             contrasts.arg = lapply(mydata[,fact], contrasts, contrasts=FALSE))
  if (length(fact) <2 ) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors])
  
  #-- from model with only batch effects, taking out these batch effects
  
  myformeq <- as.formula(paste("y", paste(colnames(mydata)[colnames(mydata)%in%batchReg], collapse=" + "), sep=" ~ "))
  glmBatch <- glm(myformeq , data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  # only glm with batch effects
  beta <- mountBeta(glmBatch,batchReg,mydata) 
  predBatch   <- AUZ %*% matrix(beta)
  
  ylogitRes1  <- ylogit - predBatch
  yRes1       <- 1/(1+exp(-ylogitRes1))
  myformeq    <- as.formula(paste("yRes1", paste(colnames(mydata)[colnames(mydata)%in%traitReg], collapse=" + "), sep=" ~ "))
  glmRes1     <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  # -- from full model, taking out the batch effects
  SStot   <- sum((mydata$y-mean(mydata$y))^2)
  SStotR1 <- sum((yRes1-mean(yRes1))^2)
  R2Full  <- round((SStot - sum((glmFull$fitted.values- mydata$y)^2))/SStot,digits=5)
  R2Batch <- round((SStot - sum((glmBatch$fitted.values- mydata$y)^2))/SStot,digits=5)
  R2Res1  <- round((SStotR1 - sum((glmRes1$fitted.values- yRes1)^2))/SStotR1,digits=5)

  tableCoefs <- data.frame(Full=glmFull$coefficients,Batch=NA,Res1=NA)
  index   <- match(names(glmBatch$coefficients),names(glmFull$coefficients))
  tableCoefs[index,2] <- glmBatch$coefficients
  
  index   <- match(names(glmRes1$coefficients),names(glmFull$coefficients))
  tableCoefs[index,3] <- glmRes1$coefficients
  
  difGlobalMod1 <- sqrt(sum((glmFull$fitted.values-glmRes1$fitted.values)^2))/length(glmFull$fitted.values)

  retObj <- list(c(R2Full=R2Full,R2Batch=R2Batch,R2Res1=R2Res1),
                 c(difGlobalMod1=difGlobalMod1),
                 tableCoefs,
                 cleanedData=data.frame(cleanMod1=yRes1))
  
  return(retObj)  
}

get_residuals <- function(event_id, mdata){
  # Model one event at a time
  y <- pmin(pmax(as.numeric(psi[event_id,]),0),1) 
  glm_data <- cbind.data.frame(mdata, y)
  
  obj <- computeResiduals(glm_data, batchReg ,traitReg)
  return(obj$cleanedData$cleanMod1)
}

# hier.part functions ####
# ---- methods derived from hier.part
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
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}

hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2]
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  if(sum(hp$IJ$I<0)>0){
    NA
  }else{
    list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
  }
}
##### END OF FUNCTIONS-------------------------------------------------------///////////////////////////////////////



# To save:
# -- PSI & TPM of alternatively spliced events (ASE)
# -- hier.part of ASE
# -- PSI residuals of ASE
# -- Differential splicing analysis (DSA) results


# Number of CPU (cores) to parallelize mclapply
n_cores <- 112


# # Data ####
# # Gene annotation ----
# # PCG and lincRNA genes with mathched biotype annotation in gencode v38
# # PAR genes excluded
# gene_annotation <- fread("../../novelannotations/241018_v47_poder_merge.gene_transcript_plusBiotypes.tsv")
# 
# # sex.biased.genes <- gene_annotation[gene_annotation$chr=="chrY" |
# #                                     gene_annotation$gene.name=="XIST", "ensembl.id"]
# 
# # Transcript annotation ----
# transcript_annotation <- fread("../../novelannotations/")
# colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "g.biotype", "transcript.name","t.biotype")
# transcript_annotation <- transcript_annotation[transcript_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
# 
# # Exon annotation ----
# exon_annotation <- read.delim("data/public/gencode.v26.GRCh38.exons.bed", header = F) #This file is too heavy to share via github, but I can be downloaded from gencode
# colnames(exon_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "transcript.name","exon.id","exon.numer","g.biotype", "t.biotype")
# exon_annotation <- exon_annotation[exon_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
# 
# # Event annotation ----
# events.info <- readRDS("SUPPA/gencode.v26.splicing_events_coordinates.rds")




### Load data
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==FALSE & merged_run_mode==TRUE]
newnames <- metadata$sample
names(newnames) <- metadata$cell_line_id

mypca <- fread(paste0("data/01_PCA_pantrx_genelvl.tsv"))
metadata <- mypca[metadata, on="cell_line_id"]

psi <- fread("../../fairlie/240903_pt/data/suppa/poder_ALL.psi")
colnames(psi)[1] <- "event_id"
colnames(psi)[grepl("_", colnames(psi))] <- gsub("_1", "", colnames(psi)[grepl("_", colnames(psi))])
colnames(psi)[-1] <- newnames[colnames(psi)[-1]]


transcript.tpm <- fread("../../fairlie/240903_pt/data/suppa/fmt_matrix.abundance.tpm.tsv")
colnames(transcript.tpm) <- gsub("_1", "", colnames(transcript.tpm))
colnames(transcript.tpm) <- newnames[colnames(transcript.tpm)]
colnames(transcript.tpm)[1] <- "isoform"
transcript.tpm <- as.data.frame(transcript.tpm)
transcript.tpm <- column_to_rownames(transcript.tpm, var="isoform")


event_annot <- fread("../../fairlie/240903_pt/data/suppa/poder.events_ALL_strict.ioe")


# Remove events with samples with NA
number_na <- rowSums(is.na(psi),na.rm=F) # vector with number of samples with NA per event
names(number_na) <-psi$event_id
noNA_events.psi <- names(number_na)[number_na==0] # events with 0 NAs
psi <- psi[event_id%in%noNA_events.psi, ]



# 1.3 Calculate event residuals ####
# covariates <- names(metadata)[!names(metadata) %in% c("Donor", "Sample", "Smoking", "Age", "Ancestry", "Sex", "BMI")]
covariates <- c("pc1", "pc2","sex", "population")
traitReg <- c("sex", "population")
batchReg <- c("pc1", "pc2")
  
# 1.3.2 Compute residuals ----
# Residuals    20 min
psi <- as.data.frame(psi)
psi <-column_to_rownames(psi, var="event_id")

# metadata <-column_to_rownames(metadata, var="sample")
# psi <- psi[1:50,]
fr <- mclapply(rownames(psi), function(event_id) get_residuals(event_id, as.data.frame(metadata)), mc.cores = n_cores )
names(fr) <- rownames(psi)


# 1.3.3 Create dataframe ----
psi_residuals <- do.call(rbind.data.frame,
                         fr)
colnames(psi_residuals) <- colnames(psi)
rownames(psi_residuals) <- rownames(psi)
psi_residuals <- round(psi_residuals, 2)

# 1.4 Exclude events with low complexity ----

# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi.complexity <- apply(psi, 1, function(x) length(unique(x)))
# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi_residuals.complexity <- apply(psi_residuals, 1, function(x) length(unique(x)))

threshold_complexity <- 15

psi <- psi[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
                              names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
),]
psi_residuals <- psi_residuals[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
                                         names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
                                         ),]



# 1.5 Exclude events with insufficient variability ----

event_freq <- apply(psi, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
event_freq_residuals <- apply(psi_residuals, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
psi <- psi[event_freq & event_freq_residuals,]
psi_residuals <- psi_residuals[event_freq & event_freq_residuals,]

##------------ NECESSITO EXPRESSIO
# 1.6 Exclude events not sufficiently expressed ----
events_exprs <- unlist(mclapply(rownames(psi), function(i) is.event.exprs(i, event_annot,0.5),  mc.cores = n_cores  ))


# Subset ASE sufficiently exprs ----
psi <- psi[events_exprs,]
psi_residuals <- psi_residuals[rownames(psi),]
# tpm <- tpm[rownames(psi),]



# 2. hier.part ----
xcan <- as.data.frame(metadata[, .SD, .SDcols = c("sample",traitReg)])
xcan <- column_to_rownames(xcan, var="sample")
hier.part.results <- mclapply(rownames(psi_residuals), function(event.id)
  hier.part.mod(y=as.numeric(psi_residuals[event.id,]), xcan=xcan,
                fam = "quasibinomial", link = "logit", gof = "Rsqu",control = list(maxit = 100)), mc.cores = n_cores)
names(hier.part.results) <- rownames(psi_residuals)

# 2.2 Exclude events with negative estimates ----
print(paste0("ASE events with unestimable contributions: ", sum(is.na(hier.part.results))))
if(length(which(is.na(hier.part.results)))>0){
  psi <- psi[-which(is.na(hier.part.results)),]
  psi_residuals <- psi_residuals[-which(is.na(hier.part.results)),]
  hier.part.results <- hier.part.results[-which(is.na(hier.part.results))]
}


# 2.3 Parse results ----
rsq <- sapply(rownames(psi_residuals), function(event)
  sum(hier.part.results[[event]]$IJ[,1]))
names(rsq) <- rownames(psi_residuals)
rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part.results), function(event)
                      as.numeric(unlist(hier.part.results[[event]]$I.perc))))
rownames(rel_perc) <- rownames(psi_residuals)
colnames(rel_perc) <- traitReg
abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part.results), function(event)
                       hier.part.results[[event]]$IJ[,1])
)
rownames(abs_perc) <- rownames(psi_residuals)
colnames(abs_perc) <- traitReg
hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(traitReg,"_rel"),paste0(traitReg,"_abs"))
hier_data <- rownames_to_column(hier_data, var="eventid")

fwrite(hier_data, "data/05_hierpart_splicing_minuspc1pc2_pluspopsex.tsv", sep="\t", row.names = F, quote = F)
