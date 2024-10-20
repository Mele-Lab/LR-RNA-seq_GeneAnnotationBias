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
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------0

# libraries ----
library(manta)
library("sQTLseekeR2")
library(parallel)
library(reshape2)


#####################------------------------------------------------------------------------------------
library(edgeR)
TYPE="pantrx"
# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,".tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  counts <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  counts <- fread("../../novelannotations/kallisto_quant/matrix.abundance.tsv")
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

# Sum transcript counts per gene
counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
colnames(counts) <- gsub("_1$", "", colnames(counts))
setnames(counts, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
colnames(counts)[1:2] <- c("trId", "geneId")

# genecounts <- counts[, lapply(.SD, sum), by = geneId, .SDcols = patterns("[a-zA-Z]{3}[0-9]")]
## keep transcripts that have at least 11 samples with more than 0.1 cpm
keeptranscripts <-rowSums(cpm(counts[,.SD, .SDcols=grep("[0-9]$", colnames(counts))]) > 0.1)>10

trxcounts <- counts[keeptranscripts][, trxpergene := uniqueN(trId) , by="geneId"]
trxcounts <- trxcounts[trxpergene>1]

genecounts <- trxcounts[, lapply(.SD, sum), by = geneId, .SDcols = patterns("[0-9]$")]
## keep genes with at least 11 samples with at least 10 counts
keepgenes <-rowSums(genecounts[,.SD, .SDcols=grep("[0-9]$", colnames(genecounts))]>=10) >10
filteringgenes <- data.table(geneId=genecounts$geneId, keep=keepgenes)
trxcountsfiltered <- filteringgenes[trxcounts, on="geneId"]
trxcountsfiltered <- trxcountsfiltered[keep==TRUE][, keep:=NULL]


meltedtrx <- melt(trxcountsfiltered, id.vars = c("geneId", "trId", "trxpergene"), value.name = "trx_counts", variable.name="sample")
setDT(meltedtrx)
meltedtrx[, trx_relative_abundance := trx_counts/sum(trx_counts), by=c("sample", "geneId")]
# relab2fairlie <-dcast(meltedtrx, trId~sample, value.var="trx_relative_abundance", fill=0)
# relab2fairlie2 <-annot[relab2fairlie, on=c("transcriptid.v"="trId")]
# fwrite(relab2fairlie2, "data/04_relative_abundances_allgenes_gencode.tsv", sep="\t", quote = F, row.names = F)


# this function proably doesn't work because it is filtering genes even when using filters at 0
# prepared_counts <- prepare.trans.exp(te.df = as.data.frame(trxcountsfiltered),
#                                      min.transcript.exp = 0,
#                                      min.gene.exp = 0,
#                                      min.prop = 0,
#                                      min.dispersion = 0,
#                                      verbose=T)
total_number_distinct_genes <- length(unique(meltedtrx$geneId))



# 2. run manta ----
# Diego's email -> Por otro lado, cuando solo queréis los covariates para corregir y testáis un único factor, usad MANTA con type I SS, y el parámetro subset.
# Ruben says: @ Diego Garrido (method author) suggests using SS type = "II" / "III" for DTU analysis

setDT(metadataraw)
metadata_premanta <-  metadataraw[,.(pc1, pc2, sex, population, sample)]
setDT(metadata_premanta)
# gene <- "ENSG00000000460.17"
# contrast <- "CEU"
# ref<-"YRI"
run_manta <- function(gene, contrast, ref) {
  # Print statements for debugging
  print(gene)
  print(contrast)
  print(ref)
  
  metadata_premanta <- metadata_premanta[population%in%c(contrast, ref)]
  setorder(metadata_premanta, sample)
  metadata <- column_to_rownames(metadata_premanta, var="sample")
  # Melt the prepared counts and reshape to wide format
  # meltedprepcounts <- melt(prepared_counts[geneId == gene], id.vars = c("trId", "geneId"), 
  #                          value.name = "rel_ab", variable.name = "sample")
  # wideprepcounts <- as.data.table(dcast(meltedprepcounts, sample ~ trId))[, pop := gsub(".$", "", sample)][pop %in% c(contrast, ref)]
  meltedtrxfilt <- meltedtrx[geneId==gene]
  meltedtrxfilt <- meltedtrxfilt[, pop := gsub(".$", "", sample)][pop %in% c(contrast, ref)]
  meltedtrxfilt <- meltedtrxfilt[, `:=`(trxpergene=NULL, geneId=NULL, trx_counts=NULL, pop=NULL)]
  wideprepcounts <- dcast(meltedtrxfilt, sample ~ trId, value.var="trx_relative_abundance")
  
  # # Step 2: Filter metadata
  # metadatafilt <- metadataraw[, .(pc1, pc2, sex, population, sample)][, population := factor(population, levels = c(ref, contrast))][wideprepcounts, on = "sample"]
  # 
  # Create input data matrix
  setorder(wideprepcounts, sample)
  input_data <- column_to_rownames(as.data.frame(wideprepcounts), var = "sample")
  # input_data <- as.matrix(input_data[, -c(1:4)])
  input_data <- as.matrix(input_data)
  
  # Step 3: Run manta model
  # manta_obj <- manta(
  #   formula = as.formula(paste("input_data ~ ", paste(c(covariates), collapse = " + "))),
  #   data = metadata,
  #   transform = "none", # previously sqrt
  #   type = "III", # previously I
  #   subset = "population"
  # )
  # Step 3: Initialize result table 'res' with 4 columns
  res <- data.table(
    gene = gene,
    contrast = contrast,
    ref = ref,
    p_value = NA  # Placeholder for p_value, will be updated
  )
  
  # Step 4: Run manta model inside tryCatch
  manta_obj <- tryCatch({
    manta(
      formula = as.formula(paste("input_data ~ ", paste(c(covariates), collapse = " + "))),
      data = metadata,
      transform = "none", # previously sqrt
      type = "III",       # previously I
      subset = "population"
    )
  }, error = function(e) {
    # If manta() fails, set p_value to "MantaError" inside res
    message("An error occurred while running manta(): ", e$message)
    res$p_value <- "MantaError"
    return(NULL)  # Return NULL so we know manta failed
  })
  
  # Step 5: Extract p-value if manta ran successfully
  if (!is.null(manta_obj)) {
    res$p_value <- tryCatch({
      # Extract the p-value from the manta object
      manta_obj$aov.tab["population", "Pr(>F)"]
    }, error = function(e) {
      message("An error occurred while extracting p-value: ", e$message)
      return(NA)  # Return NA if p-value extraction fails
    })
  }
  
  return(res)
}



#### PREPARE CONTRASTS
contrasts<- combn(levels(metadataraw$population), m=2, simplify=F)
contrasts <-unlist(lapply(contrasts, \(x) paste0(x[1], "vs", x[2])))
# Create a dataframe of all combinations of genes and contrasts
combinations <- expand.grid(gene = unique(meltedtrx$geneId), contrast = contrasts)
setDT(combinations)
combinations[, ref:=tstrsplit(contrast, "vs")[[2]]][, contrast:=tstrsplit(contrast, "vs")[[1]]]


# Apply manta to each combination
start_time <- Sys.time()
results <- apply(combinations, 1, \(mycomb) {
  run_manta(mycomb['gene'], mycomb['contrast'], mycomb['ref'])})
end_time <- Sys.time()
cat(end_time-start_time)
res <-rbindlist(results)
# multiple testing correction
res[, contrast_name:=paste(contrast, ref, sep="vs")]
res[, fdr:=p.adjust(p_value, method="BH"), by=contrast_name]

res[, sigenes:=sum(fdr<0.05, na.rm=T), by="contrast_name"]




colnames(res)[grep("^gene$", colnames(res))] <- "geneid.v"
colnames(res)[grep("ref", colnames(res))] <- "reference"
colnames(res)[grep("sigenes", colnames(res))] <- "sigDTU"
colnames(res)[grep("p_value", colnames(res))] <- "pval"
fwrite(res, paste0("data/04_DTUres_", TYPE,".tsv"), sep="\t", row.names = F, quote = F)
res <- fread(paste0("data/04_DTUres_", TYPE,".tsv"))


# Create the heatmap

melted_df <-unique(res[,.(contrast, ref, sigenes)])
melted_df <- melted_df[order(contrast, ref)]


# Step 1: Prepare the data by creating a full matrix and filling in the diagonal and lower triangle with NAs
# Combine `contrast` and `ref` to create pairs
dt_full <- rbind(
  melted_df, 
  melted_df[, .(contrast = ref, ref = contrast, sigenes = sigenes)]
)

# Remove duplicated pairs (i.e., where contrast == ref)
dt_full <- unique(dt_full)

# Step 2: Create a complete matrix of all possible combinations of contrast and ref
all_groups <- unique(c(melted_df$contrast, melted_df$ref))
dt_grid <- CJ(contrast = all_groups, ref = all_groups)

# Merge the full data with the grid
dt_heatmap <- merge(dt_grid, dt_full, by = c("contrast", "ref"), all.x = TRUE)

# Remove the lower diagonal and diagonal elements by setting sigenes to NA for one side
dt_heatmap[contrast >= ref, sigenes := NA]


# Step 3: Plot the heatmap using ggplot2
ggplot(dt_heatmap, aes(x = contrast, y = ref, fill = sigenes)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#D2932D", na.value = "grey90", name = "DTU Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", title = nametype)+
  geom_text(aes(label=sigenes, size=sigenes))+
  scale_size_continuous(range=c(3,5))+
  guides(size="none")+
  mytheme
# 3. add additional gene information ----
# 3.1 number of samples (donors) included in gene model
# 3.2 number of annotated & expressed transcripts
# 3.3 most & least contributing transcripts



##### RAQUEL'S FUNCTIONS
# # functions ----
# manta_residuals_fun <- function(formula, data, transform = "sqrt", type = "I", contrasts = NULL, subset = NULL, fit = FALSE){
#   
#   ## Checks
#   transform <- match.arg(transform, c("none", "sqrt", "log"))
#   type <- match.arg(type, c("I", "II", "III"))
#   
#   ## Save call, build model frame, obtain responses
#   cl <- match.call()
#   m <- match(c("formula", "data"), names(cl), 0L)
#   mf <- cl[c(1L, m)]
#   mf$drop.unused.levels <- TRUE
#   mf$na.action = "na.omit"
#   # The rows with at least one NA either in Y or X 
#   # (only considering variables used in the formula) 
#   # will be removed before transforming/centering
#   mf[[1L]] <- quote(stats::model.frame)
#   mf <- eval(mf, parent.frame())               
#   mt <- attr(mf, "terms")
#   response <- model.response(mf, "numeric")
#   
#   ## Checks
#   if (NCOL(response) < 2){
#     stop("The number of response variables should be >= 2")
#   }
#   if (length(attr(mt, "term.labels")) < 1) {
#     stop("The model should contain at least one predictor (excluding the intercept)")
#   }
#   
#   ## Transform and center responses, update model frame
#   if(transform == "none"){
#     Y <- response
#   } else if (transform == "sqrt"){
#     if (any(response < 0)) {
#       stop("'sqrt' transformation requires all response values >= 0")
#     }
#     Y <- sqrt(response)
#   } else if (transform == "log"){
#     if (any(response <= 0)) {
#       stop("'log' transformation requires all response values > 0")
#     }
#     Y <- log(response)
#   }
#   Y <- scale(Y, center = TRUE, scale = FALSE)
#   mf[[1L]] <- Y
#   
#   ## Define contrasts
#   if(is.null(contrasts)){
#     contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
#     dc <- attr(mt, "dataClasses")[-1]
#     contr.list <- lapply(dc, FUN = function(k){
#       # No contrast for quantitative predictors
#       # Sum contrasts for unordered categorical predictors
#       # Polynomial contrasts for ordered categorical predictors
#       contr.type <- switch(k, "factor" = contrasts$unordered,
#                            "ordered" = contrasts$ordered)
#       return(contr.type)
#     })
#     contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
#   } else {
#     contr.list <- contrasts
#   }
#   
#   ## Build model matrix
#   X <- model.matrix(mt, mf, contr.list)
#   
#   ## Fit lm
#   lmfit <- lm.fit(X, Y)
#   return(lmfit$residuals)
# }
# 
# post_hoc_analyses <- function(gene, measurement = "relative_abundance"){
#   
#   # prepare the data to test if the transcript relative abundance is different between AA and EA
#   d <- as.data.frame(manta_residuals_obj[[gene]])
#   d$sample <- rownames(d)
#   data <- melt(d)   
#   colnames(data) <- c("sample", "transcript", "value")
#   data <- merge(data, sample_metadata, by = "sample")
#   
#   # non-parametric Mann-Whitneu U test + multiple test correction
#   p_value <- sapply(unique(data$transcript), function(transcript) 
#     wilcox.test(data[data$transcript==transcript & data$inferred_ancestry == "AA", "value"],
#                 data[data$transcript==transcript & data$inferred_ancestry == "EA", "value"])$p.value)
#   fdr <- p.adjust(p_value, method = "fdr")
#   df <- cbind.data.frame(unique(data$transcript), p_value, fdr)
#   colnames(df)[1] <- "transcript"
#   
#   return(df)
# }


# 
# # 4. calculate residuals ----
# run_manta_residuals <- function(gene){
#   # relative abundances of all gene transcripts
#   input_data <- t(relative_abundances[relative_abundances$trId %in% transcript_annotation[transcript_annotation$ensembl_gene_id == gene, "ensembl_transcript_id"], 3:ncol(relative_abundances)])
#   
#   # calculate residuals (model without variable of interest)
#   residuals <- manta_residuals_fun(
#     formula = as.formula(paste("input_data ~  ", paste(c(covariates[covariates != "inferred_ancestry"]), collapse = " + "))),
#     data = sample_metadata,
#     transform = "sqrt",
#     type = "I")
#   return(residuals)
# }
# 
# # get residual values
# start_time <- Sys.time()
# #manta_residuals_obj <- mclapply(unique(relative_abundances$geneId), function(gene) {print(gene); run_manta_residuals(gene)}, mc.cores = n_cores)
# manta_residuals_obj <- lapply(unique(relative_abundances$geneId), function(gene) {print(gene); run_manta_residuals(gene)})
# names(manta_residuals_obj) <- unique(relative_abundances$geneId)
# end_time <- Sys.time()
# print("Running time for manta residuals")
# print(end_time - start_time)

# # 5. post-hoc non parametric two-tailed Mann Whitney U test ----
# #post_hoc_results <- mclapply(manta_results[manta_results$fdr < 0.05, "ensembl_gene_id"], function(gene) {print(gene); post_hoc_analyses(gene)}, mc.cores = n_cores)
# post_hoc_results <- lapply(manta_results[manta_results$fdr < 0.05, "ensembl_gene_id"], function(gene) {print(gene); post_hoc_analyses(gene)})
# names(post_hoc_results) <- manta_results[manta_results$fdr < 0.05, "ensembl_gene_id"]
# end_time <- Sys.time()
# print("Running time for post-hoc analysis")
# print(end_time - start_time)
# 
# # 6. add post-hoc to manta results ----
# rownames(manta_results) <- manta_results$ensembl_gene_id
# manta_results$number_DU_trascripts <- sapply(manta_results$ensembl_gene_id, function(gene) ifelse(manta_results[gene, "fdr"] < 0.05,
#                                                                                                   sum(post_hoc_results[[gene]]$fdr < 0.05),
#                                                                                                   NA))
# 
# # 7. add information about the DU transcipts ----
# transcripts_DU <- function(gene){
#   df <- post_hoc_results[[gene]][order(post_hoc_results[[gene]]$fdr, decreasing = F),]
#   DU_transcripts <- as.character(df[df$fdr < 0.05, "transcript"])
#   if(length(DU_transcripts)==0){
#     return(NA)
#   }else if(length(DU_transcripts)==1){
#     return(DU_transcripts)
#   }else{
#     DU_transcripts <- paste(DU_transcripts, collapse = "; ")
#     return(DU_transcripts)
#   }
# }
# manta_results$DU_trascripts <- sapply(manta_results$ensembl_gene_id, function(gene) ifelse(manta_results[gene, "fdr"] < 0.05,
#                                                                                            transcripts_DU(gene),
#                                                                                            NA))
# DU_transcripts_biotype <- function(gene){
#   df <- post_hoc_results[[gene]][order(post_hoc_results[[gene]]$fdr, decreasing = F),]
#   DU_transcripts <- as.character(df[df$fdr < 0.05, "transcript"])
#   if(length(DU_transcripts)==0){
#     return(NA)
#   }else{
#     DU_transcripts_biotypes <- sapply(DU_transcripts, function(transcript) transcript_annotation[transcript_annotation$ensembl_transcript_id==transcript, "transcript_biotype"])
#     if(length(DU_transcripts_biotypes)>1){DU_transcripts_biotypes <- paste(DU_transcripts_biotypes, collapse = "; ")}
#     return(DU_transcripts_biotypes)
#   }
# }
# manta_results$DU_transcripts_biotypes <- sapply(manta_results$ensembl_gene_id, function(gene) ifelse(manta_results[gene, "fdr"] < 0.05,
#                                                                                                      DU_transcripts_biotype(gene),
#                                                                                                      NA))
# manta_results <- manta_results[order(manta_results$fdr, decreasing = F),]
# 
# # 8. Save results ----
# 
# saveRDS(relative_abundances,
#         paste0(path_out, tissue, ".relative_abundances.rds"))
# saveRDS(sqtl.seeker.results,
#         paste0(path_out, tissue, ".sqtl.seeker.results.rds"))
# saveRDS(manta_results,
#         paste0(path_out, tissue, ".manta_results.rds"))
# saveRDS(manta_residuals_obj,
#         paste0(path_out, tissue, ".manta_residuals.rds"))
# saveRDS(post_hoc_results,
#         paste0(path_out, tissue, ".post_hoc_results.rds"))
# saveRDS(filtering_summary_chr,
#         paste0(path_out, tissue, ".genes_and_transcripts_kept.rds"))