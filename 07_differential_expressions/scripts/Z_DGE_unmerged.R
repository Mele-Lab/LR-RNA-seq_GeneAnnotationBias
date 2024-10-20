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

library("variancePartition")
library("edgeR")
library("BiocParallel")

##### PREPARE DATA ########################################################################
# Load and parse metadata
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadataraw <- metadataraw[mixed_samples==FALSE,]
metadataraw <- metadataraw[unmerged_run_mode==TRUE]

# prepare new names vector
samplesnames <- metadataraw$sample
names(samplesnames) <- metadataraw$quantification_id
samplesnames <- c(samplesnames, "transcript_id"="transcriptid.v", "geneid.v"="geneid.v")
# adapt metadata to build design matrix
metadataraw <-metadataraw[order(rownames(metadataraw)),]
metadataraw <- column_to_rownames(metadataraw, var="sample")

covariates <- c("sex", "population","captrap_new_batch", "libprep_batch", "cell_line_id")
covariates_ooa <- c("sex", "ooa","captrap_new_batch", "libprep_batch", "cell_line_id")
metadata <- metadataraw[order(rownames(metadataraw)), c(covariates)]
metadata_ooa <- metadataraw[order(rownames(metadataraw)), c(covariates_ooa)]
metadata[, population:=relevel(population, ref="CEU")]

# Load and parse counts
counts <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.tsv")
unmergedcounts <- fread("../../novelannotations/v47_unmerged_kallisto_quant/matrix.abundance.tsv")
counts <- merge(counts, unmergedcounts, by = "transcript_id", all = TRUE)
annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")

# Sum transcript counts per gene
counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
counts <- counts[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
# remove merged samples
colnames(counts) <- gsub("_1$", "", colnames(counts))
counts <- counts[, .SD, .SDcols=colnames(counts)%in%c(names(samplesnames), "geneid.v")]
setnames(counts, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)

# prepare datatable for analysis
counts <- column_to_rownames(counts, var="geneid.v")
counts <- counts[, order(colnames(counts))]


# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 1) >= 5

##### PROCESS DATA ###########################################################################
# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# QUE SIGNIFICA AQUEST 0+??????????????????????????????????????????????????
#form <- ~ 0 + population + sex + map_reads_generalmap + (1|captrap_batch) + (1|libprep_batch)
#form_ooa <- ~ 0 + ooa + sex + map_reads_generalmap + (1|captrap_batch) + (1|libprep_batch)

form <- ~ population + sex + (1|captrap_new_batch) + (1|libprep_batch) + (1|cell_line_id)
form_ooa <- ~ ooa + sex + (1|captrap_new_batch) + (1|libprep_batch) + (1|cell_line_id)
############################
# estimate weights using linear mixed model of dream
v <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

# Fit linear model
fit <- lmFit(v, design)

# Define contrasts
contrast.matrix <- makeContrasts(TreatmentVsControl = ConditionTreatment - ConditionControl, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract top differentially expressed genes
results <- topTable(fit2, coef = "TreatmentVsControl", number = Inf)



#### CREATE ALL THE CONTRASTS ############################################################
# Create all pairwise combinations
combinations <- expand.grid(unique(metadata$population), unique(metadata$population))

# Generate the contrasts
contrasts <- apply(combinations, 1,function(pair) {
  paste0(pair[1], "vs", pair[2], " = \"population", pair[1], " - population", pair[2], "\"")
})

mycon <-paste(contrasts, collapse=", ")
# Convert the contrasts into a named vector
contrast_list <- eval(parse(text = paste0("list(", mycon, ")")))

populations <- unique(metadata$population)
# Step 3: Create new contrasts where each population is compared against the average of all others
new_contrasts <- sapply(populations, function(pop) {
  other_pops <- setdiff(populations, pop)
  paste0("population", pop, " - (", paste0("population",other_pops, collapse = " + "), ")/", length(other_pops))
})

# Step 4: Convert new contrasts to a named list
new_contrasts_named <- setNames(new_contrasts, paste0(gsub(" .*", "", new_contrasts), "vsALL"))

# Step 5: Convert the original contrast string to a list
original_contrasts <- eval(parse(text = paste0("list(", mycon, ")")))
OOAcontrasts <- c(AFRvsOOA= "(populationLWK + populationMPC + populationYRI)/3-(populationITU + populationCEU + populationHAC  + populationAJI + populationPEL)/5")

# Step 6: Combine the original contrasts with the new ones
all_contrasts <- c(original_contrasts, new_contrasts_named, OOAcontrasts)


# prepare constrast
#L <- makeContrastsDream(form, metadata, contrasts = "populationCEU-populationHAC")
L <- makeContrastsDream(form, metadata,
                        contrasts = all_contrasts)
L_ooa <- makeContrastsDream(form_ooa, metadata_ooa,
                            contrasts = c(OOAvsAFR= "ooaAFR - ooaOOA"))
# Visualize contrast matrix
plotContrasts(L)
plotContrasts(L_ooa)

# fit dream model with contrasts

fit <- dream(vobjDream, form, metadata, L)
errors <- attr(fit, 'errors')
print(errors)
fit <- variancePartition::eBayes(fit)
fit_ooa <- dream(vobjDream_ooa, form_ooa, metadata_ooa, L_ooa)
fit_ooa <-  variancePartition::eBayes(fit_ooa)
errors <- attr(fit_ooa, 'errors')
print(errors)

# CHECK RESULTS........................................................

# get names of available coefficients and contrasts for testing
colnames(fit)

# extract results from first contrast
results_list <- lapply(colnames(fit)[1:37], function(contrast) {
  variancePartition::topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
})
names(results_list) <- colnames(fit)[1:37]
variancePartition::topTable(fit, coef ="populationMPC" , number = Inf, adjust.method = "BH")

# ##### TEMPORARY volcano -----------
# for(x in 1:37){
#   p<-ggplot(data.table(results_list[[x]])[, sig:=ifelse(adj.P.Val<0.05, "sig", "ns")], aes(col=sig, x=logFC, y=-adj.P.Val))+
#     geom_point(alpha=0.5, size=0.5)+
#     geom_hline(yintercept=-0.05, color="black", linetype="dashed")+
#     mytheme+
#     scale_color_manual(values=rev(c("darkred", "darkgrey")))+
#     ggtitle(names(results_list)[x])
#   print(p)}

####



#
deg_threshold <- function(result, p_value_cutoff = 0.05, logFC_cutoff = 0.5) {
  degs <- result[result$adj.P.Val < p_value_cutoff & abs(result$logFC) >= logFC_cutoff, ]
  return(degs)
}

# Apply the threshold function to each contrast result
deg_list <- lapply(results_list, deg_threshold)
pairwise_data <- lapply(deg_list, nrow)[1:(length(names(original_contrasts))-9)]



# Extract unique population names
populations <- unique(unlist(strsplit(names(pairwise_data), " - ")))
populations <- trimws(populations)  # Trim any leading or trailing whitespace

# Create an empty matrix
pairwise_matrix <- matrix(as.numeric(0), nrow = length(populations), ncol = length(populations), 
                          dimnames = list(populations, populations))
pairwise_matrix <- as.data.frame(pairwise_matrix)
# Fill the matrix with the data
for (i in seq_along(pairwise_data)) {
  pops <- unlist(strsplit(names(pairwise_data)[i], " - "))
  pops <- trimws(pops)  # Trim any leading or trailing whitespace
  
  # Correctly assign the value if both populations are found in the matrix dimensions
  if (all(pops %in% populations)) {
    pairwise_matrix[pops[1], pops[2]] <- as.numeric(pairwise_data[i])
  }
}

pairwise <- rownames_to_column(pairwise_matrix, var="contrast")
pairwiselong <- melt(pairwise,id.vars ="contrast" , value.name ="degs")
pairwiselong$contrast <- gsub("population", "", pairwiselong$contrast)
pairwiselong$variable <- gsub("population", "", pairwiselong$variable)

library(ggplot2)

# Ensure that `contrast` and `variable` are factors with ordered levels
pairwiselong$contrast <- factor(pairwiselong$contrast, levels = unique(pairwiselong$contrast))
pairwiselong$variable <- factor(pairwiselong$variable, levels = unique(pairwiselong$variable))

# Create the heatmap
ggplot(pairwiselong, aes(x = contrast, y = variable, fill = degs)) +
  geom_tile() +
  theme_minimal()+
  scale_fill_gradient(low = "white", high = "red") # Optional: Use a color scale suitable for heatmaps



# 
# 
# mygenes <-gsub("\\..*","",rownames(deg_list$`(populationLWK + populationMPC + populationYRI)/3 - (populationITU + populationCEU + populationHAC + populationAJI + populationPEL)/5`))
# fwrite(as.data.frame(mygenes),"/home/pclavell/Downloads/mygenes.tsv", quote =F, row.names = F)
# mybck <- gsub("\\..*","",rownames(dge$counts))
# fwrite(as.data.frame(mybck),"/home/pclavell/Downloads/mybck.tsv", quote =F, row.names = F)

#### run enrichments
resultsdt <- rbindlist(results_list, idcol="contrast")
resultsdt <- unique(resultsdt[,.(adj.P.Val, geneid.v)])[, geneid:=gsub("\\..*", "", geneid.v)]

passdata <-unique(data[filter=="pass" & sample=="YRI5", adj.P.Val:=fdr][, .(adj.P.Val, geneid.v)])
passdata[, geneid:=gsub("\\..*", "", geneid.v)]
ora_res <- run_ora(unique(passdata[, .(geneid, adj.P.Val)]), db=dbs, keyType="ENSEMBL")
res <-extract_ora_res_individual(ora_res)
ggplot(res, aes(x=Count, y=reorder(Description, Count), col=qvalue))+
  geom_point()+
  mytheme
eval(res$GeneRatio)