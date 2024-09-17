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


# Choose covariates
# covariates <- c("sex", "population", "map_reads_generalmap", "captrap_batch", "libprep_batch")
# covariates_ooa <- c("sex", "ooa", "map_reads_generalmap", "captrap_batch", "libprep_batch")
covariates <- c("sex", "population", "map_reads_generalmap", "pc1", "pc2")
covariates_ooa <- c("sex", "ooa", "map_reads_generalmap", "pc1", "pc2")


# Load and parse metadata
metadataraw <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv")
metadataraw <- metadataraw[mixed_samples==FALSE,]
metadataraw[, map_reads_generalmap:=scale(map_reads_generalmap)]
metadataraw[, map_reads_generalmap:=scale(map_reads_generalmap)]
metadataraw <-metadataraw[order(rownames(metadataraw)),]
# Load PCA
mypca <- fread("data/gencode/01_PCA.tsv")
metadataraw <- mypca[metadataraw, on="cell_line_id"]
metadataraw <- column_to_rownames(metadataraw, var="cell_line_id")

metadata <- metadataraw[order(rownames(metadataraw)), c(covariates)]
metadata_ooa <- metadataraw[order(rownames(metadataraw)), c(covariates_ooa)]

# Load and parse counts
counts <- fread("../06_quantification/01_isoquantify/data/gencodev47/gene_counts_matrix.tsv")
counts <- counts[!(geneid.v%in%c("__ambiguous", "__no_feature", "__not_aligned"))]
counts <- column_to_rownames(counts, var="geneid.v")
colnames(counts)  <- gsub(".*_", "", colnames(counts))
counts <- counts[, order(colnames(counts))]
counts <- counts[, colnames(counts)%in%rownames(metadataraw)]


# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 1) >= 5

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

form <- ~ 0 + population + sex + map_reads_generalmap + pc1 +pc2
form_ooa <- ~ 0 + ooa + sex + map_reads_generalmap + pc1 +pc2

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)
vobjDream_ooa <- voomWithDreamWeights(dge, form_ooa, metadata_ooa, BPPARAM = param)

#Since dream uses an estimated degrees of freedom value for each hypothsis test, the degrees of freedom is different for each gene here. Therefore, the t-statistics are not directly comparable since they have different degrees of freedom. In order to be able to compare test statistics, we report z.std which is the p-value transformed into a signed z-score. This can be used for downstream analysis.
# You can also perform a hypothesis test of the difference between two or more coefficients by using a contrast matrix. The contrasts are evaluated at the time of the model fit and the results can be extracted with topTable(). This behaves like makeContrasts() and contrasts.fit() in limma.
#Multiple contrasts can be evaluated at the same time, in order to save computation time. Make sure to inspect your contrast matrix to confirm it is testing what you intend.

## CREATE ALL THE CONTRASTS---------------------
# Create all pairwise combinations
combinations <- combn(unique(metadata$population), 2, simplify = FALSE)

# Generate the contrasts
contrasts <- sapply(combinations, function(pair) {
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
OOAcontrasts <- c(OOAvsAFR= "(populationLWK + populationMPC + populationYRI)/3-(populationITU + populationCEU + populationHAC  + populationAJI + populationPEL)/5")

# Step 6: Combine the original contrasts with the new ones
all_contrasts <- c(original_contrasts, new_contrasts_named, OOAcontrasts)
###--------------------------

# prepare constrast
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
fit <- eBayes(fit)
fit_ooa <- dream(vobjDream_ooa, form_ooa, metadata_ooa, L_ooa)
fit_ooa <- eBayes(fit_ooa)
errors <- attr(fit_ooa, 'errors')
print(errors)

# CHECK RESULTS........................................................

# get names of available coefficients and contrasts for testing
colnames(fit)

# extract results from first contrast
results_list <- lapply(colnames(fit_ooa)[1], function(contrast) {
  topTable(fit_ooa, coef = contrast, number = Inf, adjust.method = "BH")
})
names(results_list) <- colnames(fit)[1:37]

#
deg_threshold <- function(result, p_value_cutoff = 0.05, logFC_cutoff = 0.5) {
  degs <- result[result$adj.P.Val < p_value_cutoff & abs(result$logFC) >= logFC_cutoff, ]
  return(degs)
}

# Apply the threshold function to each contrast result
deg_list <- lapply(results_list, deg_threshold)
pairwise_data <- lapply(deg_list, nrow)[1:(length(names(pairwise_data))-9)]



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





mygenes <-gsub("\\..*","",rownames(deg_list$`(populationLWK + populationMPC + populationYRI)/3 - (populationITU + populationCEU + populationHAC + populationAJI + populationPEL)/5`))
fwrite(as.data.frame(mygenes),"/home/pclavell/Downloads/mygenes.tsv", quote =F, row.names = F)
mybck <- gsub("\\..*","",rownames(dge$counts))
fwrite(as.data.frame(mybck),"/home/pclavell/Downloads/mybck.tsv", quote =F, row.names = F)
