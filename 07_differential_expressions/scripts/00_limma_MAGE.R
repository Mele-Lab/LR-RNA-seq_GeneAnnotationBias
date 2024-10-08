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
library(limma)
library("variancePartition")
library("edgeR")
library("BiocParallel")

# load and parse data
counts <- fread("../../../Data/MAGE/MAGE.v1.0.data/global_trend_results/global_expression_trends/expression.counts.csv")
metadata <- fread("../../../Data/MAGE/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt")

counts <- column_to_rownames(counts, var="gene")
counts <- counts[, order(colnames(counts))]
metadata <- metadata[order(metadata$sample_coriellID),]
metadata <- column_to_rownames(metadata, var="sample_coriellID")
# filter genes by number of counts
isexpr <- rowSums(cpm(counts) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
dge <- dge[1:1000, ]

metadata$batch <- factor(metadata$batch)
metadata$population <- factor(metadata$population)
metadata$sex <- factor(metadata$sex)

design <- model.matrix(~ 0+population+sex+batch, metadata) # Assuming 'group' has your experimental conditions


#dupcor <- duplicateCorrelation(vobj_tmp, design, block = metadata$Individual)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj <- voom(dge,  design, plot = FALSE)
#vobj <- voom(dge, design, plot = FALSE, block = metadata$Individual, correlation = dupcor$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene,
#dupcor <- duplicateCorrelation(vobj, design, block = metadata$Individual)


## CREATE ALL THE CONTRASTS---------------------
# Create all pairwise combinations
combinations <- combn(unique(metadata$population), 2, simplify = FALSE)

# Generate the contrasts
contrasts <- sapply(combinations, function(pair) {
  paste0(pair[1], "-", pair[2])
})
contrasts <- sapply(combinations, function(pair) {
  paste0(pair[1], "vs", pair[2], " = \"population", pair[1], " - population", pair[2], "\"")
})
mycon <-paste(contrasts, collapse=", ")
# Convert the contrasts into a named vector
contrast_list <- eval(parse(text = paste0("list(", mycon, ")")))

# # Step 3: Create new contrasts where each population is compared against the average of all others
# new_contrasts <- sapply(populations, function(pop) {
#   other_pops <- setdiff(populations, pop)
#   paste0("population", pop, " - (", paste0("population",other_pops, collapse = " + "), ")/", length(other_pops))
# })
# 
# # Step 4: Convert new contrasts to a named list
# new_contrasts_named <- setNames(new_contrasts, paste0(gsub(" .*", "", new_contrasts), "vsALL"))
# 
# # Step 5: Convert the original contrast string to a list
# original_contrasts <- eval(parse(text = paste0("list(", mycon, ")")))
# 
# # Step 6: Combine the original contrasts with the new ones
# all_contrasts <- c(original_contrasts, new_contrasts_named)
# ###--------------------------

# make contrasts
L <- makeContrasts(contrasts = contrast_list, levels=design)


# But this step uses only the genome-wide average for the random effect
fit <- lmFit(vobj, design)
fit2 <- contrasts.fit(fit, L)
#fitDupCor <- lmFit(vobj, design, block = metadata$Individual, correlation = dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics
fit2 <- eBayes(fitDupCor)

# extract results from first contrast
results_list <- lapply(colnames(fit2)[1:37], function(contrast) {
  variancePartition::})
names(results_list) <- colnames(fit)[1:37]
topTable(fit, coef ="populationACB-populationASW" , number = Inf, adjust.method = "BH")
################333




# make this vignette faster by analyzing a subset of genes
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# QUE SIGNIFICA AQUEST 0+??????????????????????????????????????????????????
#form <- ~ 0 + population + sex + map_reads_generalmap + (1|captrap_batch) + (1|libprep_batch)
#form_ooa <- ~ 0 + ooa + sex + map_reads_generalmap + (1|captrap_batch) + (1|libprep_batch)

form <- ~ 0+ population+ sex + (1|batch)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)


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

# Step 6: Combine the original contrasts with the new ones
all_contrasts <- c(original_contrasts, new_contrasts_named)
###--------------------------

# prepare constrast
#L <- makeContrastsDream(form, metadata, contrasts = "populationCEU-populationHAC")
L <- makeContrastsDream(form, metadata,
                        contrasts = all_contrasts)

# Visualize contrast matrix
plotContrasts(L)

# fit dream model with contrasts

fit <- dream(vobjDream, form, metadata, L)
errors <- attr(fit, 'errors')
print(errors)
fit <- variancePartition::eBayes(fit)

# CHECK RESULTS........................................................

# get names of available coefficients and contrasts for testing
colnames(fit)

# extract results from first contrast
results_list <- lapply(colnames(fit)[1:325], function(contrast) {
  variancePartition::topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
})
names(results_list) <- colnames(fit)[1:325]



# plot
#
deg_threshold <- function(result, p_value_cutoff = 0.05, logFC_cutoff = 0.5) {
  degs <- result[result$adj.P.Val < p_value_cutoff & abs(result$logFC) >= logFC_cutoff, ]
  return(degs)
}

# Apply the threshold function to each contrast result
deg_list <- lapply(results_list, deg_threshold)
pairwise_data <- lapply(deg_list, nrow)



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



