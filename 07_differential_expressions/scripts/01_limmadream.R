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
TYPE<-"gencode"
library("variancePartition")
library("edgeR")
library("BiocParallel")

# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")
# covariates_ooa <- c("sex", "ooa", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,".tsv"))
if(TYPE=="GENCODE"){
  counts <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
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
metadataraw[, ooa:=factor(ooa, levels=c("OOA", "AFR"))]

# prepare new names vector
samplesnames <- metadataraw$sample
names(samplesnames) <- metadataraw$quantification_id
samplesnames <- c(samplesnames, "transcript_id"="transcriptid.v", "geneid.v"="geneid.v")
metadataraw <- column_to_rownames(metadataraw, var="sample")
metadata <- metadataraw[order(rownames(metadataraw)), c(covariates)]
metadata_ooa <- metadataraw[order(rownames(metadataraw)), c(covariates_ooa)]
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

# make this vignette faster by analyzing a subset of genes
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

form <- ~population+sex+pc1+pc2
# form_ooa <- ~ooa+sex+pc1+pc2

design <- model.matrix(form, metadata)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)
# vobjDream_ooa <- voomWithDreamWeights(dge, form_ooa, metadata_ooa, BPPARAM = param)

# PREPARE CONTRASTS #########################################################################........................................................
# combinations <- expand.grid(levels(metadata$population), 
#                             levels(metadata$population))
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
# fit_ooa <- dream(vobjDream_ooa, form_ooa, metadata_ooa, L_ooa)
# fit_ooa <-  variancePartition::eBayes(fit_ooa)

# extract results from first contrast
results_list <- lapply(colnames(fit), function(contrast) {
  variancePartition::topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
})
names(results_list) <- colnames(fit)



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
  up <- nrow(result[result$adj.P.Val < p_value_cutoff & result$logFC >= logFC_cutoff, ])
  down <- nrow(result[result$adj.P.Val < p_value_cutoff & result$logFC <= -logFC_cutoff, ])
  return(c(up, down))
}

# Apply the threshold function to each contrast result
deg_list <- lapply(results_list, deg_threshold)
# Create an empty data.table
result <- data.table()

# Iterate over each item in the list
for (name in names(deg_list)) {
  # Split the name at "vs"
  parts <- strsplit(name, "vs")[[1]]
  
  # Extract the up and down values
  values <- deg_list[[name]]
  
  # Create a new row and bind it to the result
  result <- rbind(result, data.table(
    up = values[1],
    down = values[2],
    before_vs = trimws(parts[1]),
    after_vs = trimws(parts[2])
  ))
}
result[, down:=-down]
# Reorder columns to desired format
setcolorder(result, c("before_vs", "after_vs", "up", "down"))

# Create a unique list of contrasts
contrasts <- unique(c(result$before_vs, result$after_vs))

# Initialize a matrix with NA values
matrix_size <- length(contrasts)
result_matrix <- matrix(NA, nrow = matrix_size, ncol = matrix_size)

# Set row and column names
rownames(result_matrix) <- contrasts
colnames(result_matrix) <- contrasts

# Fill the matrix with up and down values
for (i in 1:nrow(result)) {
  before <- result$before_vs[i]
  after <- result$after_vs[i]
  up_value <- result$up[i]
  down_value <- result$down[i]
  
  # Assign values to the matrix
  result_matrix[before, after] <- up_value      # Upregulated in upper triangle
  result_matrix[after, before] <- down_value     # Downregulated in lower triangle
}

# Melt the data frame to long format
melted_df <- melt(result_matrix, value.name = "DEGs")

# Remove NA values for better visualization
melted_df <- melted_df[!is.na(melted_df$DEGs), ]
colnames(melted_df) <- c( "Contrast","Ref", "DEGs")
# Create the heatmap
ggplot(melted_df, aes(x = Contrast, y = Ref, fill = DEGs)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#A0000D", midpoint = 0, na.value = "grey") +
  theme_minimal() +
  labs(title="PODER annotation",
       y = "Comparison (vs Ref)",
       x = "Reference",
       fill = "# DEGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=abs(DEGs), size=abs(DEGs)))+
  coord_fixed()+
  scale_size_continuous(range=c(3,5))+
  guides(size="none")+
  mytheme










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