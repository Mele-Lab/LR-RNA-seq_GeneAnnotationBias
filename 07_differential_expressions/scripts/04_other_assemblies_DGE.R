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


TYPES <-c("gencode", "pantrx")
ASSEMBLIES <- c("HG002_maternal", "HG02717_maternal")

TYPE <- TYPES[2]
ASSEMBLY <- ASSEMBLIES[1]
# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,"_genelvl.tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47"
  counts <- fread(paste0("../../novelannotations/quantifications/v47_personal_kallisto_quant/", ASSEMBLY,"/matrix.abundance.tsv"))
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER"
  counts <- fread(paste0("../../novelannotations/quantifications/personal_kallisto_quant/", ASSEMBLY,"/matrix.abundance.tsv"))
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
if(TYPE=="pantrx"){
  fwrite(rownames_to_column(counts, var="geneid.v"), 
         paste0("../../novelannotations/quantifications/personal_kallisto_quant/",ASSEMBLY,"/matrix.abundance.genelevel.tsv"), row.names = F, quote = F, sep="\t")
}else if(TYPE=="gencode"){
  fwrite(rownames_to_column(counts, var="geneid.v"), 
         paste0("../../novelannotations/quantifications/v47_personal_kallisto_quant/",ASSEMBLY,"/matrix.abundance.genelevel.tsv"), row.names = F, quote = F, sep="\t")
}
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

resultsdt <- rbindlist(lapply(results_list, rownames_to_column, var="geneid.v"), idcol="contrast")
resultsdt$comparison <- tstrsplit(resultsdt$contrast, "vs")[[1]]
resultsdt$reference <- tstrsplit(resultsdt$contrast, "vs")[[2]]

resultsdt[, upDEG := sum(logFC >= 0.5 & adj.P.Val < 0.05), by = "contrast"]
resultsdt[, downDEG := sum(logFC <= -0.5 & adj.P.Val < 0.05), by = "contrast"]


colnames(resultsdt)[grep("adj.P.Val", colnames(resultsdt))] <- "fdr"
colnames(resultsdt)[grep("contrast", colnames(resultsdt))] <- "contrast_name"
colnames(resultsdt)[grep("comparison", colnames(resultsdt))] <- "contrast"
colnames(resultsdt)[grep("P.Value", colnames(resultsdt))] <- "pval"
fwrite(resultsdt, paste0("data/other_assemblies/02_DEGres_", ASSEMBLY,"_", TYPE,".tsv"), quote = F, sep = "\t", row.names = F)




resultsdt <- fread(paste0("data/other_assemblies/02_DEGres_", ASSEMBLY,"_",TYPE,".tsv"))

# resultsdt[, upDEG:=-upDEG]

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
# 
# #
# deg_threshold <- function(result, p_value_cutoff = 0.05, logFC_cutoff = 0.5) {
#   up <- nrow(result[result$adj.P.Val < p_value_cutoff & result$logFC >= logFC_cutoff, ])
#   down <- nrow(result[result$adj.P.Val < p_value_cutoff & result$logFC <= -logFC_cutoff, ])
#   return(c(up, down))
# }
# 
# # Apply the threshold function to each contrast result
# deg_list <- lapply(results_list, deg_threshold)
# # Create an empty data.table
# result <- data.table()
# 
# # Iterate over each item in the list
# for (name in names(deg_list)) {
#   # Split the name at "vs"
#   parts <- strsplit(name, "vs")[[1]]
#   
#   # Extract the up and down values
#   values <- deg_list[[name]]
#   
#   # Create a new row and bind it to the result
#   result <- rbind(result, data.table(
#     up = values[1],
#     down = values[2],
#     before_vs = trimws(parts[1]),
#     after_vs = trimws(parts[2])
#   ))
# }
# result[, down:=-down]
# Reorder columns to desired format
result <- unique(resultsdt[, .(contrast, reference, upDEG, downDEG)])
# setcolorder(result, c("before_vs", "after_vs", "up", "down"))

# Create a unique list of contrasts
contrasts <- unique(c(result$contrast, result$reference))

# Initialize a matrix with NA values
matrix_size <- length(contrasts)
result_matrix <- matrix(NA, nrow = matrix_size, ncol = matrix_size)

# Set row and column names
rownames(result_matrix) <- contrasts
colnames(result_matrix) <- contrasts

# Fill the matrix with up and down values
for (i in 1:nrow(result)) {
  before <- result$contrast[i]
  after <- result$reference[i]
  up_value <- result$upDEG[i]
  down_value <- result$downDEG[i]
  
  # Assign values to the matrix
  result_matrix[before, after] <- up_value      # Upregulated in upper triangle
  result_matrix[after, before] <- down_value     # Downregulated in lower triangle
}

# Melt the data frame to long format
melted_df <- melt(result_matrix, value.name = "DEGs")

# Remove NA values for better visualization
melted_df <- melted_df[!is.na(melted_df$DEGs), ]
colnames(melted_df) <- c( "Contrast","Ref", "DEGs")

melted_df<-rbind.data.frame(melted_df, data.frame(cbind("Contrast"=levels(melted_df$Contrast),
                                                        "Ref"=levels(melted_df$Contrast), 
                                                        "DEGs"=rep(NA, 8))))
setDT(melted_df)
melted_df[,DEGs:=as.numeric(DEGs)]


### PLOT
TYPES <-c("gencode", "pantrx")
ASSEMBLIES <- c("HG002_maternal", "HG02717_maternal")

TYPE <- TYPES[2]
if(TYPE=="gencode"){
  nametype <- "GENCODE"
}else if(TYPE=="pantrx"){
  nametype <- "PODER"
}
ASSEMBLY <- ASSEMBLIES[2]
melted_df <- fread(paste0("data/other_assemblies/02_DEGres_", ASSEMBLY, "_",TYPE, "_sigGenesMatrix.tsv"))

# Create the heatmap
ggplot(melted_df, aes(x = Ref, y = Contrast, fill = DEGs,)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#A0000D", midpoint = 0,limits=c(0,210), na.value = "#9AC7CB") +
  theme_minimal() +
  labs(title=paste0(ASSEMBLY, " with ",nametype),
       x = "Against Reference",
       y="Upregulated Genes in",
       fill = "# DEGs") +
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title=element_text(face="bold", hjust=0.5)) +
  geom_text(data=melted_df[!is.na(DEGs)], aes(label=abs(DEGs), size=abs(DEGs)),  na.rm =T)+
  coord_fixed()+
  scale_size_continuous(range=c(4.5,8)*0.35)+
  guides(size="none")+
  scale_y_discrete(limits = levels(melted_df$Contrast))
ggsave(paste0("../10_figures/01_plots/supp/28_deg_assemblies/heatmap_DGE_",ASSEMBLY, "_",TYPE,".pdf"), dpi=700, width = 2.75, height = 2.75,  units = "in")



fwrite(melted_df, paste0("data/other_assemblies/02_DEGres_", ASSEMBLY, "_",TYPE, "_sigGenesMatrix.tsv"), sep = "\t", quote = F, row.names = F)

