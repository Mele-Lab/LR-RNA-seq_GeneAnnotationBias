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
LEVEL <-"genelvl"

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
metadata <- metadata[merged_run_mode==TRUE & mixed_samples==FALSE]
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
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





## plot pca
# Extract PCA coordinates
res.pca <- mypca
pca_coords <- as.data.frame(mypca$ind$coord)  # PC coordinates for individuals
pca_coords <- rownames_to_column(pca_coords, var="cell_line_id")
# Add your factor variables to the PCA data frame
# Assuming 'factor_data' contains a column named 'Group' that categorizes the samples
pca_coords <- metadata[, .(captrap_batch, population, sex, libprep_batch, cell_line_id, sample)][pca_coords, on="cell_line_id"]
# Plot PC1 vs PC2 with ggplot2
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = captrap_batch)) +
  geom_point(size = 1.5) +  # Adjust point size if needed
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)"),
    color = "CapTrap\nBatch"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  mytheme
ggsave("../10_figures/01_plots/supp/23_variance_part/scatter_PCA_PODER_GENE_captrapbatch.pdf", dpi=700, width = 3, height = 2.7,  units = "in")
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = population)) +
  geom_point(size = 1.5) +  # Adjust point size if needed
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)"),
    color = "Population"
  ) +
  scale_color_manual(values=popcols)+
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  mytheme
ggsave("../10_figures/01_plots/supp/23_variance_part/scatter_PCA_PODER_GENE_population.pdf", dpi=700, width = 3, height = 2.7,  units = "in")
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = libprep_batch)) +
  geom_point(size = 1.5) +  # Adjust point size if needed
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)"),
    color = "Library\nPreparation\nBatch"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  mytheme
ggsave("../10_figures/01_plots/supp/23_variance_part/scatter_PCA_PODER_GENE_libprep.pdf", dpi=700, width = 3, height = 2.7,  units = "in")
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = sex)) +
  geom_point(size = 1.5) +  # Adjust point size if needed
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)"),
    color = "Sex"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))+
  mytheme
ggsave("../10_figures/01_plots/supp/23_variance_part/scatter_PCA_PODER_GENE_sex.pdf", dpi=700, width = 3, height = 2.7,  units = "in")



# Now asses variation explained
eta2_matrix <- mypca[["quali.sup"]][["eta2"]]  # Extract eta-squared values
eta2_matrix <- rownames_to_column(as.data.frame(eta2_matrix), var="factor")
setDT(eta2_matrix)
# Convert eta2_matrix into a tidy format for ggplot2
eta2_long <- melt(eta2_matrix, variable.name = "PC", value.name = "eta_squared")

# Plot the eta2 matrix using ggplot2
ggplot(eta2_long, aes(x = PC, y = factor, fill = eta_squared)) +
  geom_tile(color = "white") +  # Add gridlines between tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +  # Adjust colors
  labs( x = "Principal Components",
    y = "",
    fill = "Variance\nExplained") +
  mytheme+
  theme( axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_discrete(labels=c("trizol_batchB"="Cell Collection 2",
                            "trizol_batchA"="Cell Collection 1",
                            "sex"="Sex",
                            "rna_extraction_batchB"="RNA Extraction 2",
                            "rna_extraction_batchA"="RNA Extraction 1",
                            "population"="Population",
                            "libprep_batch"="Library\nPreparation Batch",
                            "captrap_batch"="CapTrap Batch"))+
  geom_text(aes(label=round(eta_squared, digits = 2), size=eta_squared))+
  scale_size_continuous(range=c(4,7)*0.35)+
  guides(size="none")
ggsave("../10_figures/01_plots/supp/23_variance_part/heatmap_PCA_variance.pdf", dpi=700, width = 4, height = 2.7,  units = "in")
