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
isexpr <- rowSums(cpm(counts) > 1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

metadata$batch <- factor(metadata$batch)
metadata$population <- factor(metadata$population)
metadata$population <- factor(metadata$population, levels=c("CEU",levels(metadata$population)[!grepl("CEU",levels(metadata$population))]))
metadata$sex <- factor(metadata$sex)


param <- SnowParam(4, "SOCK", progressbar = TRUE)

form <- ~population+sex+(1|batch)
form2 <- ~population+sex+batch
design <- model.matrix(form2, metadata)
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
                        list(
                          "ACBvsCEU" = "populationACB",
                          "ASWvsCEU" = "populationASW",
                          "BEBvsCEU" = "populationBEB",
                          "CDXvsCEU" = "populationCDX",
                          "CHBvsCEU" = "populationCHB",
                          "CHSvsCEU" = "populationCHS",
                          "CLMvsCEU" = "populationCLM",
                          "ESNvsCEU" = "populationESN",
                          "FINvsCEU" = "populationFIN",
                          "GBRvsCEU" = "populationGBR",
                          "GIHvsCEU" = "populationGIH",
                          "GWDvsCEU" = "populationGWD",
                          "IBSvsCEU" = "populationIBS",
                          "ITUvsCEU" = "populationITU",
                          "JPTvsCEU" = "populationJPT",
                          "KHVvsCEU" = "populationKHV",
                          "LWKvsCEU" = "populationLWK",
                          "MSLvsCEU" = "populationMSL",
                          "MXLvsCEU" = "populationMXL",
                          "PELvsCEU" = "populationPEL",
                          "PJLvsCEU" = "populationPJL",
                          "PURvsCEU" = "populationPUR",
                          "STUvsCEU" = "populationSTU",
                          "TSIvsCEU" = "populationTSI",
                          "YRIvsCEU" = "populationYRI"
                        ))
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

fwrite(resultsdt, "07_differential_expressions/data/02_DEGres_MAGE.tsv", quote = F, sep = "\t", row.names = F)
resultsdt <- fread("07_differential_expressions/data/MAGE_validations/02_DEGres_MAGE.tsv")

resultsdt[, upDEG:=-upDEG]


# Reorder columns to desired format
result <- unique(resultsdt[, .(comparison, reference, upDEG, downDEG)])
# setcolorder(result, c("before_vs", "after_vs", "up", "down"))

# Create a unique list of contrasts
contrasts <- unique(c(result$comparison, result$reference))

# Initialize a matrix with NA values
matrix_size <- length(contrasts)
result_matrix <- matrix(NA, nrow = matrix_size, ncol = matrix_size)

# Set row and column names
rownames(result_matrix) <- contrasts
colnames(result_matrix) <- contrasts

# Fill the matrix with up and down values
for (i in 1:nrow(result)) {
  before <- result$comparison[i]
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

melted_df<-rbind.data.frame(melted_df, data.frame(cbind("Contrast"=levels(metadata$population),
                                                        "Ref"=levels(metadata$population), 
                                                        "DEGs"=rep(NA, 8))))
setDT(melted_df)
melted_df[,DEGs:=as.numeric(DEGs)]
# generate annotations
# Define which populations are in Group A and Group B
eur <- c("CEU", "FIN", "GBR", "IBS", "TSI")  # Example populations in Group A
afr <- c("YRI", "LWK", "ESN", "GWD", "MSL", "ASW", "ACB" )  # Example populations in Group B
amr <- c("MXL", "CLM", "PEL", "PUR")
sas <- c("PJL", "GIH", "ITU", "STU", "BEB")
eas <- c("CHB", "CHS","JPT", "KHV", "CDX")

# Create group annotations for x-axis (Ref)
x_annotations <- data.frame(
  Ref = unique(melted_df$Ref),
  Superpopulation = ifelse(unique(melted_df$Ref) %in% eur, "EUR",
                           ifelse(unique(melted_df$Ref) %in% afr, "AFR",
                                  ifelse(unique(melted_df$Ref) %in% amr, "AMR",
                                         ifelse(unique(melted_df$Ref) %in% sas, "SAS",
                                                ifelse(unique(melted_df$Ref) %in% eas, "EAS", NA))))))

superpopcols <- c("#A58F39", "#BF4B4B", "#3A90A8", "#3E7F42", "#6444AD")
names(superpopcols)<- c("AFR", "AMR", "EUR", "EAS", "SAS")
# Create group annotations for y-axis (Contrast)
y_annotations <- data.frame(
  Contrast = unique(melted_df$Contrast),
  Superpopulation = ifelse(unique(melted_df$Contrast) %in% eur, "EUR",
                           ifelse(unique(melted_df$Contrast) %in% afr, "AFR",
                                  ifelse(unique(melted_df$Contrast) %in% amr, "AMR",
                                         ifelse(unique(melted_df$Contrast) %in% sas, "SAS",
                                                ifelse(unique(melted_df$Contrast) %in% eas, "EAS", NA))))))

setDT(y_annotations)
setDT(x_annotations)

meltednew <- y_annotations[melted_df[order(Contrast, Ref)], on="Contrast"][order(Superpopulation)][, Contrast:=factor(Contrast, levels=unique(Contrast))]
meltednew <- x_annotations[meltednew, on="Ref"][order(Superpopulation)][, Ref:=factor(Ref, levels=unique(Ref))]

meltednew <- fread(paste0("07_differential_expressions/data/MAGE_validations/02_DEGres_MAGE_sigGenesMatrix.tsv"))
meltednew[Ref==Contrast, DEGs:=NA]
# Create the heatmap
ggplot(meltednew, aes(x = Ref, y = Contrast, fill = abs(DEGs),)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#A0000D", midpoint = 0, na.value = "#9AC7CB") +
  theme_minimal() +
  labs(title="MAGE data - GENCODEv38",
       x = "Against Reference",
       y="Upregulated Genes in",
       fill = "# DEGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=melted_df[!is.na(DEGs)], aes(label=abs(DEGs), size=abs(DEGs)),  na.rm =T)+
  coord_fixed()+
  scale_size_continuous(range=c(3.5,7)*0.35)+
  guides(size="none")+
  mytheme+
  ggnewscale::new_scale_fill()+
  # Add color annotations next to the x-axis
  geom_tile(data = x_annotations, aes(x = Ref, y = -0.5, fill = Superpopulation), inherit.aes = FALSE, height = 1) +
  
  # Add color annotations next to the y-axis
  geom_tile(data = y_annotations, aes(x = -0.5, y = Contrast, fill = Superpopulation), inherit.aes = FALSE, width = 1) +
  
  scale_fill_manual(values = superpopcols)+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(paste0("10_figures/01_plots/supp/26_mage/heatmap_DGE_MAGE.pdf"), dpi=700, width = 7, height = 7,  units = "in")
ggsave(paste0("10_figures/02_panels/png/supp/26_mage.png"), dpi=700, width = 7, height = 7,  units = "in")
ggsave(paste0("10_figures/02_panels/svg/supp/26_mage.svg"), dpi=700, width = 7, height = 7,  units = "in")

fwrite(meltednew, paste0("07_differential_expressions/data/02_DEGres_MAGE_sigGenesMatrix.tsv"), sep = "\t", quote = F, row.names = F)

