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
metadata$batch <- factor(metadata$batch)
metadata$population <- factor(metadata$population)
metadata$population <- factor(metadata$population, levels=c("CEU",levels(metadata$population)[!grepl("CEU",levels(metadata$population))]))
metadata$sex <- factor(metadata$sex)


param <- SnowParam(4, "SOCK", progressbar = TRUE)



# We now have the counts, gender of each sample and annotation (gene symbol and chromosome) for each Ensemble gene. We can form a DGElist object using the edgeR package.
library(missMethyl)
library(edgeR)
y <- DGEList(counts=counts) # potser em demanar que li doni els gens també

# drop lowly expressed genes by keeping genes with at least 1 count per million reads in at least 20 samples. Finally we perform scaling normalisation.
isexpr <- rowSums(cpm(y)>1) >= 85
y <- y[isexpr,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# PREPARE CONTRASTS #########################################################################........................................................
combinations <- expand.grid(levels(metadata$population),
                            levels(metadata$population))
# combinations<- combn(levels(metadata$population)[!grepl("CEU", levels(metadata$population))], m=2, simplify=F)

# Generate the contrasts
contrasts <- apply(combinations,1, function(pair) {
  paste0(pair[1], "vs", pair[2], " = \"population", pair[1], " - population", pair[2], "\"")
})
mycon <-paste(contrasts, collapse=", ")

# Convert the contrasts into a named vector
form <- ~0+population+sex+batch

design <- model.matrix(form, metadata)

contrast_list <- eval(parse(text = paste0("list(", mycon, ")")))
contrast_expr <- do.call(makeContrasts, c(contrast_list, list(levels = design)))

#We set up the design matrix and test for differential variability.


fitvar.contr <- varFit(y, design=design, coef=grep("population", colnames(design)))
fitvar.contr <- contrasts.varFit(fitvar.contr, contrasts=contrast_expr)
# extract results
res <-data.frame(summary(decideTests(fitvar.contr)))
setDT(res)
res[, `:=`(comparison=tstrsplit(Var2,"vs")[[1]], reference=tstrsplit(Var2, "vs")[[2]])]
reswide <- dcast(res, Var2+comparison+reference~..., value.var="Freq")
reswide[, `:=`(Var2=NULL, NotSig=NULL)]
colnames(reswide)[3:4] <- c("down", "up")
# Extract sig genes if desires
allgenes <-lapply(colnames(contrast_expr), \(x) topVar(fitvar.contr, number=nrow(y$counts), coef=x))
allgenes <- lapply(allgenes, \(x) rownames_to_column(x, var="geneid.v"))
names(allgenes) <- colnames(contrast_expr)
allgenes <- rbindlist(allgenes, idcol = "contrast")
allgenes <- allgenes[!is.na(Adj.P.Value)]
allgenes[, `:=`(comparison=tstrsplit(contrast, "vs")[[1]], reference=tstrsplit(contrast, "vs")[[2]])]
allgenes[, hypervariableDVG := sum(DiffLevene>0 & Adj.P.Value < 0.05), by = "contrast"]
allgenes[, hypovariableDVG := sum(DiffLevene < 0 & Adj.P.Value < 0.05), by = "contrast"]



colnames(allgenes)[grep("Adj.P.Val", colnames(allgenes))] <- "fdr"
colnames(allgenes)[grep("contrast", colnames(allgenes))] <- "contrast_name"
colnames(allgenes)[grep("comparison", colnames(allgenes))] <- "contrast"
colnames(allgenes)[grep("P.Value", colnames(allgenes))] <- "pval"
# fwrite(allgenes, paste0("data/03_DVGres_", TYPE,".tsv"), quote = F, sep = "\t", row.names = F)
#allgenes <- fread( paste0("data/03_DVGres_", TYPE,".tsv"))

### PARSE RESULTS TO CREATE HEATMAP WITH PAIRWISE COMPARISONS
# Reorder columns to desired format
result <- unique(reswide[, .(comparison, reference, up, down)])
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
  up_value <- result$up[i]
  down_value <- result$down[i]
  
  # Assign values to the matrix
  result_matrix[before, after] <- up_value      # Upregulated in upper triangle
  result_matrix[after, before] <- down_value     # Downregulated in lower triangle
}

# Melt the data frame to long format
melted_df <- melt(result_matrix, value.name = "DVGs")

# Remove NA values for better visualization
melted_df <- melted_df[!is.na(melted_df$DVGs), ]
colnames(melted_df) <- c( "Contrast","Ref", "DVGs")

melted_df<-rbind.data.frame(melted_df, data.frame(cbind("Contrast"=levels(melted_df$Contrast),
                                                        "Ref"=levels(melted_df$Contrast), 
                                                        "DVGs"=rep(NA, 8))))
setDT(melted_df)
melted_df[,DVGs:=as.numeric(DVGs)]
# Create the heatmap
ggplot(melted_df, aes(x = Ref, y = Contrast, fill = DVGs+1,)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#576A34",trans="log", midpoint = 0, na.value = "#9AC7CB") +
  theme_minimal() +
  labs(title="MAGE GENCODEv38",
       x = "Against Reference",
       y="Hypervariable Genes in",
       fill = "# DVGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data=melted_df[Contrast!=Ref], aes(label=DVGs, size=DVGs))+
  coord_fixed()+
  scale_size_continuous(range=c(3,5))+
  guides(size="none")+
  mytheme
# fwrite(melted_df[!is.na(DVGs)], paste0("data/03_DVGres_", TYPE, "_sigGenesMatrix.tsv"), sep = "\t", quote = F, row.names = F)


meltcounts <- melt(rownames_to_column(counts, var="geneid.v"), id.vars="geneid.v", variable.name="sample", value.name = "count")
setDT(metadata)
setDT(meltcounts)
meltcountspop <- metadata[, .(internal_libraryID, population)][meltcounts, on=c("internal_libraryID"= "sample")]
ggplot(meltcountspop[geneid.v%in%], 
       aes(x=population, y=count))+
  geom_jitter()+
  mytheme+
  labs(x="", y="Counts")+
  scale_y_continuous(trans="log10")+
  geom_text(aes(label=internal_libraryID))
