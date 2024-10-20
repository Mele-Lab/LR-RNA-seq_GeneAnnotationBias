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
pcas <- list()
pcanames <- c()
for(file in list.files(path="data", pattern = "01_PCA*", full.names = T)){
  pcanames <- c(pcanames, gsub(".tsv", "", gsub("data/01_PCA_", "",file)))
  pcas <- append(pcas, list(fread(file)))
}
names(pcas) <- pcanames

data <- rbindlist(pcas, idcol="annotation_feature")
data <- data[, .(annotation_feature, cell_line_id, pc1, pc2)]
data$annotation_feature <- gsub("gencode_genelvl", "GENCODEv47_Genes", data$annotation_feature)
data$annotation_feature <- gsub("gencode_transcriptlvl", "GENCODEv47_Transcripts", data$annotation_feature)
data$annotation_feature <- gsub("pantrx_genelvl", "PODER_Genes", data$annotation_feature)
data$annotation_feature <- gsub("pantrx_transcriptlvl", "PODER_Transcripts", data$annotation_feature)
# Correlate
# Step 1: Pivot the data to wide format based on 'annotation_feature'
# Assuming 'annotation_feature' has unique levels
wide_dt <- data %>%
  pivot_wider(names_from = annotation_feature, 
              values_from = starts_with("pc"))

# Step 2: Calculate the correlation for each PC between the 'annotation_feature' levels

# Define the PC columns
pc_columns <- paste0("pc", 1:2)

# Create a list to store correlation results
cor_results <- list()

# Loop through each PC to calculate correlations across annotation_feature levels
for (pc in pc_columns) {
  # Extract the columns corresponding to this PC for different annotation features
  pc_data <- wide_dt %>%
    select(starts_with(pc)) %>%
    as.matrix()
  
  # Calculate pairwise correlations between annotation features for this PC
  cor_matrix <- cor(pc_data, use = "pairwise.complete.obs")
  
  # Store the result
  cor_results[[pc]] <- cor_matrix
}

# View the correlation results for each PC
cor_results
pc1 <-cor_results$pc1
pc2 <-cor_results$pc2

# Step 2: Reshape the data from wide to long format
long_cor_pc1 <- melt(as.matrix(pc2))

# Step 3: Rename the columns for better clarity
colnames(long_cor_pc1) <- c("Row", "Column", "Value")
setDT(long_cor_pc1)
long_cor_pc1[, `:=`(Row=gsub("pc._", "", Row), Column=gsub("pc._", "", Column))]
# Step 4: Plot the heatmap using ggplot2
ggplot(long_cor_pc1, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "darkred", mid = "white", 
                       midpoint = 0.99, limit = c(0.98, 1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() +
  mytheme+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(title = "PC2", x = "Annotation-Feature Pair", y = "Annotation-Feature Pair")+
  geom_text(aes(label=round(Value,digits=3)))
