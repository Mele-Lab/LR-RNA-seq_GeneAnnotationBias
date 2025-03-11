setwd("[PATH]/mapping_minimap/isoquant_analysis/mapIsoTrToRef/gffcompare/results/")
source("~/template/pau_theme.R")

tracking.list <- list.files(".", recursive = T, pattern = ".tracking")

## import
file <- tracking.list[1]

all <- NULL

for (file in tracking.list){
    
    gffcompare <-  read.table(file, header = FALSE, sep = "\t")
    colnames(gffcompare) <- c("query_transfrag_id", "query_locus_id", "reference_gene_id", "class_code", "lifted_id")
    count <- as.data.frame(table(gffcompare$class_code), stringsAsFactors = F)
    colnames(count) <- c("category","count")
    
    
    collapse <- read.csv("[PATH]/tools/gffcompare/class_conversion_gffcmp_GENCODE.csv", header = TRUE)
    data <- merge(count, collapse, by = "category")
    
    final <- data %>% group_by(collapse) %>% summarize(count = sum(count))
    final$genome <- str_split(str_split(file, "/", simplify = T)[1], "_", simplify = T)[1]
    
    all <- rbind(all, final)
}


tmp <- aggregate(all$count, list(all$collapse), median)
tmp$prct <- (tmp$x / sum(tmp$x))*100


ggplot(all) +
    aes(x = collapse, y = count, fill = collapse) +
    geom_boxplot(outliers = F) +
    scale_fill_hue(direction = 1) +
    pauTheme + 
    geom_jitter(aes(color = genome)) + 
    scale_y_log10() +  # Logarithmic y-axis
    annotation_logticks(sides ="l") +
    theme(legend.position = "none")





