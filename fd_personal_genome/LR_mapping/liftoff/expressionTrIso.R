setwd("/[PATH]/mapping_minimap/isoquant_analysis/mapIsoTrToRef/expression")
genomes <- c("HG002", "HG00621", "HG01928", "NA18906", "NA19240", "HG01952")

genome <- genomes[1]

all <- NULL

for (genome in genomes){
    print(genome)
    expr_counts <- read.delim(paste0("[PATH]/mapping_minimap/RNAseqLR_mapping/",
                                     genome,".fq_",genome,".fa/isoquant/OUT/OUT.transcript_model_counts.tsv"))
    colnames(expr_counts) <- c("trId", "counts")
    expr_TPM <- read.delim(paste0("[PATH]/mapping_minimap/RNAseqLR_mapping/",
                                  genome,".fq_",genome,".fa/isoquant/OUT/OUT.transcript_model_tpm.tsv"))
    colnames(expr_TPM) <- c("trId", "TPM")
    
    expr <- merge(expr_counts, expr_TPM, by = "trId")
    
    
    sqanti <-read.delim(paste0("[PATH]/mapping_minimap/isoquant_analysis/mapIsoTrToRef/expression/",
                                     genome,"_hg38/test_classification.txt_tmp"))
    expr$structural_category <-sqanti$structural_category[match(expr$trId, sqanti$isoform )]
    expr <- expr[complete.cases(expr), ]
    expr$genome <- genome
    all <- rbind(all, expr)
    
}

tmp <- aggregate(all[,2:3], list(all$structural_category, all$genome), median)

ggplot(expr) +
    aes(x = structural_category, y = TPM, fill = structural_category) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 3))

ggplot(expr) +
    aes(x = structural_category, y = counts, fill = structural_category) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 20))    

ggplot(tmp) +
    aes(x = Group.1, y = TPM, fill = Group.1) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    geom_jitter(aes(color = Group.2)) +
    theme_minimal()

