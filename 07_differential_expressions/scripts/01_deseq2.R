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
TYPE="gencode"

# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,"_genelvl.tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODE"
  countsprefilt <- fread("../../novelannotations/quantifications/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER "
  countsprefilt <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
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
metadataunfiltered <- metadataraw[order(rownames(metadataraw)), c(covariates)]

# Sum transcript counts per gene
countsprefilt <- annot[, .(transcriptid.v, geneid.v)][countsprefilt, on=c("transcriptid.v"= "transcript_id")]
countsprefilt <- countsprefilt[, lapply(.SD, sum), by = geneid.v, .SDcols = patterns("_")]
colnames(countsprefilt) <- gsub("_1$", "", colnames(countsprefilt))
setnames(countsprefilt, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
countsprefilt <- column_to_rownames(countsprefilt, var="geneid.v")
colnames(countsprefilt)  <- gsub(".*_", "", colnames(countsprefilt))
countsprefilt <- countsprefilt[, order(colnames(countsprefilt))]


# filtering

isexpr <- rowSums(cpm(countsprefilt) > 1) >= 5
dge <- DGEList(countsprefilt[isexpr, ])


library(DESeq2)

resultslist <- list()
for(population in unique(metadataraw$population)){
  metadataunfiltered$population <- factor(metadataunfiltered$population, levels=c(population, levels(metadataraw$population)[!grepl(population,levels(metadataraw$population))]))


  dds <- DESeqDataSetFromMatrix(countData = round(dge$counts, digits=0),
                                colData = metadataunfiltered,
                                design= ~population+sex+pc1+pc2)
  dds <- DESeq(dds)
  contrastsOI <- resultsNames(dds)[grep("population", resultsNames(dds))]
  contrastvecname <- c()
  reslist <- list()
  for(contrast in contrastsOI){
    res <- results(dds, name=contrast)
    res <- as.data.table((as.data.frame(res)))
    res <- res[padj<0.05]
    reslist <- append(reslist, list(res))
    contrastvecname <- c(contrastvecname, contrast)
  }
  names(reslist) <- contrastvecname
  results <- rbindlist(reslist, idcol="contrast")
  resultslist <- append(resultslist, list(results))
}

data <- rbindlist(resultslist)

# Classify genes as upregulated or downregulated based on log2FoldChange
data$regulation <- ifelse(data$log2FoldChange > 0, "Upregulated", "Downregulated")

# Count upregulated and downregulated genes for each contrast
result <- aggregate(regulation ~ contrast, data, function(x) {
  table(factor(x, levels = c("Upregulated", "Downregulated")))
})
res <- as.data.table(result)

res[, `:=`(Contrast=gsub("vs...", "", contrast), Reference=gsub("...vs", "", contrast))]


# count number of degs per contrast
data[, sigenes := .N, by="contrast"]
data[, contrast:=gsub("_","",gsub("population_", "", contrast))]
subdata <-unique(data[, .(contrast, sigenes)])
subdata[, comparison:=tstrsplit(contrast, "vs")[[1]]]
subdata[, reference:=tstrsplit(contrast, "vs")[[2]]]

ggplot(subdata, aes(x=reference, y=comparison, fill=sigenes))+
  geom_tile()+
  mytheme+
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#A0000D", midpoint = 0, na.value = "#9AC7CB")+
  geom_text(data=subdata, aes(label=sigenes, size=sigenes),  na.rm =T)+
  guides(size="none")
ggplot(res, aes(x = Reference, y = Contrast, fill = regulation.Upregulated,)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0080AF", mid = "white", high = "#A0000D", midpoint = 0,limits=c(0,300), na.value = "#9AC7CB") +
  theme_minimal() +
  labs(title=paste0(nametype, " using DESeq2"),
       x = "Against Reference",
       y="Upregulated Genes in",
       fill = "# DEGs") +
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title=element_text(face="bold", hjust=0.5)) +
  geom_text(data=res[!is.na(regulation.Upregulated)], aes(label=abs(regulation.Upregulated), size=abs(regulation.Upregulated)),  na.rm =T)+
  coord_fixed()+
  scale_size_continuous(range=c(4.5,8)*0.35)+
  guides(size="none")+
  scale_y_discrete(limits = levels(res$Contrast))
ggsave(paste0("../10_figures/01_plots/supp/25_deseqdropout/heatmap_DEG_",nametype,"_DESeq2.pdf"), dpi=700, width = 3, height = 3,  units = "in")
