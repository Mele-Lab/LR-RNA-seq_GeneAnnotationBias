
library(edgeR)
TYPE="pantrx"
# Choose covariates
covariates <- c("sex", "population", "pc1", "pc2")

# LOAD DATA
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
mypca <- fread(paste0("data/01_PCA_", TYPE,".tsv"))
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  counts <- fread("../../novelannotations/v47_kallisto_quant/matrix.abundance.tsv")
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
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

# prepare new names vector
samplesnames <- metadataraw$sample
names(samplesnames) <- metadataraw$quantification_id
samplesnames <- c(samplesnames, "transcript_id"="transcriptid.v", "geneid.v"="geneid.v")

# Sum transcript counts per gene
counts <- annot[, .(transcriptid.v, geneid.v)][counts, on=c("transcriptid.v"= "transcript_id")]
colnames(counts) <- gsub("_1$", "", colnames(counts))
setnames(counts, old = names(samplesnames), new = samplesnames,skip_absent=TRUE)
colnames(counts)[1:2] <- c("trId", "geneId")

keeptranscripts <-rowSums(cpm(counts[,.SD, .SDcols=grep("[0-9]$", colnames(counts))]) > 0.1)>10

trxcounts <- counts[keeptranscripts][, trxpergene := uniqueN(trId) , by="geneId"]
trxcounts <- trxcounts[trxpergene>1]

genecounts <- trxcounts[, lapply(.SD, sum), by = geneId, .SDcols = patterns("[0-9]$")]
keepgenes <-rowSums(genecounts[,.SD, .SDcols=grep("[0-9]$", colnames(genecounts))]> 10) >10




meltedgenes <-melt(genecounts, id.vars = "geneId", value.name = "counts", variable.name = "sample" )

ggplot(meltedgenes, aes(x=counts, fill=sample))+
  geom_density(alpha=0.2)+
  scale_x_continuous(trans="log10")+
  mytheme+
  guides(fill="none")+
  xlab("Gene Counts")+
  geom_vline(xintercept=2, linetype="dashed", color="darkgrey")




## keep transcripts that have at least 9 samples with more than 0.1 counts
keeptrx <- counts[rowSums(counts[,3:ncol(counts)]>0.1)>=9,]
# keep genes with at least 4 counts in 9 samples
prepared_counts <- prepare.trans.exp(te.df = as.data.frame(keeptrx),
                                     min.transcript.exp = 0.1,
                                     min.gene.exp = 2,
                                     min.prop = 0.2,
                                     min.dispersion = 0.1,
                                     verbose=T)


