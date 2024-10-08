## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/08_allele_specifics"
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
catch_args(3, "SAMPLE", "CELLINE", "TYPE") # celline not used
##
## 0----------------------------END OF HEADER----------------------------------0
data <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_quant.tsv"))
if(TYPE=="gencode"){
  genes <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcripidv_geneidv_match.tsv")
}else if(TYPE=="pantrx"){
    genes <- fread("../04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_transcript_associatedgene_correspondence.tsv", header = F)
    colnames(genes) <- c("geneid.v", "transcriptid.v")
    }
    
    
data <- genes[data, on=c("transcriptid.v"="transcript")]
data[, totalCount := refCount+altCount]

# filter and check if they have more than 1 transcript
data <- data[totalCount>=10]
data[, SNPPerGene:=.N , by="geneid.v" ]
new_data <- unique(data[totalCount>36, .(transcriptid.v, geneid.v)])[, TrxPerGene:=.N, by="geneid.v"]
new_data <- new_data[TrxPerGene>1][, filter:="pass"]

# merge with original data to filter it out
data <- unique(new_data[, .(transcriptid.v, filter)])[data, on="transcriptid.v"]


# Step 1: Aggregate the data for each geneid.v and transcriptid.v
agg_data <- data[filter=="pass"][, .(totalRefCount = sum(refCount), totalAltCount = sum(altCount)), by = .(geneid.v, transcriptid.v)]

# Step 2: Perform the Chi-square test for each geneid.v to check for non-random distribution of counts
# Prepare for contingency table
test_results <- list()


for (gene in unique(agg_data$geneid.v)) {
  gene_data <- agg_data[geneid.v == gene]
  
  # Create a contingency table: rows are refCount and altCount, columns are transcriptid.vs
  contingency_table <- as.matrix(gene_data[, .(totalRefCount, totalAltCount)])
  
  # Perform the Chi-square test (or Fisher's exact test depending on your data)
    chisq_test <- chisq.test(contingency_table)
    test_results[[gene]] <- chisq_test

}
astsresults <- unique(rbindlist(test_results, idcol="geneid.v"))
astsresults[, FDR:=p.adjust(p.value, method="BH")]

fwrite(astsresults[, .(geneid.v, statistic, p.value, observed, expected, residuals, stdres, FDR)], paste0("data/", TYPE, "/05_calc_asts/", SAMPLE,"_asts_annotated.tsv"),
       quote = F, row.names = F, sep="\t")
