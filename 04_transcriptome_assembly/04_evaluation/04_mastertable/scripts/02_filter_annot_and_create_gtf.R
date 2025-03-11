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
data <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo.tsv")

## FILTER ANNOTATION
gtf <- fread("../novelannotations/merged/240917_merge_geneEntry_correctedScaffolds_nochrEBV.sorted.gtf")
gtf$V2 <- "Cerberus"
data[, filter:=ifelse((structural_category=="FSM" |
                          structural_category=="ISM" & existsFSMinTranscript==F |
                          structural_category=="NIC" & sample_sharing>1 & flair_max_counts>=3 & length >=300|
                          structural_category%in%c("NNC", "Genic", "Fusion", "Antisense", "Intergenic") & sample_sharing>1 & flair_max_counts>=3 & sj_less_recountsupported_counts>=25 & tool_sharing>1 & length >=300) &
                         associated_gene_biotype%in%unique(data$associated_gene_biotype)[grepl("lncRNA|^$|protein_coding", unique(data$associated_gene_biotype))] &
                        contig!="chrM", "pass", "fail")]
fwrite(data, "04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv", sep="\t", row.names = F, quote=F)
# keep only filtered isoforms
filtereddata <- data[filter=="pass" & structural_category!="ISM"]

# Step 1: Extract transcripts
gtf[, id:=tstrsplit(V9, "\"")[[2]]]
transcripts <- unique(filtereddata$isoform)

# Step 2: Filter gtf
filtered_gtf<- gtf[V3!="gene" & id%in%transcripts]

# Change the old name of some FSM genes (that had an ambiguous geneid) to the new one
newgeneid <- filtereddata[, .(old_associated_geneid.v, associated_geneid.v, isoform)][filtered_gtf, on=c("isoform"="id")]
newgeneid[, V9:= ifelse(old_associated_geneid.v!=associated_geneid.v,
                        gsub(old_associated_geneid.v, associated_geneid.v, V9), V9), by = 1:nrow(newgeneid)]
newgeneid[, `:=`(geneid.v=NULL, old_associated_geneid.v=NULL, associated_geneid.v=NULL, isoform=NULL)]
# Identify ISM without FSM
ism <- unique(data[filter=="pass" & structural_category=="ISM", .(associated_geneid.v, associated_transcriptid.v)])

# load gencode annotation
gencode <- fread("../../../Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf")
gencode <- gencode[V3%in%c("transcript", "exon")]
gencode[, id:=tstrsplit(V9, "\"")[[4]]]
# keep only those features related to our transcripts of interest
gencode <- gencode[id %in%ism$associated_transcriptid.v][, id:=NULL]
colnames(gencode)[1] <- "V1"

# concatenate
final_annot <- rbind(gencode, newgeneid)
fwrite(final_annot, "../novelannotations/merged/240925_filtered.gtf", quote = F, row.names = F, col.names = F, sep="\t")


