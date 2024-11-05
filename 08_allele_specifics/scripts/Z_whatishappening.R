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
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------
gencode <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
gencode <-gencode[, .(geneid.v, transcriptid.v, gene_biotype)]
poder <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene_with_biotypes.tsv")

# load data of calc asts
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
array <- rbind.data.frame(arraygen, arraypan)

rawasts <- list()

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  temp <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_quant.tsv"))
  temp[, `:=`(annot=TYPE, samplecode=SAMPLE)]
  rawasts <- append(rawasts, list(temp))}

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
astsoriginal <- rbindlist(rawasts)
astsoriginal <- astsoriginal[annot==TYPE & samplecode==SAMPLE]
astsraw <- astsoriginal
astsraw[, `:=`(annot=NULL, samplecode=NULL)]
# prepare raw data
astsraw <- annot[astsraw, on=c("transcriptid.v"="transcript")]
astsraw[, variant := paste(contig, position, refAllele, altAllele, sep="_")]   
astsraw[, `:=`(contig=NULL, position=NULL, refAllele=NULL, altAllele=NULL)]

# Transcripts with at least 10 counts
astsraw[, Count:=refCount+altCount][, trx_tencounts:=Count>=10]
# Gene-variant pairs with less than 20 counts
astsraw[, gene_twentycount:=sum(Count)>=20, by=c("variant", "geneid.v")]

# Keep only genes with more than 1 transcript
astsraw[trx_tencounts==TRUE, gene_numtrx := uniqueN(transcriptid.v) , by=c("variant","geneid.v")]

# Remove genes where all the counts are on the same allele because they are most likely wrongly genotyped (or not trustable)
astsraw[trx_tencounts==TRUE & gene_numtrx>1, gene_heterozygous:=(sum(refCount)!=0 & sum(Count)!=sum(refCount)), by=c("variant","geneid.v")]


astsraw[gene_numtrx>1 & gene_heterozygous==TRUE , gene_testable:=TRUE]
astsraw[, gene_variant:=paste(geneid.v, variant, sep=":")]

testable <- unique(astsraw[gene_testable==TRUE, geneid.v])



###################### now the original code
data <- astsoriginal
data[, variant := paste(contig, position, refAllele, altAllele, sep="_")]   
data <- annot[data, on=c("transcriptid.v"="transcript")]
data[, `:=`(contig=NULL, position=NULL, refAllele=NULL, altAllele=NULL)]
data[, totalCount := refCount+altCount]

# remove transcript-variant pairs with less than 10 counts
data <- data[totalCount>=10]

# Remove those gene-variant pairs with less than 20 counts
data <- data[, .SD[sum(totalCount)>=20], by=c("variant", "geneid.v")]

# Keep only genes with more than 1 transcript
data <- data[, .SD[uniqueN(transcriptid.v)>1], by=c("variant", "geneid.v")]

# Create gene-variant pairs to test all the variants in a gene
data[, gene_variant:=paste(geneid.v, variant, sep=":")]
data <- data[, .SD[(sum(refCount)!=0 & sum(totalCount)!=sum(refCount))], by="gene_variant"]

testable_original <- unique(data[, gene_variant])
perce <- (sum(testable%in%testable_original)+sum(testable_original%in%testable))/length(c(testable, testable_original))*100

cat("Gene-variant overlap of", SAMPLE, TYPE, perce, "%")
}