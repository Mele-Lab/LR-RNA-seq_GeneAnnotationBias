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

# load array to run all samples
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
array <- rbind.data.frame(arraygen, arraypan)

# load annotations
poder <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene_with_biotypes.tsv")
gencode <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
gencode <- gencode[, .(geneid.v, transcriptid.v, gene_biotype)]



for(i in 1:nrow(array)){
  # prepare loop
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  
  # prepare annot
  if(TYPE=="gencode"){
    annot <- gencode
  }else if(TYPE=="pantrx"){
    annot <- poder
  }
  
  # load data
  data <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_quant.tsv"))

  # arrange info for testing
  
  data <- annot[data, on=c("transcriptid.v"="transcript")]
  data[, variant := paste(contig, position, refAllele, altAllele, sep="_")]   
  data[, `:=`(contig=NULL, position=NULL, refAllele=NULL, altAllele=NULL)]
  data[, gene_variant :=paste(geneid.v, variant, sep=":")]
  data[, transcript_variant :=paste(transcriptid.v, variant, sep=":")]
  
  # Tag transcripts with at least 10 counts
  data[, Count:=refCount+altCount][, trx_tencounts:=Count>=10]  
  # Count transcripts (per gene) with min 10 counts 
  data[trx_tencounts==TRUE, gene_numtrx := uniqueN(transcriptid.v) , by="gene_variant"]
  # Tag Gene-variant pairs with >= 20 counts considering only trx with at least 10 counts
  data[trx_tencounts==TRUE, gene_twentycount:=sum(Count)>=20, by="gene_variant"]
  # Tag genes where all the counts are on the same allele only considering transcripts with min 10 counts
    #  Justification: because they are most likely wrongly genotyped (or not trustable)
  data[trx_tencounts==TRUE, gene_heterozygous:=(sum(refCount)!=0 & sum(Count)!=sum(refCount)), by="gene_variant"]
  
  # Tag gene-variants as testeable if they are heterozygous and have more than 1 trx with at least 10 counts and have more than 19 counts
  data[trx_tencounts==TRUE & gene_numtrx>1 & gene_twentycount==TRUE & gene_heterozygous==TRUE , gene_testable:=TRUE]


# Step 1: Aggregate the data for each geneid.v and transcriptid.v
  agg_data <- data[gene_testable==TRUE, .(totalRefCount = sum(refCount), totalAltCount = sum(altCount)), by = .(gene_variant, transcript_variant)]
  agg_data <- column_to_rownames(agg_data, var="transcript_variant")
  # Step 2: Perform the Chi-square test for each geneid.v to check for non-random distribution of counts
  # Prepare for contingency table
  test_results <- list()
  
  
  for (gene in unique(agg_data$gene_variant)) {
    gene_data <- agg_data[agg_data$gene_variant == gene,]
    
    # Create a contingency table: rows are refCount and altCount, columns are transcriptid.v
    contingency_table <- as.matrix(gene_data[,c("totalRefCount", "totalAltCount")])
    
    # Perform the Chi-square test 
    chisq_test <- chisq.test(contingency_table, simulate.p.value=T) # with this option on we do not get warning and the pval correlate 0.989
    test_results[[gene]] <- chisq_test
    
  }
  
  
  astsresults <- unique(rbindlist(test_results, idcol="gene_variant"))
  astsresults[, FDR:=p.adjust(p.value, method="BH")]
  astsresults <- unique(astsresults[, .(gene_variant, statistic, p.value, FDR)])
  finaldata <- astsresults[data, on="gene_variant"]
  
  fwrite(finaldata, paste0("data/", TYPE, "/05_calc_asts/", SAMPLE,"_asts_annotated_taggedresults.tsv"),
         quote = F, row.names = F, sep="\t")}
