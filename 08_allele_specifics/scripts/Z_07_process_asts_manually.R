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


arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
array <- rbind.data.frame(arraygen, arraypan)

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  data <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_quant.tsv"))
  if(TYPE=="gencode"){
    genes <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcripidv_geneidv_match.tsv")
  }else if(TYPE=="pantrx"){
      genes <- fread("../04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_transcript_associatedgene_correspondence.tsv", header = F)
      colnames(genes) <- c("geneid.v", "transcriptid.v")
      }
      
  data[, variant := paste(contig, position, refAllele, altAllele, sep="_")]   
  data <- genes[data, on=c("transcriptid.v"="transcript")]
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
  data[, transcript_variant := paste(transcriptid.v, variant, sep=":")]

  
  # Step 1: Aggregate the data for each geneid.v and transcriptid.v
  agg_data <- data[, .(totalRefCount = sum(refCount), totalAltCount = sum(altCount)), by = .(gene_variant, transcript_variant)]
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
  
  fwrite(finaldata, paste0("data/", TYPE, "/05_calc_asts/", SAMPLE,"_asts_annotated_variantgenefilt.tsv"),
         quote = F, row.names = F, sep="\t")}
