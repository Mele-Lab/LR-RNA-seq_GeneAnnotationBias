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
gtf <- fread("../novelannotations/merged/poder_v1.gtf")
gtf$V2 <- "Cerberus"
data[, filter:=ifelse((structural_category=="NIC" & sample_sharing>1 & flair_max_counts>=3 & length >=300 & sj_less_recountsupported_counts>=50|
                         structural_category%in%c("NNC", "Genic", "Fusion", "Antisense", "Intergenic") & sample_sharing>1 & flair_max_counts>=3 & sj_less_recountsupported_counts>=50 & tool_sharing>1 & length >=300) &
                        associated_gene_biotype%in%unique(data$associated_gene_biotype)[grepl("lncRNA|^$|protein_coding", unique(data$associated_gene_biotype))] &
                        contig!="chrM", "pass", "fail")]
# keep only filtered isoforms
filtereddata <- data[filter=="pass" & structural_category!="FSM"]
fwrite(filtereddata, "04_transcriptome_assembly/04_evaluation/05_mastertable/data/subset_PODER/241203_PODER_noveltrx_recountsupp50_trx_mastertable.tsv", quote = F, row.names = F, col.names = T, sep="\t")

# Step 1: Extract transcripts
gtf[, id:=tstrsplit(V9, "\"")[[4]]]
transcripts <- unique(filtereddata$isoform)

# Step 2: Filter gtf
filtered_gtf<- gtf[V3!="gene" & id%in%transcripts]



filtereddata <-fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/subset_PODER/241203_PODER_noveltrx_recountsupp50_trx_mastertable.tsv")

# remove some columns that are not interesting

filtereddata<- filtereddata[, .SD, .SDcols=!grepl("trx_per_asstrx_count", colnames(filtereddata))]
filtereddata<- filtereddata[, .SD, .SDcols=!grepl("trx_pergene_count", colnames(filtereddata))]
filtereddata<- filtereddata[, .SD, .SDcols=!grepl("old", colnames(filtereddata))]
filtereddata <- filtereddata[, filter:=NULL]
filtereddata<- filtereddata[, .SD, .SDcols=!grepl("existsFSMin", colnames(filtereddata))]
filtereddata<- filtereddata[, discovered_transcripts_per_gene:=NULL]
filtereddata<- filtereddata[, associated_transcript_biotype:=NULL]


# # Change the old name of some FSM genes (that had an ambiguous geneid) to the new one
# newgeneid <- filtereddata[, .(old_associated_geneid.v, associated_geneid.v, isoform)][filtered_gtf, on=c("isoform"="id")]
# newgeneid[, V9:= ifelse(old_associated_geneid.v!=associated_geneid.v,
#                         gsub(old_associated_geneid.v, associated_geneid.v, V9), V9), by = 1:nrow(newgeneid)]
# newgeneid[, `:=`(geneid.v=NULL, old_associated_geneid.v=NULL, associated_geneid.v=NULL, isoform=NULL)]
# # Identify ISM without FSM
# ism <- unique(data[filter=="pass" & structural_category=="ISM", .(associated_geneid.v, associated_transcriptid.v)])
# 
# # load gencode annotation
# gencode <- fread("../../../Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf")
# gencode <- gencode[V3%in%c("transcript", "exon")]
# gencode[, id:=tstrsplit(V9, "\"")[[4]]]
# # keep only those features related to our transcripts of interest
# gencode <- gencode[id %in%ism$associated_transcriptid.v][, id:=NULL]
# colnames(gencode)[1] <- "V1"
# 
# # concatenate
# final_annot <- rbind(gencode, newgeneid)
fwrite(filtered_gtf[, id:=NULL], "04_transcriptome_assembly/04_evaluation/05_mastertable/data/subset_PODER/241203_PODER_noveltrx_recountsupp50.gtf", quote = F, row.names = F, col.names = F, sep="\t")



## CHECK if population specific trx have lower recount support

popsp <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/241128_PODER_pop_specific_transcripts.tsv")
popsp[, popsp_log:=TRUE]
popsp[, popsp_eur:=eur]
popsp[, popsp:=population]

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

# filter data
data[, filter:=ifelse(structural_category=="FSM" | (structural_category=="NIC" & sample_sharing>1 & flair_max_counts>=3 & length >=300|
                         structural_category%in%c("NNC", "Genic", "Fusion", "Antisense", "Intergenic") & sample_sharing>1 & flair_max_counts>=3 & sj_less_recountsupported_counts>=25 & tool_sharing>1 & length >=300) &
                        associated_gene_biotype%in%unique(data$associated_gene_biotype)[grepl("lncRNA|^$|protein_coding", unique(data$associated_gene_biotype))] &
                        contig!="chrM", "pass", "fail")]
# keep only filtered isoforms
filtereddata <- data[filter=="pass"]



datapop <- unique(popsp[, .(isoform, popsp_log, popsp_eur, popsp)])[filtereddata, on="isoform"]
datapop[is.na(popsp_eur), popsp_eur:="Population\nShared"]
ggplot(datapop[structural_category%in%c("FSM", "NIC", "NNC")], aes(x=popsp_eur, y=sj_less_recountsupported_counts+1, fill=popsp_eur))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.05, outliers = F)+
  mythemen+
  scale_y_continuous(trans = "log10", limits = c(NA, 1e13)) +
  ggpubr::stat_compare_means(comparisons = list(c("European", "Non-European"), 
                                                c("Population\nShared","European"), 
                                                c("Population\nShared","Non-European")),method.args = list(alternative="two.sided"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1))+
  labs(x="", y="Recount3 Support Counts\nof Least Supported SJ per Transcript")+
  scale_fill_manual(values=c("#466995", "#A53860", "darkgrey"))+
  guides(fill="none")+
  facet_wrap(~structural_category)+
  scale_x_discrete(limits=c("Population\nShared","European", "Non-European"))



ggplot(datapop[structural_category%in%c("NIC")], aes(x=popsp_eur, y=sj_less_recountsupported_counts+1, fill=popsp_eur))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.05, outliers = F)+
  mythemen+
  scale_y_continuous(trans = "log10", limits = c(NA, 1e13)) +
  ggpubr::stat_compare_means(comparisons = list(c("European", "Non-European"), 
                                                c("Population\nShared","European"), 
                                                c("Population\nShared","Non-European")),method.args = list(alternative="two.sided"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1))+
  labs(x="", y="Recount3 Support Counts\nof Least Supported SJ per Transcript")+
  scale_fill_manual(values=c("#466995", "#A53860", "darkgrey"))+
  guides(fill="none")+
  facet_wrap(~subcategory)+
  scale_x_discrete(limits=c("Population\nShared","European", "Non-European"))



ggplot(datapop, aes(x=popsp_eur, y=sj_less_recountsupported_counts+1, fill=popsp_eur))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.05, outliers = F)+
  mythemen+
  scale_y_continuous(trans = "log10", limits = c(NA, 1e13)) +
  ggpubr::stat_compare_means(comparisons = list(c("European", "Non-European"), 
                                                c("Population\nShared","European"), 
                                                c("Population\nShared","Non-European")),method.args = list(alternative="two.sided"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1))+
  labs(x="", y="Recount3 Support Counts\nof Least Supported SJ per Transcript")+
  scale_fill_manual(values=c("#466995", "#A53860", "darkgrey"))+
  guides(fill="none")+
  scale_x_discrete(limits=c("Population\nShared","European", "Non-European"))
