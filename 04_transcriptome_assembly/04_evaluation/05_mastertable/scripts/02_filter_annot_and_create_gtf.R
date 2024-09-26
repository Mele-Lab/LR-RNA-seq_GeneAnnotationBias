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


##### PLOT ANNOTATION
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C", "#B5B5B5", "#4F4F4F", "#6E5353", "darkred")
names(colsqanti) <- unique(data$structural_category)[c(1,2,3,5,8,6,7,4)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

##### PLOTS

ggplot(data, aes(x=structural_category, fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Transcripts")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=0)+
  facet_wrap(~filter)
ggplot(data, aes(x=structural_category, fill=filter))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("darkred","#356CA1" ))+
  labs(x="", y="# Transcripts")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))
ggplot(unique(data[,.(associated_geneid.v, structural_category, filter)]), aes(x=structural_category, fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Genes")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=0)+
  facet_wrap(~filter)
ggplot(unique(data[,.(associated_geneid.v, structural_category, filter)]), aes(x=structural_category, fill=filter))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("darkred","#356CA1" ))+
  labs(x="", y="# Genes")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))

data$associated_gene_biotype <- ifelse(data$associated_gene_biotype=="", "novel/ambiguous gene", data$associated_gene_biotype )
ggplot(unique(data[filter=="pass",.(associated_geneid.v, associated_gene_biotype, structural_category, filter)]), aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#F79D5C","darkgrey","#297373" ))+
  labs(x="", y="# Genes")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))
ggplot(unique(data[filter=="pass",.(isoform,associated_geneid.v, associated_gene_biotype, structural_category, filter)]), aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#F79D5C","darkgrey","#297373"))+
  labs(x="", y="# Transcripts")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))

ggplot(unique(data[filter=="pass",.(old_associated_geneid.v, old_associated_gene_biotype, structural_category, filter)]), aes(x=structural_category, fill=old_associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=rev(c("#297373","#F79D5C", "darkgrey" )))+
  labs(x="", y="# Genes")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))




