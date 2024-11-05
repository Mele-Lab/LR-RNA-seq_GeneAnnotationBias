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
## 0----------------------------END OF HEADER----------------------------------0
library(ggpmisc)
if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
  annot <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
}else if(TYPE=="pantrx"){
  nametype <- "PODER annotation"
  annot <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene_with_biotypes.tsv")
}

namevec <- c()
all_samples <- list()
for(file in list.files(paste0("data/", TYPE,"/04_calc_ase"), pattern="ase_annotated.tsv", full.names=T)){
  tempdata <- fread(file)
  all_samples <- append(all_samples, list(tempdata[,geneid.v:=tstrsplit(GENE_ID,",")[[1]]]))
  namevec <- c(namevec, gsub(".*_","",gsub("_ase_annotated.tsv","",file)))
}
names(all_samples) <- namevec
data <- rbindlist(all_samples, idcol="sample")
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata<- metadata[mixed_samples==F]
metadata<- metadata[merged_run_mode==T]
data <- metadata[, .(cell_line_id, sample, population, map_reads_assemblymap)][data, on=c("cell_line_id"="sample")]


# compute fdr per sample
data[, filter:=ifelse(GENOTYPE_WARNING==0 & BLACKLIST==0 & MULTI_MAPPING==0 & OTHER_ALLELE_WARNING==0 & HIGH_INDEL_WARNING==0 & rawDepth>=20, "pass", "fail")]
data[filter=="pass", fdr:=p.adjust(BINOM_P, method="BH"), by="sample"]
data[, significant:=ifelse(fdr<0.05, "FDR<0.05", "ns")]
hist(data[sample=="LWK4" & filter=="pass"]$BINOM_P  )

# sigificant-tested plot
ggplot(unique(data[filter=="pass", .(significant_genes, tested_genes, population, map_reads_assemblymap, annot, sample)]), aes(x=tested_genes, y=significant_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASTU Significant Genes", x="# ASTU Tested Genes", size="Reads (M)", col="Population")+
  scale_color_manual(values=popcols)+
  facet_wrap(~annot)


# plot number of significant variants and variants passing threshold
ggplot(data, aes(x=reorder(sample, map_reads_assemblymap), fill=filter))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c("darkred","darkgreen"))+
  labs(x="", y="# SNPs")+
  facet_wrap(~significant, scale="free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# compute number of unique genes per sample
data[, tested_genes := uniqueN(geneid.v[filter == "pass"]), by="sample"]
data[, significant_genes := uniqueN(geneid.v[filter == "pass" & significant!="ns"]), by="sample"]

# save a list with all the sig variants
datapass <- data[filter=="pass"][, .(contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, GENOTYPE, geneid.v, BINOM_P, fdr, tested_genes, significant_genes, cell_line_id, sample, population, map_reads_assemblymap)]
setnames(datapass, old = c("variantID", "BINOM_P", "fdr"), new = c("variant", "p.value", "FDR"))

fwrite(datapass, paste0("data/ASE_hits_", TYPE, ".tsv"), sep="\t", quote = F, row.names = F)
datapass <- fread(paste0("data/ASE_hits_", TYPE, ".tsv"))

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
ggplot(unique(data[, .(significant_genes, tested_genes, population, map_reads_assemblymap)]), aes(x=tested_genes, y=significant_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASE Significant Genes", x="# ASE Tested Genes", size="Reads (M)", title="PODER annotation")+
  scale_color_manual(values=popcols)


# add annotation data
datapassannot <- unique(annot[, .(geneid.v, gene_biotype)])[datapass, on="geneid.v"]
datapassannotsig <- datapassannot[FDR<0.05]

datapassannotsig[, countperbiotype := uniqueN(geneid.v), by=c("gene_biotype", "sample")]

ggplot(unique(datapassannotsig[gene_biotype%in%c("protein_coding", "lncRNA"), .(sample, gene_biotype, countperbiotype)]), aes(x=gene_biotype, col=gene_biotype, y=countperbiotype))+
  geom_jitter()+
  mytheme+
  scale_color_manual(values=c("#F79D5C","#297373" ))+
  labs(y="ASE genes", x="", title=nametype)+
  guides(color="none")






# # compute enrichment of significatn genes in total
# passdata <-unique(data[filter=="pass", adj.P.Val:=fdr][, .(adj.P.Val, geneid.v)])
# passdata[, geneid:=gsub("\\..*", "", geneid.v)]
# ora_res <- run_ora(unique(passdata[, .(geneid, adj.P.Val)]), db=dbs, keyType="ENSEMBL")
# res <-extract_ora_res_individual(ora_res)
# ggplot(res, aes(x=Count, y=reorder(Description, Count), col=qvalue))+
#   geom_point()+
#   mytheme

library(ggpmisc)
library(ggpubr)
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadataraw <- metadataraw[mixed_samples==FALSE]
metadata <- metadataraw[merged_run_mode==TRUE]
metadata[, samplecode:=paste(lab_number_sample, lab_sampleid, cell_line_id, sep="_")]

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
# load data
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
arrayenh <- fread("array_enhanced_gencode", header = F)
array <- rbind.data.frame(arraygen, arraypan)
array <- rbind.data.frame(array, arrayenh)

asts <- list()

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  temp <- fread(paste0("data/", TYPE,"/04_calc_ase/", SAMPLE,"_ase_annotated_filtered.tsv"))
  temp[, `:=`(annot=TYPE, samplecode=SAMPLE)]
  ase <- append(asts, list(temp))}

aseraw <- rbindlist(ase, use.names=TRUE)
aseraw <- metadata[, .(samplecode, sample, population, map_reads_assemblymap)][aseraw, on="samplecode"]
aseraw[, annot := ifelse(annot=="gencode", "GENCODEv47", 
                          ifelse(annot=="enhanced_gencode", "Enhanced\nGENCODEv47", "PODER"))]
aseraw[, annot:=factor(annot, levels=c("GENCODEv47", "PODER", "Enhanced\nGENCODEv47"))]

ase <- aseraw[filter=="pass"]

# computed number of tested genes
ase[, tested_genes:=uniqueN(geneid.v), by=c("sample", "annot")]
ase[FDR<0.05, significant_genes:=uniqueN(geneid.v), by=c("sample", "annot")]
ase[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "OOA")]
ggplot(unique(ase[, .(significant_genes, tested_genes, population, map_reads_assemblymap, annot, sample)]), aes(x=tested_genes, y=significant_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASTU Significant Genes", x="# ASTU Tested Genes", size="Reads (M)", col="Population")+
  scale_color_manual(values=popcols)+
  facet_wrap(~annot)
