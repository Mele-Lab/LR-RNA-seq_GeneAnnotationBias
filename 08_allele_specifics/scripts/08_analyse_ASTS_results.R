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

TYPE<-"pantrx"

if(TYPE=="gencode"){
  nametype <- "GENCODEv47 annotation"
}else{nametype <- "PODER annotation"}

namevec <- c()
all_samples <- list()
for(file in list.files(paste0("data/", TYPE,"/05_calc_asts"), pattern="asts_annotated_variantgenefilt.tsv", full.names=T)){
  tempdata <- fread(file)
  all_samples <- append(all_samples, list(tempdata[,geneid.v:=tstrsplit(geneid.v,",")[[1]]]))
  namevec <- c(namevec, gsub(".*_","",gsub("_asts_annotated_variantgenefilt.tsv","",file)))
}
names(all_samples) <- namevec
data <- rbindlist(all_samples, idcol="sample")
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata<- metadata[mixed_samples==F]
metadata<- metadata[merged_run_mode==T]
data <- metadata[, .(cell_line_id, sample, population, map_reads_assemblymap)][data, on=c("cell_line_id"="sample")]


# is fdr significant?
data[, significant:=ifelse(FDR<0.05, "FDR<0.05", "ns")]

# compute number of unique genes per sample
data[, tested_genes := uniqueN(geneid.v), by="sample"]
data[, significant_genes := uniqueN(geneid.v[significant!="ns"]), by="sample"]


# Plot
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
ggplot(unique(data[, .(significant_genes, tested_genes, population, map_reads_assemblymap)]), aes(x=tested_genes, y=significant_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASTS Significant Genes", x="# ASTS Tested Genes", size="Reads (M)", title=nametype)+
  scale_color_manual(values=popcols)






# save results
fwrite(data, paste0("data/ASTS_results_",TYPE,"_variantgenefilt.tsv"), row.names = F, quote = F, sep="\t")


# # compute enrichment of significatn genes in total
# passdata <-unique(data[filter=="pass", adj.P.Val:=fdr][, .(adj.P.Val, geneid.v)])
# passdata[, geneid:=gsub("\\..*", "", geneid.v)]
# ora_res <- run_ora(unique(passdata[, .(geneid, adj.P.Val)]), db=dbs, keyType="ENSEMBL")
# res <-extract_ora_res_individual(ora_res)
# ggplot(res, aes(x=Count, y=reorder(Description, Count), col=qvalue))+
#   geom_point()+
#   mytheme
# library(rrvgo)
# # input enrichments to reduce the number of terms
# simMatrix <- calculateSimMatrix(res$ID,
#                                 orgdb="org.Hs.eg.db",
#                                 ont="BP",
#                                 method="Rel")
# simMatrix2 <- calculateSimMatrix(res$ID,
#                                  annoDb="org.Hs.eg.db",
#                                  ont="BP",
#                                  method="Rel")
# scores <- setNames(-log10(res$qvalue), res$ID)
# reducedTerms <- reduceSimMatrix(simMatrix,
#                                 scores,
#                                 threshold=0.7,
#                                 orgdb="org.Hs.eg.db")
# setDT(reducedTerms)
# reducedTerms[, countTerms := .N, by="parentTerm"]
# reducedTerms[, percentTerms := countTerms/.N]
# 
# ggplot(reducedTerms, aes(x=percentTerms, y=reorder(parentTerm,percentTerms), size=countTerms))+
#   geom_point()



