## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
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
# load data
res <- fread()
res <- data.table(results_list[[1]]%>%rownames_to_column(var="geneid.v")) # temporary
res[, geneid:=gsub("\\..*", "", geneid.v)]

# load packages
ora.packages <- list("WebGestaltR",
                     "clusterProfiler",
                     "org.Hs.eg.db",
                     "ReactomePA",
                     "DOSE",
                     "vroom",
                     "tidyr",
                     "msigdbr",
                     "GSEABase",
                     "reshape2")
lapply(ora.packages, library, character.only = TRUE)


# new version
# modification of Raquel's function: keyType and ...
ora_fun <- function(gene.list, 
                     bg.list, 
                     db, 
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 0.05, 
                     keyType, ...){ if(length(gene.list)==0){return("gene.list is empty")}
  output <- list()
  if(keyType %in% c("SYMBOL", "ENSEMBL", "ENTREZID")){
    if("GO:BP" %in% db){
      print(paste0("Computing GO:BP"))
      gobp <- enrichGO(gene     = gene.list,
                       universe =    bg.list,
                       keyType = keyType, #"ENSEMBL" "SYMBOL" "ENTREZID"
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = F)
      output <- append(output, gobp)
      names(output)[length(output)] <- "GO:BP"}
    if("GO:MF" %in% db){
      print(paste0("Computing GO:MF"))
      gomf <- enrichGO(gene     = gene.list,
                       universe = bg.list,
                       keyType = keyType, #"ENSEMBL" "SYMBOL" "ENTREZID"
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = F)
      output <- append(output, gomf)
      names(output)[length(output)] <- "GO:MF"}
    if("GO:CC" %in% db){
      print(paste0("Computing GO:CC"))
      gcc <- enrichGO(gene     = gene.list,
                      universe = bg.list,
                      keyType = keyType, #"ENSEMBL" "SYMBOL" "ENTREZID"
                      OrgDb        = org.Hs.eg.db,
                      ont          = "CC",
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff,
                      readable = F)
      output <- append(output, gcc)
      names(output)[length(output)] <- "GO:CC"}
    if("Hallmark" %in% db){
      print(paste0("Loading Hallmark"))
      hallmark_geneset <- msigdbr(species = "human", category = "H")
      setDT(hallmark_geneset)
      hallmark_geneset[, `:=`(term=gs_name, gene=ensembl_gene)][, .(term, gene)]
      print(paste0("Computing Hallmark"))
      hall <- enricher(gene     = gene.list,   # #"ENSEMBL" disabled("SYMBOL" "ENTREZID")
                       universe = bg.list,
                       TERM2GENE = hallmark_geneset,
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff)
      if(length(hall)>0){output <- append(output, hall)
      names(output)[length(output)] <- "Hallmark"}}}
  if(keyType=="ENTREZID"){
    if("KEGG" %in% db){
      print(paste0("Computing KEGG"))
      kegg <- enrichKEGG(gene     = gene.list, # only supports entrez.id
                         universe = bg.list,
                         keyType = "kegg", 
                         organism = "hsa",
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         pAdjustMethod = "BH",
                         pvalueCutoff = pvalueCutoff,
                         qvalueCutoff = qvalueCutoff,
      )
      output <- append(output, kegg)
      names(output)[length(output)] <- "KEGG"}
    if("ReactomePA" %in% db){
      print(paste0("Computing ReactomePA"))
      reactome <- enrichPathway(gene     = gene.list, # only supports entrez.id
                                universe = bg.list,
                                organism = "human",
                                minGSSize    = 10,
                                maxGSSize    = 500,
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff,
                                readable = F)
      output <- append(output, reactome)
      names(output)[length(output)] <- "ReactomePA"}
    if("DO" %in% db){
      print(paste0("Computing DO"))
      DO <- enrichDO(gene          = gene.list,  # only supports entrez.id
                     universe      = bg.list,
                     ont           = "DO",
                     pAdjustMethod = "BH",
                     minGSSize     = 10,
                     maxGSSize     = 500,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     readable      = FALSE)
      output <- append(output, DO)
      names(output)[length(output)] <- "DO"}
    if("DisGeNET" %in% db){
      print(paste0("Computing DisGeNET"))
      DisGeNET <- enrichDGN(gene          = gene.list,  # only supports entrez.id
                            universe      = bg.list,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500,
                            pvalueCutoff = pvalueCutoff,
                            qvalueCutoff = qvalueCutoff,
                            readable      = F)
      output <- append(output, DisGeNET)
      names(output)[length(output)] <- "DisGeNET"}
    if("OMIM" %in% db){
      print(paste0("Computing OMIM"))
      OMIM <- enricher(gene          = gene.list, # only supports entrez.id
                       universe      = bg.list,
                       TERM2GENE     = OMIM.term2gene,
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff)
      if(is.null(OMIM)){OMIM <- "No gene sets have size between 10 and 500"}
      output <- append(output, OMIM)
      names(output)[length(output)] <- "OMIM"}
    if("human_phenotype" %in% db){
      print(paste0("Computing human_phenotype"))
      human_pheno <- enricher(gene          = gene.list, # only supports entrez.id
                              universe      = bg.list,
                              TERM2GENE     = human_pheno.term2gene,
                              minGSSize    = 10,
                              maxGSSize    = 500,
                              pAdjustMethod = "BH",
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
      output <- append(output, human_pheno)
      names(output)[length(output)] <- "human_phenotype"}}
  else if(keyType=="SYMBOL"){output <- list()
  output.names <- character()
  
  if("H" %in% db){
    print(paste0("Computing H"))
    ora.H <- enricher(gene          = gene.list, # only supports gene symbol
                      universe      = bg.list,
                      TERM2GENE     = H,
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff)
    output <- append(output, ora.H)
    names(output)[length(output)] <- "H"
  }
  
  if("cell_marker" %in% db){
    print(paste0("Computing cell_marker"))
    ora.cell_marker <- enricher(gene       = gene.list,  # only supports gene symbol
                                universe      = bg.list,
                                TERM2GENE    = cells,
                                minGSSize    = 5,
                                maxGSSize    = 500,
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff)
    output <- append(output, ora.cell_marker)
    names(output)[length(output)] <- "cell_marker"
  }
  
  if("gwas" %in% db){
    print(paste0("Computing gwas"))
    ora.gwas <- enricher(gene         = gene.list, # SYMBOL
                         universe     = bg.list, # SYMBOL
                         TERM2GENE    = gwas.term2gene,
                         minGSSize    = 5,
                         maxGSSize    = 500,
                         pAdjustMethod= "BH",
                         pvalueCutoff = pvalueCutoff,
                         qvalueCutoff = qvalueCutoff)
    output <- append(output, ora.gwas)
    names(output)[length(output)] <- "gwas"}}
  return(output)}

# create a function to apply the enrichment analysis with all genes as background
run_ora <- function(x, db,...){
  significant <- unique(x[adj.P.Val<0.05 & !is.na(adj.P.Val), geneid])
  background <- unique(x[!is.na(geneid), geneid])
  ora_obj <- ora_fun(significant, background, db, ...) #ora.fun is a function in this Functional_enrichments.RData
  return(ora_obj)}

# apply.ora.fun <- function(x, db,...){
#   x <- as.data.frame(x)
#   sig.DGE.genes <- unique(unlist(x %>% dplyr::filter(fdr<=0.05) %>% dplyr::select(entrez.id)%>% na.omit()))
#   background <- unique(x$entrez.id)[!is.na(unique(x$entrez.id))]
#   ora.obj <- ora.fun2(sig.DGE.genes, background, db, ...) #ora.fun is a function in this Functional_enrichments.RData
#   return(ora.obj)}

# execute enrichment and save Rdata
dbs <- c("GO:BP")


# Sys.time()
ora_res <- run_ora(res, db=dbs, keyType="ENSEMBL")
#ora.res <- lapply(res, run_ora, db=dbs, keyType="ENSEMBL")
# Sys.time()


# create a function to extract the significant terms of all db in a single df
extract_ora_res_individual <- function(x){if(is.list(x)){
  if(sum(sapply(x, \(y) sum(y@result$p.adjust<=0.05)))>0){ # if there is a significant path
    print("There are significant terms")
    extraction <- lapply(x[lapply(x, \(y) sum(y@result$p.adjust <=0.05)) >0],
                         \(z) z@result %>% filter(p.adjust <=0.05)) # subset those db and significant terms
    extraction <- lapply(names(extraction), \(a){transform(extraction[[a]], pathway = a)}) # add a column with db name
    extraction <- do.call(rbind.data.frame, extraction)
    return(extraction)
  }else{print("No enriched terms")
    return("No enriched terms")}}else{return("NoSigGenes")}}

# # create a generalized version of the previous function so it can be applied to list of tissues
# extract.ora.res.general <- function(x, nested=F){
#   if(nested&is.list(x)){
#     without.nulls <- lapply(x, \(y){sapply(y, \(z) z)[!sapply(y, is.character)]}) # drop db without results
#     res.not.parsed <- lapply(without.nulls, extract.ora.res.individual) # extract db results by tissue
#     res.not.parsed.sub <- res.not.parsed[sapply(res.not.parsed, is.data.frame)] # drop tissues without enrichments
#     res.with.tissue.name <- lapply(names(res.not.parsed.sub), \(a){transform(res.not.parsed.sub[[a]], tissue = a)}) # add a column with tissue name
#     output <- do.call(rbind.data.frame, res.with.tissue.name)
#     return(output)
#   }else if(!nested){
#     return(extract.ora.res.individual(x))
#   }else{return("notlist")}
# }

extract_ora_res_individual(ora_res)

# plot
ora.extracted %>%
  filter(p.adjust <= 0.05) %>%
  mutate("odds.ratio" = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio)) %>%
  group_by(Description) %>%
  mutate(count = n()) %>%
  group_by(tissue)%>%
  mutate(count2=n())%>%
  #filter(pathway%in%c("GO:BP"))%>%
  #filter(tissue=="WholeBlood")%>%
  filter(grepl("[Rr]ib", Description) | grepl("[Tt]ranslation", Description))%>%
  ggplot(aes(y = reorder(tissue, -count2), x = reorder(Description, -count))) +
  geom_point(aes(size = odds.ratio, col = as.numeric(p.adjust))) +
  theme(axis.text.x = element_text(angle = 45,  hjust = 1, size=5))+
  xlab("")+
  ylab("")
ggsave("figures/dotplot.ORA_DGE_results.png", dpi = 96, height= 8, width= 14)


##### Now rerun ORA but only with genes that are DGE across at least 5 tissues
sig.genes <- c()
all.genes <- c()
for(tissue in DGEres.entrez){
  genes <- tissue[tissue$fdr<=0.05, "entrez.id"]
  sig.genes <- c(sig.genes, genes)
  all.genes <- c(all.genes, tissue[,"entrez.id"])
}

shared.sig.entrez <- names(table(sig.genes)[table(sig.genes)>=5])
all.sig.entrez <- unique(sig.genes)
all.entrez <- unique(all.genes)

shared.sig.entrez%>%
  as.data.frame()%>%
  write_delim("output/list.shared_sig_genes.txt", quote="none", col_names=F)
all.entrez%>%
  as.data.frame()%>%
  write_delim("output/list.all_genes_background.txt", quote="none", col_names=F)

# execute function twice using 2 different backgrounds: all the significant DGE genes and all the tested genes
dbs <- c("GO:BP", "GO:MF", "GO:CC", "Hallmark","KEGG", "ReactomePA","DO","DisGeNET", "OMIM", "human_phenotype")

ora.res.shared.allsig <- ora.fun2(shared.sig.entrez, all.sig.entrez, dbs, keyType = "ENTREZID")
ora.res.shared.all <- ora.fun2(shared.sig.entrez, all.entrez, dbs, keyType = "ENTREZID")


# extract results
parsed.ora.res.allsig <-extract.ora.res.general(ora.res.shared.allsig, nested=F)
parsed.ora.res.all <- extract.ora.res.general(ora.res.shared.all, nested=F)

# plot results
parsed.ora.res.allsig%>%
  mutate("odds.ratio" = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio)) %>%
  group_by(Description) %>%
  ggplot(.)+
  geom_point(aes(y=reorder(Description, odds.ratio), x=odds.ratio, size=Count, color=p.adjust))+
  ylab("")
ggsave("figures/dotplot.ORA_DGE_5tissues_allsigbackground.png", dpi = 96, height= 5, width= 8)

parsed.ora.res.all%>%
  mutate("odds.ratio" = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio)) %>%
  group_by(Description) %>%
  ggplot(.)+
  geom_point(aes(y=reorder(Description, odds.ratio), x=odds.ratio, size=Count, color=p.adjust))+
  ylab("")
ggsave("figures/dotplot.ORA_DGE_5tissues_allbackground.png", dpi = 96, height= 5, width= 8)










##########################
#--- FISHER ---#
##########################
rpannot <- read.delim("../global_data/ribosomal_proteins_annotation_HGNC.txt", sep="\t")
rpannot$RP <- "yes"

DGEres <- lapply(DGEres, \(x) {x$ensembl <- gsub("\\..*", "", x$ensembl_gene_id)
x})

DGEres.annot <- lapply(DGEres, left_join, rpannot, by=c("ensembl"="GeneID"))
DGEres.annot <- lapply(DGEres.annot, \(x) replace(x, is.na(x), "no"))
DGEres.annot <- lapply(DGEres.annot, \(x){x$RP <- factor(x$RP, c("yes", "no"))
x}) 


# Run Fisher test on RPs
apply.fisher.test <- function(x){
  x <- as.data.frame(x)
  contingencytable <- as.matrix(table(x$fdr>0.05, x$RP))
  print(contingencytable)
  if(nrow(contingencytable)==1){contingencytable <- rbind(c(0,0), contingencytable)
  print(contingencytable)}
  colnames(contingencytable) <- c("RP", "nonRP")
  rownames(contingencytable) <- c("DGE", "nonDGE")
  return(fisher.test(contingencytable))
}

# compute enrichment on all tissues
fisherres <- lapply(DGEres.annot, apply.fisher.test)

# extract results and build data.frame
fisherres.ext <- lapply(fisherres, \(x) data.frame("p.value"=unlist(x$p.value),
                                                   "odds.ratio"=unlist(x$estimate)))
fisher.res.df <- do.call(rbind.data.frame, fisherres.ext)
fisher.res.df$p.adjust <- p.adjust(fisher.res.df$p.value, method="BH")
fisher.res.df <- rownames_to_column(fisher.res.df, "tissue")
fisher.res.df$nomsig <- ifelse(fisher.res.df$p.value<=0.05, "sig", "non.sig")
fisher.res.df$adjsig <- ifelse(fisher.res.df$p.adjust<=0.05, "sig", "non.sig")

fisher.res.df.long <- pivot_longer(fisher.res.df, c(2,4), values_to = "statistic", names_to = "p.value")
fisher.res.df.long$p.value <- factor(fisher.res.df.long$p.value, levels=c("p.value","p.adjust"))


ggplot(fisher.res.df.long)+
  geom_col(aes(x=odds.ratio, y=reorder(tissue, statistic)))+
  geom_point(aes(x=-0.3, y=tissue, col=nomsig), shape=15, size=4)+
  scale_color_manual(values=c("blue", "red"))+
  guides(col=guide_legend(title="P-value"))+
  ggnewscale::new_scale_color()+
  geom_point(aes(x=-0.7, y=tissue, col=adjsig), shape=15, size=4)+
  scale_color_manual(values=c("darkblue", "darkred"))+
  guides(col=guide_legend(title="FDR"))+
  ylab("")+
  scale_x_continuous(expand = c(0, 0))+
  coord_cartesian(xlim=c(-1.5,23))
ggsave("figures/barplot.fisher_DGE.png", dpi = 96, height= 8, width= 10)


