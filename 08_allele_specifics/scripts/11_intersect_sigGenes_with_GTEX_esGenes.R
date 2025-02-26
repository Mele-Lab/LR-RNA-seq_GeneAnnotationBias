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


# load metadata
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
# load eGenes and sGenes
eGenes <- fread("../../../../Data/GTExv8/GTEx_Analysis_v8_eQTL/Cells_EBV-transformed_lymphocytes.v8.egenes.txt.gz")
eGenes <- eGenes[, eGene:=ifelse(qval<0.05, "eGene", "not eGene")][, geneid:=gsub("\\..*", "", gene_id)]
eGenes <- unique(eGenes[, .(geneid, eGene)])
sGenes <- fread("../../../../Data/GTExv8/GTEx_Analysis_v8_sQTL/Cells_EBV-transformed_lymphocytes.v8.sgenes.txt.gz")
sGenes <- sGenes[, sGene:=ifelse(qval<0.05, "sGene", "not sGene")][, geneid:=gsub("\\..*", "", gene_id)]
sGenes <- unique(sGenes[, .(geneid, sGene)])

# load ASE and ASTS
ase_gencode <- fread("data/ASE_results_gencode.tsv")[, annot:="GENCODEv47"]
ase_poder <- fread("data/ASE_results_pantrx.tsv")[, annot:="PODER"]
ase <-rbind.data.frame(ase_gencode, ase_poder)[, geneid:=gsub("\\..*", "", geneid.v)]
rm(ase_gencode, ase_poder)
asts <- fread("data/ASTS_results_bothannots.tsv")[, geneid:=gsub("\\..*", "", geneid.v)]
asts <- asts[annot%in%c("GENCODEv47", "PODER")]
asts <- asts[gene_testable==TRUE]
ase <- unique(ase[, ASE:= ifelse(any(FDR < 0.05), "ASE", "not ASE"), by=c("geneid", "annot")][, .(geneid, ASE, annot)])
asts <- unique(asts[, ASTS:= ifelse(any(FDR < 0.05), "ASTS", "not ASTS"), by=c("geneid", "annot")][, .(geneid, ASTS, annot)])

# full join everything
ase_asts <- merge(ase, asts, by = c("geneid", "annot"), all = TRUE)
ase_asts <- ase_asts[geneid!=""]
egenes_sgenes <- merge(eGenes, sGenes, by = c("geneid"), all = TRUE)
asqtl <- merge(ase_asts, egenes_sgenes, by = c("geneid"), all = TRUE)
asqtl[is.na(asqtl)] <- "not Tested"
asqtl[, sGene:=factor(sGene, levels=c("sGene", "not sGene", "not Tested"))]
# asemer <- eGenes[unique(ase[, .(geneid, annot, population, sample)]), on=c("geneid")]
# asemer[is.na(eGene), eGene:=FALSE]
# astsmer <- sGenes[unique(asts[, .(geneid, annot, population, sample)]), on=c("geneid")]
# astsmer[is.na(sGene), sGene:=FALSE]
# 
# asemer[, `:=`(
#   total_genes = uniqueN(geneid),
#   validated_genes = uniqueN(geneid[eGene == TRUE]),
#   validation_rate = (uniqueN(geneid[eGene == TRUE]) / uniqueN(geneid)) * 100
# ), by = c("annot", "sample")]
# astsmer[, `:=`(
#   total_genes = uniqueN(geneid),
#   validated_genes = uniqueN(geneid[sGene == TRUE]),
#   validation_rate = (uniqueN(geneid[sGene == TRUE]) / uniqueN(geneid)) * 100
# ), by = c( "annot", "sample")]
# asemer[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR", "OOA")]
# astsmer[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR", "OOA")]


# overlap between ASE and ASTS
ggplot(asqtl[ASTS!="not Tested"], aes(x=ASTS, fill=ASE))+
  geom_bar(alpha=0.75)+
  mytheme+
  facet_wrap(~annot)+
  scale_fill_manual(values=c("#356CA1","darkred", "darkgrey"))+
  labs(x="", y="# Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)
ggsave("../10_figures/01_plots/supp/35_as_validation/barplot_ASTS_ASE_overlap.pdf", dpi=700, width=3.25, height = 3, units="in")

# overlap betwen ASE and eGenes
ggplot(asqtl[ASE!="not Tested"], aes(x=ASE, fill=eGene))+
  geom_bar(alpha=0.75)+
  mytheme+
  facet_wrap(~annot)+
  scale_fill_manual(values=c("#356CA1","darkred", "darkgrey"))+
  labs(x="", y="# Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))

# overlap betwwwn ASTS and sGenes
ggplot(asqtl[ASTS!="not Tested"], aes(x=ASTS, fill=sGene))+
  geom_bar(alpha=0.75)+
  mytheme+
  facet_wrap(~annot)+
  scale_fill_manual(values=c("#356CA1","darkred", "darkgrey"))+
  labs(x="", y="# Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))


##### FISHERS prepare fisher function----------------------------------------------------------
run_fisher <- function(data, annot, analysis1, analysis2){
  analysis <- c(analysis1, analysis2)
  data <- data[annot==annot, .SD, .SDcols = c(analysis, "geneid")]
  con_table <-table(data[, .SD, .SDcols=analysis])
  con_table <- con_table[rownames(con_table)!="not Tested", colnames(con_table)!="not Tested"]
  test_result <- fisher.test(con_table)
  res <- data.frame(
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate,
    conf_low = test_result$conf.int[1],
    conf_high = test_result$conf.int[2]
  )
  res$annot <- annot
  res$analysis <- analysis1
  return(res)
}
# run fishers and merge
all_fishers <-list(run_fisher(asqtl, "PODER", "ASTS", "sGene"),
                   run_fisher(asqtl, "PODER", "ASE", "eGene"),
                   run_fisher(asqtl, "GENCODEv47", "ASTS", "sGene"),
                   run_fisher(asqtl, "GENCODEv47", "ASE", "eGene"))
all_fishers <- rbindlist(all_fishers)
all_fishers[, FDR:=p.adjust(p_value, method="BH")][, fdr_sig:=ifelse(FDR<0.05, "FDR<0.05", "FDR>=0.05")]
# plot fisher results

ggplot(all_fishers, aes(x=odds_ratio, y=annot, color=annot, alpha=fdr_sig))+
  geom_point(size=5) +  # Add points
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.1) + 
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  labs(x="Odds Ratio", y="", color="", alpha="")+
  scale_color_manual(values=c( "#7F675B", "#A2AD59"))+
  facet_wrap(~analysis)+
  scale_alpha_manual(values=c(1, 0.5))
ggsave("../10_figures/01_plots/supp/35_as_validation/dotplot_ASTS_ASE_enrichmentInESGenes.pdf", dpi=700, width=4, height = 3, units="in")


### Plot the validation
# load ASE and ASTS
ase_gencode <- fread("data/ASE_results_gencode.tsv")[, annot:="GENCODEv47"]
ase_poder <- fread("data/ASE_results_pantrx.tsv")[, annot:="PODER"]
ase_sample <-rbind.data.frame(ase_gencode, ase_poder)[, geneid:=gsub("\\..*", "", geneid.v)]
rm(ase_gencode, ase_poder)
asts_sample <- fread("data/ASTS_results_bothannots.tsv")[, geneid:=gsub("\\..*", "", geneid.v)]
asts_sample <- asts_sample[gene_testable==TRUE]
ase_sample <- unique(ase_sample[, ASE:= ifelse(any(FDR < 0.05), "ASE", "not ASE"), by=c("geneid", "annot", "sample")][, .(geneid, ASE, annot, sample, population)])
asts_sample <- unique(asts_sample[, ASTS:= ifelse(any(FDR < 0.05), "ASTS", "not ASTS"), by=c("geneid", "annot", "sample")][, .(geneid, ASTS, annot, sample, population)])

# full join everything
ase_asts_sample <- merge(ase_sample, asts_sample, by = c("geneid", "annot", "sample", "population"), all = TRUE)
ase_asts_sample <- ase_asts_sample[geneid!=""]
egenes_sgenes <- merge(eGenes, sGenes, by = c("geneid"), all = TRUE)
asqtl_sample <- merge(ase_asts_sample, egenes_sgenes, by = c("geneid"), all = TRUE)
asqtl_sample<- asqtl_sample[!is.na(annot)]
asqtl_sample[, sGene:=factor(sGene, levels=c("sGene", "not sGene", "not Tested"))]
asqtl_sample[, `:=`(ase_validated=ifelse(ASE=="ASE" & eGene=="eGene", "validated",
                                         ifelse(is.na(ASE), "not tested", "not validated")),
                    asts_validated=ifelse(ASTS=="ASTS" & sGene=="sGene", "validated",
                                          ifelse(is.na(ASTS), "not tested", "not validated")))]
asqtl_sample[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR","OOA")]
asqtl_sample[, trx_per_cat:=uniqueN(geneid), by=c("sample", "annot", "asts_validated")]
asqtl_sample[, percent_validation:=uniqueN(geneid[ASTS=="ASTS" & asts_validated=="validated"])/uniqueN(geneid[ASTS=="ASTS"])*100, by=c("sample", "annot")]
asqtl_sample[, percent_validation_ase:=uniqueN(geneid[ASE=="ASE" & ase_validated=="validated"])/uniqueN(geneid[ASE=="ASE"])*100, by=c("sample", "annot")]


ggplot(unique(asqtl_sample[!is.na(annot) & ase_validated!="not tested", .(ooa, annot,  population,percent_validation, sample)]), 
       aes(x=ooa, y=percent_validation, fill=ooa))+
  geom_violin(alpha = 0.5) +
  geom_boxplot(outliers = F, show.legend = F, width=0.1)+
  ggbeeswarm::geom_quasirandom(aes(col=population), alpha=0.75, size=1)+
  facet_wrap(~annot)+
  mytheme+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  scale_color_manual(values=popcols)+
  guides(fill="none")+
  labs(x="", y="% ASTU genes validated by sGenes", col="Population")+
  ggpubr::stat_compare_means( comparisons=list(c("AFR", "OOA")),method = "wilcox.test", method.args = list(alternative="two.sided"), size=6*0.35)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=33), vjust=0.5, size=6*0.35)+
  ylim(c(20,77))
ggsave("../10_figures/01_plots/supp/35_as_validation/violin_ASTS_validation.pdf", dpi=700, width = 3.2, height = 2.5,  units = "in")

ggplot(unique(asqtl_sample[!is.na(annot) & ase_validated!="not tested", .(ooa, annot,  population,percent_validation_ase, sample)]), 
       aes(x=ooa, y=percent_validation_ase, fill=ooa))+
  geom_violin(alpha = 0.5) +
  geom_boxplot(outliers = F, show.legend = F, width=0.1)+
  ggbeeswarm::geom_quasirandom(aes(col=population), alpha=0.75, size=1)+
  facet_wrap(~annot)+
  mytheme+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  scale_color_manual(values=popcols)+
  guides(fill="none")+
  labs(x="", y="% ASE genes validated by sGenes", col="Population")+
  ggpubr::stat_compare_means( comparisons=list(c("AFR", "OOA")),method = "wilcox.test", method.args = list(alternative="two.sided"), size=6*0.35)+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=20), vjust=0.5, size=6*0.35)+
  ylim(c(20,77))
ggsave("../10_figures/01_plots/supp/35_as_validation/violin_ASE_validation.pdf", dpi=700, width = 3.2, height = 2.5,  units = "in")

###### OVERLAP WITH GWAS CATALOG--------------------------------------------------------

gwassets <- fread("../../../../Data/GWAScatalog/modified/full_gwas_catalog.parsed4_enrichments.geneids.tsv")
library(clusterProfiler)



# tag hla
hla <-  fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.HLA_proteincoding_genes.tsv", header=F)
hla[, gene:=gsub("\\..*","", V1)][, V1:=NULL]
hla <- hla$gene

asqtl[, hla:=fifelse(geneid%in%hla, "HLA", "notHLA")]

res_ase_gencode <- enricher(gene=asqtl[annot=="GENCODEv47" & ASE=="ASE", geneid],
                            universe=asqtl[annot=="GENCODEv47" & ASE!="not Tested", geneid],
                            TERM2GENE= gwassets)
res_ase_poder <- enricher(gene=asqtl[annot=="PODER" & ASE=="ASE", geneid],
                          universe=asqtl[annot=="PODER" & ASE!="not Tested", geneid],
                          TERM2GENE= gwassets)
res_astu_gencode <- enricher(gene=asqtl[annot=="GENCODEv47" & ASTS=="ASTS" , geneid],
                             universe=asqtl[annot=="GENCODEv47" & ASTS!="not Tested" & hla=="notHLA", geneid],
                             TERM2GENE= gwassets,
                             qvalueCutoff=0.02)
res_astu_poder <- enricher(gene=asqtl[annot=="PODER" & ASTS=="ASTS" , geneid],
                           universe=asqtl[annot=="PODER" & ASTS!="not Tested" & hla=="notHLA", geneid],
                           TERM2GENE= gwassets,
                           qvalueCutoff=0.02)

library(org.Hs.eg.db)
resgo_astu_poder <- enrichGO(gene=asqtl[annot=="PODER" & ASTS=="ASTS" & hla=="notHLA" , geneid],
                           universe=asqtl[annot=="PODER" & ASTS!="not Tested" & hla=="notHLA", geneid],
                           keyType = "ENSEMBL",
                           OrgDb = "org.Hs.eg.db",
                           ont="BP" ,
                           qvalueCutoff=0.05)
resgo_ase_poder <- enrichGO(gene=asqtl[annot=="PODER" & ASE=="ASE" & hla=="notHLA" , geneid],
                             universe=asqtl[annot=="PODER" & ASE!="not Tested" & hla=="notHLA", geneid],
                             keyType = "ENSEMBL",
                             OrgDb = "org.Hs.eg.db",
                             ont="BP" ,
                             qvalueCutoff=0.05)
as.data.frame(resgo_astu_poder)
as.data.frame(resgo_ase_poder)


enhanced <- fread("data/ASTS_results_threeannots.tsv", sep="\t")
enhanced <- unique(enhanced[annot=="EnhancedGENCODE" & !is.na(FDR), .(FDR, geneid.v)])
enhanced[, geneid:=gsub("\\..*", "", geneid.v)]
res_astu_enhanced <- enricher(gene=enhanced[FDR<0.05, geneid],
                           universe=enhanced[, geneid],
                           TERM2GENE= gwassets,
                           qvalueCutoff=0.02)
# PLOTS
ggplot(as.data.frame(res_astu_gencode), aes(x=FoldEnrichment, size = Count, y=reorder(Description, Count), color = p.adjust)) +
  geom_point(stat = "identity") +
  scale_color_gradient2(low = "#101108", mid="#A09384", high = "#F0E4CC", name = "FDR",    limits = c(0, 0.02), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="GENCODE ASTU") +
  theme(legend.position = "top")+
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))
ggsave("../10_figures/01_plots/supp/36_gwas/dotplot.GWASgencode.pdf", dpi=700, width = 4, height = 2.5,  units = "in")



res_astu_poder<-data.table(as.data.frame(res_astu_poder))
fwrite(res_astu_poder, "data/pantrx/ASTU_GWAS_PODER.tsv", sep="\t", row.names = F)
res_astu_poder[, newdescription:=fifelse(grepl("DL|gly|cho",Description), "Lipoprotein Levels/Composition (x11 traits)", as.character(Description))]
ggplot(res_astu_poder[p.adjust<0.02], aes(x=FoldEnrichment, size = Count, y=reorder(newdescription, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#929c4d", high = "#cfd5aa", name = "FDR",    limits = c(0, 0.02), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="PODER ASTU") +
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))
ggsave("../10_figures/01_plots/main/fig_04/dotplot.GWASpoder.pdf", dpi=700, width = 4, height = 2.5,  units = "in")




res_astu_poder[, newdescription:=fifelse(grepl("DL|gly|cho",Description), "Lipoprotein Levels/Composition (x11 traits)", as.character(Description))]

dt <- res_astu_poder[newdescription=="Lipoprotein Levels/Composition (x11 traits)",]


res_astu_enhanced<-data.table(as.data.frame(res_astu_enhanced))
res_astu_enhanced[, newdescription:=fifelse(grepl("DL|gly|cho",Description), "Lipoprotein Levels/Composition (x15 traits)", as.character(Description))]
ggplot(res_astu_enhanced[p.adjust<0.02], aes(x=FoldEnrichment, size = Count, y=reorder(newdescription, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#929c4d", high = "#cfd5aa", name = "FDR",    limits = c(0, 0.02), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="Enhanced GENCODE ASTU") +
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))
ggsave("../10_figures/01_plots/supp/36_gwas/dotplot.GWASenhanced.pdf", dpi=700, width = 4, height = 2.5,  units = "in")



# Step 1: Split the `geneID` column into individual genes
gene_counts <- dt[, .(gene = unlist(strsplit(geneID, "/"))), by = .(ID)]

# Step 2: Count occurrences of each gene
gene_summary <- gene_counts[, .N, by = gene]
setnames(gene_summary, "N", "count")  # Rename the column for clarity

# Step 3: Sort by count (optional: filter for top genes)
gene_summary <- gene_summary[order(-count)]

# Step 4: Create the bar plot
ggplot(gene_summary, aes(x = reorder(gene, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene",
       y = "# Lipoprotein-related GWAS Traits\nincluding a Gene") +
  mytheme+  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
ggsave("../10_figures/01_plots/supp/36_gwas/barplot.GWAStrait_redundancy.pdf", dpi=700, width = 4, height = 2.5,  units = "in")



##########################----------------------------------------------------------
### Are my ASTU genes enriched in HLA genes?
hla <-  fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.HLA_proteincoding_genes.tsv", header=F)
hla[, gene:=gsub("\\..*","", V1)][, V1:=NULL]
hla <- hla$gene
bckg <- asqtl[annot=="GENCODEv47" & ASTS!="not Tested", geneid]
sig <- asqtl[annot=="GENCODEv47" & ASTS=="ASTS", geneid]


hla <- intersect(hla, bckg)
sig <- intersect(sig, bckg)

# Contingency table
in_hla_and_sig <- length(intersect(hla, sig))
in_hla_not_sig <- length(setdiff(hla, sig))
not_hla_and_sig <- length(setdiff(sig, hla))
not_hla_not_sig <- length(setdiff(bckg, union(hla, sig)))

contingency_table <- matrix(
  c(in_hla_and_sig, in_hla_not_sig, not_hla_and_sig, not_hla_not_sig),
  nrow = 2,
  dimnames = list(
    "Sig" = c("In Sig", "Not in Sig"),
    "HLA" = c("In HLA", "Not in HLA")
  )
)

# Perform Fisher's Exact Test
fisher_test <- fisher.test(contingency_table)

#### Are my ASE or ASTU enriched in any BP?
library(clusterProfiler)
library(org.Hs.eg.db)
raw_gen_ase <- enrichGO(keyType = "ENSEMBL",
                        gene=asqtl_sample[ASE=="ASE" & annot=="GENCODEv47", geneid],
                        universe = asqtl_sample[!is.na(ASE) & annot=="GENCODEv47", geneid],
                        ont="BP",
                        OrgDb = "org.Hs.eg.db",
                        qvalueCutoff = 0.05)
raw_gen_astu <- enrichGO(keyType = "ENSEMBL",
                        gene=asqtl_sample[ASTS=="ASTS" & annot=="GENCODEv47", geneid],
                        universe = asqtl_sample[!is.na(ASTS) & annot=="GENCODEv47", geneid],
                        ont="BP",
                        OrgDb = "org.Hs.eg.db",
                        qvalueCutoff = 0.05)
raw_pod_ase <- enrichGO(keyType = "ENSEMBL",
                        gene=asqtl_sample[ASE=="ASE" & annot=="PODER", geneid],
                        universe = asqtl_sample[!is.na(ASE) & annot=="PODER", geneid],
                        ont="BP",
                        OrgDb = "org.Hs.eg.db",
                        qvalueCutoff = 0.05)
raw_pod_astu <- enrichGO(keyType = "ENSEMBL",
                         gene=asqtl_sample[ASTS=="ASTS" & annot=="PODER", geneid],
                         universe = asqtl_sample[!is.na(ASTS) & annot=="PODER", geneid],
                         ont="BP",
                         OrgDb = "org.Hs.eg.db",
                         qvalueCutoff = 0.05)

# ggplot(as.data.frame(raw_gen_astu), aes(x=FoldEnrichment, size = Count, y=reorder(Description, Count), color = p.adjust)) +
#   geom_point(stat = "identity", alpha=0.85) +
#   scale_color_gradient2(low = "#101108", mid="#A09384", high = "#F0E4CC", name = "FDR",    limits = c(0, 0.05), midpoint=0.01) +
#   labs(x = "Odds Ratio",
#        y = "",
#        size="Gene Count",
#        title="GENCODE ASTU") +
#   mythemen+
#   geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
#   theme(legend.key.size = unit(0.2, "cm"),
#         legend.margin = margin(0, 0, 0, 0),
#         legend.box.margin = margin(-10, 3, -10, -7),
#         plot.title = element_text(family="Helvetica", face="bold"))+
#   scale_size_continuous(range = c(0.5, 4))

##### REMOVE REDUNDANCY
library(GOSemSim)
library(AnnotationHub)
library(rrvgo)
# GENCODE ASE
simMatrix <- calculateSimMatrix(raw_gen_ase$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

res_gen_ase <- as.data.table(reducedTerms)[as.data.table(raw_gen_ase), on=c("term"="Description")]
ggplot(res_gen_ase, aes(x=FoldEnrichment, size = Count, y=reorder(parentTerm, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#A09384", high = "#F0E4CC", name = "FDR",    limits = c(0, 0.05), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="GENCODE ASE") +
  mythemen+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))

# PODER ASE
simMatrix <- calculateSimMatrix(raw_pod_ase$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

res_pod_ase <- as.data.table(reducedTerms)[as.data.table(raw_pod_ase), on=c("term"="Description")]
ggplot(res_pod_ase, aes(x=FoldEnrichment, size = Count, y=reorder(parentTerm, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#929c4d", high = "#cfd5aa", name = "FDR",    limits = c(0, 0.05), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="PODER ASE") +
  mythemen+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))

# GENCODE ASTU
simMatrix <- calculateSimMatrix(raw_gen_astu$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

res_gen_astu <- as.data.table(reducedTerms)[as.data.table(raw_gen_astu), on=c("term"="Description")]
ggplot(res_gen_astu, aes(x=FoldEnrichment, size = Count, y=reorder(parentTerm, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#A09384", high = "#F0E4CC", name = "FDR",    limits = c(0, 0.05), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="GENCODE ASTU") +
  mythemen+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))

# PODER ASTU
simMatrix <- calculateSimMatrix(raw_pod_astu$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

res_pod_astu <- as.data.table(reducedTerms)[as.data.table(raw_pod_astu), on=c("term"="Description")]
ggplot(res_pod_astu, aes(x=FoldEnrichment, size = Count, y=reorder(parentTerm, Count), color = p.adjust)) +
  geom_point(stat = "identity", alpha=0.85) +
  scale_color_gradient2(low = "#101108", mid="#929c4d", high = "#cfd5aa", name = "FDR",    limits = c(0, 0.05), midpoint=0.01) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="PODER ASTU") +
  mythemen+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        plot.title = element_text(family="Helvetica", face="bold"))+
  scale_size_continuous(range = c(0.5, 4))

