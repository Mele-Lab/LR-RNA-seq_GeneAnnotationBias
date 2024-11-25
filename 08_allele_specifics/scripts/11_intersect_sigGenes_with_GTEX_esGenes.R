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
metadata <- fread("..//00_metadata/data/pantranscriptome_samples_metadata.tsv")
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
# asemer <- eGenes[][unique(ase[, .(geneid, annot, population, sample, sig)]), on=c("geneid")]
# asemer[is.na(eGene), eGene:=FALSE]
# astsmer <- sGenes[unique(asts[, .(geneid, annot, population, sample, sig)]), on=c("geneid")]
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
# 
# asemer[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR", "OOA")]
# astsmer[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR", "OOA")]


# overlap between ASE and ASTS
ggplot(asqtl[ASTS!="not Tested"], aes(x=ASTS, fill=ASE))+
  geom_bar(alpha=0.75)+
  mytheme+
  facet_wrap(~annot)+
  scale_fill_manual(values=c("#356CA1","darkred", "darkgrey"))+
  labs(x="", y="# Genes", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))

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
asqtl_sample[, sGene:=factor(sGene, levels=c("sGene", "not sGene", "not Tested"))]
asqtl_sample[, `:=`(ase_validated=ifelse(ASE=="ASE" & eGene=="eGene", "validated",
                                         ifelse(is.na(ASE), "not tested", "not validated")),
                    asts_validated=ifelse(ASTS=="ASTS" & sGene=="sGene", "validated",
                                         ifelse(is.na(ASTS), "not tested", "not validated")))]
asqtl_sample[, ooa:=ifelse(population%in%c("YRI", "LWK"), "AFR","OOA")]
asqtl_sample[, trx_per_cat:=uniqueN(geneid), by=c("sample", "annot", "ase_validated")]

ggplot(unique(asqtl_sample[!is.na(annot) & ase_validated!="not tested", .(ooa, trx_per_cat, annot, ase_validated, population)]), 
       aes(x=ooa, y=trx_per_cat, fill=ase_validated))+
  geom_violin(position = position_dodge(width = 0.2), alpha = 0.5, width = 0.2) +
  ggbeeswarm::geom_quasirandom(aes(col=population), width=0.05)+
  facet_wrap(~annot)+
  mytheme+
  scale_fill_manual(values=c("#356CA1","darkred", "darkgrey"))+
  scale_color_manual(values=popcols)
  

###### OVERLAP WITH GWAS CATALOG--------------------------------------------------------

gwassets <- fread("../../../../Data/GWAScatalog/modified/full_gwas_catalog.parsed4_enrichments.geneids.tsv")
library(clusterProfiler)
res_ase_gencode <- enricher(gene=asqtl[annot=="GENCODEv47" & ASE=="ASE", geneid],
         universe=asqtl[annot=="GENCODEv47" & ASE!="not Tested", geneid],
         TERM2GENE= gwassets)
res_ase_poder <- enricher(gene=asqtl[annot=="PODER" & ASE=="ASE", geneid],
                            universe=asqtl[annot=="PODER" & ASE!="not Tested", geneid],
                            TERM2GENE= gwassets)
res_astu_gencode <- enricher(gene=asqtl[annot=="GENCODEv47" & ASTS=="ASTS", geneid],
                            universe=asqtl[annot=="GENCODEv47" & ASTS!="not Tested", geneid],
                            TERM2GENE= gwassets,
                            qvalueCutoff=0.02)
res_astu_poder <- enricher(gene=asqtl[annot=="PODER" & ASTS=="ASTS", geneid],
                          universe=asqtl[annot=="PODER" & ASTS!="not Tested", geneid],
                          TERM2GENE= gwassets,
                          qvalueCutoff=0.02)

# PLOTS
ggplot(as.data.frame(res_astu_gencode), aes(x=FoldEnrichment, size = Count, y=reorder(Description, Count), color = p.adjust)) +
  geom_point(stat = "identity") +
  scale_color_gradient(low = "#7F675B", high = "grey", name = "Adjusted P-value",    limits = c(0, 0.02)) +
  labs(x = "Odds Ratio",
    y = "",
    size="Gene Count",
    title="GENCODEv47 ASTU") +
  theme(legend.position = "top")+
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")
ggsave("../10_figures/suppfig/dotplot.GWASgencode.pdf", dpi=700, width = 18, height = 10,  units = "cm")



res_astu_poder<-data.table(as.data.frame(res_astu_poder))
fwrite(res_astu_poder, "data/pantrx/ASTU_GWAS_PODER.tsv", sep="\t", row.names = F)
res_astu_poder[, newdescription:=fifelse(grepl("DL|gly|cho",Description), "Lipoprotein Levels/Composition (x11 traits)", as.character(Description))]
ggplot(res_astu_poder[p.adjust<0.02], aes(x=FoldEnrichment, size = Count, y=reorder(newdescription, Count), color = p.adjust)) +
  geom_point(stat = "identity") +
  scale_color_gradient(low = "#A2AD59", high = "grey", name = "Adjusted P-value",    limits = c(0, 0.02)) +
  labs(x = "Odds Ratio",
       y = "",
       size="Gene Count",
       title="PODER ASTU") +
  theme(legend.position = "top")+
  mytheme+
  geom_vline(xintercept = 1, linetype="dashed", color="darkgrey")
ggsave("../10_figures/fig_04/dotplot.GWASpoder.pdf", dpi=700, width = 20, height = 10,  units = "cm")


##########################----------------------------------------------------------
# # Per sample validation
# ggplot(unique(astsmer[, .(validated_genes,ooa, annot, total_genes, population, sample)]), aes(x=total_genes, 
#                                                                                               y=validated_genes))+
#   stat_poly_line(color="darkgrey")+
#   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5)+
#   geom_point(aes(col=population), size=2)+
#   mytheme+
#   labs(x="# ASTS genes", y="# sGene-validated ASTS genes")+
#   scale_color_manual(values=popcols)+
#   facet_wrap(~annot)
# 
# ggplot(unique(asemer[, .(validated_genes,ooa, annot, total_genes, population, sample)]), aes(x=total_genes, 
#                                                                                               y=validated_genes))+
#   stat_poly_line(color="darkgrey")+
#   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5)+
#   geom_point(aes(col=population), size=2)+
#   mytheme+
#   labs(x="# ASE genes", y="# eGene-validated ASE genes")+
#   scale_color_manual(values=popcols)+
#   facet_wrap(~annot)
# 
# ggplot(unique(asemer[, .(validation_rate,ooa, annot, population, sample)]), aes(x=ooa, y=validation_rate, fill=ooa))+
#   geom_violin(alpha=0.75)+
#   geom_boxplot(outliers=F, width=0.1)+
#   ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.2)+
#   mytheme+
#   facet_wrap(~annot)+
#   scale_fill_manual(values=c("#F7D257", "#496F5D"))+
#   guides(fill="none")+
#   labs(x="", y="% ASE genes validated by eGenes", title="ASE gene validation")+
#   ggpubr::stat_compare_means( comparisons=list(c("AFR", "OOA")),method = "wilcox.test", method.args = list(alternative="two.sided"))+
#   stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=23), vjust=0.5)
# ggplot(unique(astsmer[, .(validation_rate,ooa, annot, population, sample)]), aes(x=ooa, y=validation_rate, fill=ooa))+
#   geom_violin(alpha=0.75)+
#   geom_boxplot(outliers=F, width=0.1)+
#   ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.2)+
#   mytheme+
#   facet_wrap(~annot)+
#   scale_fill_manual(values=c("#F7D257", "#496F5D"))+
#   guides(fill="none")+
#   labs(x="", y="% ASTS genes validated by sGenes", title="ASTS gene validation")+
#   ggpubr::stat_compare_means( comparisons=list(c("AFR", "OOA")),method = "wilcox.test", method.args = list(alternative="two.sided"))+
#   stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=23), vjust=0.5)


