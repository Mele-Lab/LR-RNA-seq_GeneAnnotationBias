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
# arraygen <- fread("array_gencode", header = F)
# arraypan <- fread("array_pantrx", header = F)
# arrayenh <- fread("array_enhanced_gencode", header = F)
# array <- rbind.data.frame(arraygen, arraypan)
# array <- rbind.data.frame(array, arrayenh)

ase <- list()

for(TYPE in c("pantrx", "gencode", "enhanced_gencode")){
  temp <- fread(paste0("data/ASE_results_", TYPE,".tsv"))
  temp[, `:=`(annot=TYPE)]
  temp <- temp[, .(annot,contig,position,variant,refAllele,altAllele,refCount,altCount,totalCount,GENOTYPE,geneid.v,p.value,FDR,tested_genes,significant_genes,cell_line_id,sample,population,map_reads_assemblymap)]
  ase <- append(ase, list(temp))}

aseraw <- rbindlist(ase, use.names=TRUE)
aseraw[, annot := ifelse(annot=="gencode", "GENCODE", 
                          ifelse(annot=="enhanced_gencode", "Enhanced\nGENCODE", "PODER"))]
aseraw[, annot:=factor(annot, levels=c("GENCODE", "PODER", "Enhanced\nGENCODE"))]
fwrite(aseraw, "data/ase_results_threeannots.tsv", quote=F, row.names = F, sep="\t")
aseraw <- fread("data/ase_results_threeannots.tsv")
# ase <- aseraw[gene_testable==TRUE]

# # computed number of tested genes
# ase[, tested_genes:=uniqueN(geneid.v), by=c("sample", "annot")]
# unique_sig_count <- ase[FDR < 0.05, .(sig_genes=uniqueN(geneid.v)), by=c("sample", "annot")]
# ase <- unique_sig_count[ase, on=c("sample", "annot")]
ase <- aseraw
ase[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "OOA")]
ggplot(unique(ase[, .(significant_genes, tested_genes, population, map_reads_assemblymap, annot, sample)]), 
       aes(x=tested_genes, y=significant_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 6*0.35) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASTU Significant Genes", x="# ASTU Tested Genes", size="Reads (M)", col="Population")+
  scale_color_manual(values=popcols)+
  facet_wrap(~annot)+
  scale_size_continuous(range = c(0.5, 4))
ggsave("../10_figures/01_plots/supp/29_as_test_disc/scatter.ASEdisc_test.pdf", dpi=700, width = 6.5, height = 3,  units = "in")

ggplot(unique(ase[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=afr, y=tested_genes, fill=afr))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  scale_size_continuous(range = c(0.5, 1.5))+
  mytheme+
  labs(x="", y="# Tested ASE Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE),size=6*0.35)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  guides(fill="none")+
  facet_wrap(~annot)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("../10_figures/01_plots/supp/30_as_afr/violin_testedASE_afrOOA.pdf", dpi=700, width = 4, height = 4,  units = "in")

ggplot(unique(ase[, .(sample, annot, significant_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=afr, y=significant_genes, fill=afr))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  scale_size_continuous(range = c(0.5, 1.5))+
  mytheme+
  labs(x="", y="# Significant ASE Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("African", "OOA")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE),size=6*0.35)+
  scale_fill_manual(values=c("#F7D257", "#496F5D"))+
  guides(fill="none")+
  facet_wrap(~annot)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("../10_figures/01_plots/supp/30_as_afr/violin_significantASE_afrOOA.pdf", dpi=700, width = 4, height = 4,  units = "in")




ggplot(unique(ase[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=tested_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  scale_size_continuous(range = c(0.5, 1.5))+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Tested ASE Genes", size="Mapped Reads (M)", col="Population", fill="")+
  geom_pwc(method="t_test", label.size=7*0.35, p.adjust.by ="panel",label = "p.adj.format",step.increase = 0.16)+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("GENCODEv47"="GENCODE", "Enhanced\nGENCODEv47"="Enhanced\nGENCODE"))
ggsave("../10_figures/01_plots/supp/31_as_annots/violin_testedASE_byAnnot.pdf", dpi=700, width = 4, height = 4,  units = "in")


ggplot(unique(ase[, .(sample, annot, significant_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=significant_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.5)+
  scale_size_continuous(range = c(0.5, 1.5))+
  geom_boxplot(outliers = F, width=0.1)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Significant ASE Genes", size="Mapped\nReads (M)", col="Population", fill="")+
  geom_pwc(method="t_test", label.size=7*0.35, p.adjust.by ="panel",label = "p.adj.format",step.increase = 0.16)+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_discrete(labels=c("GENCODEv47"="GENCODE", "Enhanced\nGENCODEv47"="Enhanced\nGENCODE"))
ggsave("../10_figures/01_plots/supp/31_as_annots/violin_sigASE_byAnnot.pdf", dpi=700, width = 4, height = 4,  units = "in")







