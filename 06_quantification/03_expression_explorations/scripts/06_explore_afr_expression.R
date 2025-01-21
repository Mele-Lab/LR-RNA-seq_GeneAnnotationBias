## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/06_quantification"
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

library(ggpubr)
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

# load data
tpm <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tpm.tsv")
data <- fread("../04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
metadata[, total_pop_throughput:=sum(map_reads_generalmap), by="population"]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,4,2,6,7,3,5)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
data[, associated_gene_biotype := factor(associated_gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene"))]


# prepare tpm to plot
tpm_long <- melt(tpm, id.vars="transcript_id", variable.name="cell_line_id", value.name = "TPM")
setDT(tpm_long)
tpm_long[, cell_line_id:=gsub("_1", "", cell_line_id)]
tpm_long <- unique(data[, .(associated_gene_biotype, geneid.v, isoform,structural_category)])[tpm_long, on=c("isoform"="transcript_id")]
tpm_long <- metadata[, .(sample, cell_line_id, population, map_reads_generalmap,total_pop_throughput)][tpm_long, on="cell_line_id"]
tpm_long[, Expressed:=fifelse(TPM>0.01, "Expressed", "Not Expressed")]
tpm_long[, count_expressed:=uniqueN(isoform), by=c("sample", "Expressed")]
tpm_long[, count_expressed_cat_biotype:=uniqueN(isoform), by=c("sample", "Expressed", "associated_gene_biotype", "structural_category")]
tpm_long[, count_expressed_cat:=uniqueN(isoform), by=c("sample", "Expressed", "structural_category")]
tpm_long[, eur:=fifelse(population%in%c("CEU", "AJI"), "European", "non-European")]
tpm_long[, afr:=fifelse(population%in%c("YRI", "MPC", "LWK"), "African", "OOA")]

# plot distribution of tpm
ggplot(tpm_long, aes(x=TPM, color=afr))+
  stat_ecdf(geom="step")+
  scale_color_manual(values=c("#F7D257", "#496F5D"))+
  mytheme+
  labs(x="# TPM", y="Cumulative Proportion", color="")+
  theme(legend.position="top",
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))+
  scale_x_continuous(trans="log10")+
  geom_vline(xintercept=0.01, linetype="dashed", color="darkgrey")+
  coord_cartesian(xlim = c(0.01, 1e5))+
  ggmagnify::geom_magnify(from = c(xmin = 0.03, xmax = 0.2, ymin = 0, ymax = 0.25), 
                          to = c(xmin =10, xmax = 100000, ymin = 0.1, ymax = 0.85))+
  mythemen
  # annotate(geom="text", label=paste0("p<2.22e-16"), y=0.3, x=0.75, size=6*0.35) # pvalue obtained from below plot
# ggsave(plot=p,"10_figures/01_plots/supp/22_afr_trx/ecdf_ExpressedTrx_PerGene_AFROOA.pdf", dpi=700, width = 3, height = 3,  units = "in")

t.test(tpm_long[afr=="African" & Expressed=="Expressed", TPM], tpm_long[afr=="OOA" & Expressed=="Expressed", TPM])
