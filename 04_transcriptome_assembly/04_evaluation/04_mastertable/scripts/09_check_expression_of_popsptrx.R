## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/05_mastertable"
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

popsptrx <- fread("data/241128_PODER_pop_specific_transcripts.tsv")[, popsptrx:=TRUE]
tpm <- fread("../../../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tpm.tsv")
tpm[, trx_mean_tpm := rowMeans(.SD, na.rm = TRUE), .SDcols = is.numeric]

popsptrx_expression <- popsptrx[tpm, on=c("isoform"="transcript_id")]
popsptrx_expression[is.na(popsptrx), popsptrx:=FALSE]

ggplot(popsptrx_expression, aes(x=popsptrx, y=trx_mean_tpm, fill=popsptrx))+
  geom_violin(alpha=0.5)+
  geom_boxplot(outliers = F, show.legend = F, width=0.1)+
  scale_y_continuous(trans="log10")+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -8), size=6*0.35)+
  ggpubr::stat_compare_means(comparisons = list(c("FALSE", "TRUE")), method.args = list(alternative="two.sided"), label.y=4, size=6*0.35)+
  scale_fill_manual(values=c("#c44536", "#457b9d"))+
  guides(fill="none")+
  labs(x="", y="Mean TPM across samples")+
  annotation_logticks(sides="l")+
  scale_x_discrete(labels=c("Population-Shared\nTranscripts", "Population-Specific\nTranscripts"))
ggsave("../../../10_figures/01_plots/supp/16_popsp_description/violin_PODER_TPM_popSpecificTrx_expressions.pdf", dpi=700, width = 3, height = 2.25,  units = "in")

  
    
