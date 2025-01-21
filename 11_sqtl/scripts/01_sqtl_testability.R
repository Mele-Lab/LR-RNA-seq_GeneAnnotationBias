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
library(vroom)

final<-data.frame()
for(annotname in c("v47","enh_v47", "poder")){
  file_path <- paste0("../novelannotations/mage/", annotname, "_kallisto/testable/matrix.abundance_testable.tpm.tsv")
  # data <- vroom(file_path, col_types = cols(), n_max = 100)
  # setDT(data)
  data<-fread(file_path)
  data[, annot:=annotname]
  final <- rbind.data.frame(final, data)
}

final[, testable:=fifelse(t_passed_exp_filt==TRUE & g_passed_n_t_filt==TRUE & g_passed_exp_filt==TRUE & g_passed_exp_min_samples_filt==TRUE, TRUE, FALSE)]
final_testable <- final[testable==TRUE]
fwrite(final_testable, "11_sqtl/data/testable_genes.tsv", sep="\t", quote = F, row.names = F)
print(unique(final_testable[, testable_genes_per_annot :=uniqueN(gid), by="annot"][, .(annot, testable_genes_per_annot)]))













