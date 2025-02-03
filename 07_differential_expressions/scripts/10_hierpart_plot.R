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

gene <- fread("07_differential_expressions/data/05_hierpart_geneexpression_minuspc1pc2_pluspopsex.tsv")[, type:="Gene Expression"]
trx <- fread("07_differential_expressions/data/05_hierpart_transcriptexpression_minuspc1pc2_pluspopsex.tsv")[, type:="Transcript Expression"]
psi <- fread("07_differential_expressions/data/05_hierpart_splicing_minuspc1pc2_pluspopsex.tsv")[, type:="Splicing Event"]
colnames(psi)[1] <- "V1"
psi[, population_abs:=population_abs*100] # PSI ends up having an absolute value from 0 to 1
data <- rbind.data.frame(gene, trx)
data <- rbind.data.frame(data, psi)

ggplot(data, aes(x=population_abs, fill=type))+
  geom_density(alpha=0.5)+
  mythemen+
  labs(y="Density", x="% Variation explained by Population", fill="Data Type")


