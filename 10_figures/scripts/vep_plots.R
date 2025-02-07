## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/10_figures"
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

# Plot number and proportions of VEP coding consequences
classification <- fread("data/vep/VEP_severityClassification.tsv")

enhanced <- fread("data/vep/GencodePoder_coding.tsv")[, source:="Enhanced\nGENCODE"]
gencode <- fread("data/vep/Gencode_coding.tsv")[, source:="GENCODE"]
poder <- fread("data/vep/Poder_coding.tsv")[, source:="PODER"]

data <- rbind.data.frame(enhanced, gencode)
data <- rbind.data.frame(data, poder)

data <- classification[data, on=c("SO_term"="type")]

data <- data[, Impact:= factor(Impact, levels=c("HIGH", "MODERATE","LOW", "MODIFIER"))]
impactcol <- unique(data$Color)
names(impactcol) <- unique(data$Impact)
data[, count :=sum(freq), by=.(Impact, source)]
data[, totalcount :=sum(freq), by=.(source)]
data[, percent:=paste0(round(count/totalcount*100, 0), "%")]


ggplot(unique(data[, .(Impact, source, count,percent)]), aes(x=source, y=count/10^6, fill=Impact))+
  geom_col()+
  mytheme+
  scale_fill_manual(values=impactcol)+
  labs(x="", y="# Coding Consequences (Million)")+
  geom_text(aes(label=percent), position=position_stack(vjust=0.5), size=6*0.35)
ggsave("../10_figures/01_plots/supp/40_vep/barplot_vep.pdf", dpi=700, width = 3, height = 3,  units = "in")



getnum <- unique(data[, .(Impact, source, count,percent)])
getnum[, sum:=sum(count), by="source"]
### Plot change of consequences

data <- fread("data/vep/mostSeverConsChange.tsv")

newclass <- unique(classification$Display_term)
names(newclass) <- unique(classification$SO_term)



data <- data[, `:=`(Consequence.gencode.new=newclass[Consequence.gencode],
                    Consequence.merged.new=newclass[Consequence.merged])]
ggplot(data, 
       aes(x=reorder(Consequence.gencode.new, Consequence.gencode.orderSeverity), 
           y=reorder(Consequence.merged.new,Consequence.merged.orderSeverity), 
           fill=freq.log))+
  geom_tile()+
  mytheme+
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=-1),
        legend.position = c(0.7,0.25))+
  labs(x="Consequence in GENCODE", y="Consequence in enhanced GENCODE", fill="# Consequences\n(log10)")+
  scale_x_discrete(position = "top")+
  viridis::scale_fill_viridis()
ggsave("01_plots/supp/40_vep/heatmap_vep.pdf", dpi=700, width = 4, height = 4,  units = "in")
