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
### From PG to ref

all <-fread("10_figures/data/blastSeqIdentity.tsv")

# toPlot <- "T2T"
toPlot <- "hg38"


all.genome <- all[all$ref == toPlot,]

all.genome[, meanprct:=mean(prct), by=.(ref, thr)]
annotation_data <- unique(all.genome[thr%in%c(1, 90, 99, 100), .(thr , meanprct)])
annotation_data[, xend:=c(10, 80, 90, 90)]
annotation_data[, yend:=c(90, 90, 84, 53)]


ggplot(all.genome) +
    aes(x = thr, y = prct, colour = genome) +
    geom_line() +
    scale_color_hue(direction = 1) +
    mytheme +
    xlab("% Sequence identity") +
    ylab("% of Transripts\nover threshold")+
  labs(color="", title="Personal Genome to GRCh38")+
  theme(legend.position = c(0.3, 0.3))+
  geom_point(data = annotation_data, aes(x = thr, y = meanprct), size = 1, color = "black")+
  geom_segment(data = annotation_data, aes(x = thr, y = meanprct, xend = xend, yend = yend),
               color = "black")+  # Stick lines
  geom_text(data = annotation_data, aes(x = xend, y = yend-2, label = paste0(round(meanprct, digits=1), "%")),
            size = 6*0.35, hjust = 0.5, vjust=1, color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept=90, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept=100, linetype="dashed", color="darkgrey")
ggsave("10_figures/01_plots/supp/38_map_pg/lineplot_transcriptBLAST_identity_threshold_PG2ref.pdf", dpi=700, width = 2.5, height = 2.25,  units = "in")



### From ref to PG

all <-fread("10_figures/data/blastSeqIdentity_fromRef.tsv")

# toPlot <- "T2T"
toPlot <- "hg38"

all.genome <- all[all$genome == toPlot,]

all.genome[, meanprct:=mean(prct), by=.(thr)]
annotation_data <- unique(all.genome[thr%in%c(1, 90, 99, 100), .(thr , meanprct)])
annotation_data[, xend:=c(10, 80, 90, 90)]
annotation_data[, yend:=c(90, 90, 84, 53)]


ggplot(all.genome) +
  aes(x = thr, y = prct, colour = ref) +
  geom_line() +
  scale_color_hue(direction = 1) +
  mytheme +
  xlab("% Sequence identity") +
  ylab("% of Transripts\nover threshold")+
  labs(color="", title="GRCh38 to Personal Genome")+
  theme(legend.position = c(0.3, 0.3))+
  geom_point(data = annotation_data, aes(x = thr, y = meanprct), size = 1, color = "black")+
  geom_segment(data = annotation_data, aes(x = thr, y = meanprct, xend = xend, yend = yend),
               color = "black")+  # Stick lines
  geom_text(data = annotation_data, aes(x = xend, y = yend-2, label = paste0(round(meanprct, digits=1), "%")),
            size = 6*0.35, hjust = 0.5, vjust=1, color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept=90, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept=100, linetype="dashed", color="darkgrey")
ggsave("10_figures/01_plots/supp/38_map_pg/lineplot_transcriptBLAST_identity_threshold_Ref2PG.pdf", dpi=700, width = 2.5, height = 2.25,  units = "in")


