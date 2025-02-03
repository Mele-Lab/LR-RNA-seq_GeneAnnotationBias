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

data <- fread("10_figures/data/gene_liftoff/combined_output.txt")
colnames(data) <- c("gtf", "gene", "coverage", "identity")

data[, `:=`(sample=gsub("_.*", "", gtf),
            ref=gsub("\\.gtf","",gsub(".*_", "", gtf)))]
data[ref=="hg38", ref:="GRCh38"]


data_cdf <- data[, .(count = .N), by = .(sample, ref, identity)][order(sample, ref, identity)]
data_cdf[, cum_pct := (rev(cumsum(rev(count))) / sum(count)) * 100, by = .(sample, ref)]

data_cdf[, min:=min(cum_pct), by=.(sample, ref)]
mean(unique(data_cdf[, min]))
# Plot cumulative percentage
main <- ggplot(data_cdf[ref=="GRCh38"], aes(x = identity, y = cum_pct, col=sample)) +
  geom_step() +
  labs(col="",
       x = "Gene Sequence Identity",
       y = "% Genes over Identity Threshold") +
  mythemen+
  scale_color_manual(values=rev(c("#3C9EB6","#3D8FBC","#597EB9", "#7869A9", "#90528D", "#9C3C67")))+
  coord_cartesian(xlim=c(0.95, 1), ylim=c(0,100))

# Do boxplots
data[, `:=`(
  total = .N,
  thresh1=sum(identity >= 1),
  thresh999 = sum(identity >= 0.999),  # Count of minparam >= 0.999 in each group
  thresh99 = sum(identity >= 0.99),    # Count of minparam >= 0.99 in each group
  thresh95 = sum(identity >= 0.95)     # Count of minparam >= 0.95 in each group
), by = .(sample, ref)]

datalong <-melt(unique(data[, .(ref, sample, thresh999,thresh1, thresh99, thresh95, total)]), id.vars = c("ref", "sample", "total"))
datalong[, percent:=round(value/total*100,2)]

sub <-ggplot(datalong, aes(x=variable, y=percent))+
  geom_boxplot()+
  mythemen+
  labs(x="Identity Threshold", y="% Lifted off Genes")+
  scale_x_discrete(labels=c( "thresh999"="99.9%", "thresh99"="99%", "thresh95"="95%", "thresh1"="100%"), limits=rev(c("thresh1","thresh999", "thresh99", "thresh95")))+
  coord_cartesian(ylim=c(0,100))

# put plot inside
library(cowplot)
ggdraw() +
  draw_plot(main) +
  draw_plot(sub, x = 0.15, y = .15, width = .6, height = .6)
