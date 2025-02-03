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
vec <-list.files("data/fabien_data/read_diffs/")
filtered_vec <- vec[!grepl("coverage", vec)]


epsilon <- 1e-10  # A very small value
breaks <- c(-Inf, -8, -1, 0, 8, Inf)  # Adjusted to create 5 intervals
# Create a factor by cutting the numeric vector into intervals
factor_vector <- cut(num_vector,
                     breaks = breaks,
                     labels = c("(-Inf, -8]", "(-8, -1]", "(-1, 0)", "[0, 0]", "(0, 1)", "(1, 8]", "(8, Inf)"),
                     include.lowest = TRUE)
for(myfile in filtered_vec){
  print(myfile)
  data <- fread(paste0("data/fabien_data/read_diffs/", myfile, "/missmatchesGapCovDiffTable.tsv"), sep="\t")
  data[, `:=`(count=as.numeric(gsub(" .*", "", V1)), misgap=as.numeric(gsub(".* ", "", V1)), cov=V2)][, V1:=NULL][, V2:=NULL]

  # # create categories on the diff of coverage %
  # data[, newcoverage :=fifelse(cov==0, "No Differences",
  #                                  fifelse(cov<0 & cov>=-0.01, "Up to 1% less",
  #                                          fifelse(cov<(-0.01), "More than 1% less",
  #                                                  fifelse(cov>0 & cov<=0.01, "Up to 1% more",
  #                                                          fifelse(cov>0.01, "More than 1% more", "Error")))))]
  # data[, newcoverage:=factor(newcoverage, levels=c("More than 1% less", "Up to 1% less", "No Differences", "Up to 1% more","More than 1% more"))]
# 
# 
#   # Create a factor by cutting the numeric vector into intervals
#   data[, newmisgap:=cut(misgap, breaks = breaks,
#                         labels = c("(-Inf, -8]", "[-8, -1]",  "[0, 0]",  "[1, 8]", "[8, Inf)"),
#                         include.lowest = T)]
#   data[, newcount:=sum(count), by=.(newcoverage, newmisgap)]
  data[, misgaps:=fifelse(misgap==0, "Equal", 
                         fifelse(misgap>0, "More", "Less"))]
  data[, coverage:=fifelse(cov==0, "Equal", 
                            fifelse(cov>0, "Higher", "Lower"))]
  data[, data:=myfile]
  fwrite(data, paste0("data/fabien_data/read_diffs/", myfile, "coverage_mismatches_gaps_summary.tsv"))
}

# Summaries have been generated now merged them
files <- list.files("data/fabien_data/read_diffs/", pattern="summary")
data <- list()
for(file in files){
  subdata <- fread(paste0("data/fabien_data/read_diffs/",file))
  data <- append(data, list(as.data.frame(subdata)))
}
data <- rbindlist(data)
data[, `:=`(sample=gsub("_.*", "", data), ref=gsub(".*_", "", data))]
data[ref=="hg38", ref:="GRCh38"]

data[, newcount:=sum(count), by=.(misgaps, coverage, sample, ref)]
data[, total:=sum(count), by=.(sample, ref)]
data[, percent:=round(newcount/total*100, 2)]
mydata <- unique(data[, .(percent, misgaps, coverage, sample, ref)])
mydata[, mean_percent:=round(mean(percent), 2), by=.(misgaps, coverage, ref)]


mydata[, `:=`(coverage=factor(coverage, levels=c("Lower", "Equal", "Higher")),
              misgaps=factor(misgaps, levels=c("Less", "Equal", "More")))]
ggplot(unique(mydata[ref=="GRCh38", .(misgaps, coverage, mean_percent, ref)]), 
       aes(x=coverage, y=misgaps, fill=mean_percent))+
  geom_tile()+
  geom_text(aes(label=paste0(mean_percent, " %")), size=6*0.35)+
  scale_fill_gradient(trans="log10", low="white", high = "#3c9eb6ff")+
  mytheme+
  labs(x="Read Coverage\nin personal assembly\ncompared to GRCh38",
       y="# of Mismatches and Gaps\nin personal assembly\ncompared to GRCh38",
       fill="Mean %\nof Reads")
ggsave("01_plots/main/fig_05/heatmap_diff_in_coverage_and_misgaps_reads_GRCh38.pdf", dpi=700, width = 3, height = 2.25,  units = "in")
ggplot(unique(mydata[ref=="T2T", .(misgaps, coverage, mean_percent, ref)]), 
       aes(x=coverage, y=misgaps, fill=mean_percent))+
  geom_tile()+
  geom_text(aes(label=paste0(mean_percent, " %")), size=6*0.35)+
  scale_fill_gradient(trans="log10", low="white", high = "#3c9eb6ff")+
  mytheme+
  labs(x="Read Coverage\nin personal assembly\ncompared to T2T",
       y="# of Mismatches and Gaps\nin personal assembly\ncompared to T2T",
       fill="Mean %\nof Reads")
ggsave("01_plots/supp/38_map_pg/hetmap_diff_in_coverage_and_misgaps_reads_T2T.pdf", dpi=700, width = 3, height = 2.25,  units = "in")




# data[, meancount:=paste0(round(mean(newcount)/1000, digits=0), "k"), by=.(newcoverage, newmisgap)]
# data[, sd:=sd(perc), by=.(newcoverage, newmisgap)]
# data[, meanper:=paste0(round(mean, 1), "%")]
# mydata <- unique(data[, .(newcoverage, newmisgap, meancount, mean, meanper, sd)])

# mydata[, newcoverage:=factor(newcoverage, levels=rev(c("More than 1% less", "Up to 1% less", "No Differences", "Up to 1% more","More than 1% more")))]

# ggplot(unique(mydata[, .(newcoverage, mean, newmisgap, meanper, sd,meancount)]), aes(x=as.factor(newcoverage), y=mean, fill=newmisgap))+
#   geom_col(position=position_dodge(width=0.85))+
#   geom_errorbar(aes(ymin = mean-sd, ymax=mean+sd), width=0.25, position=position_dodge(width=0.85))+
#   mythemen+
#   scale_fill_manual(values=c("darkred", "red", "grey", "blue", "darkblue"))+
#   scale_y_continuous(trans="log10")+
#   geom_text(aes(label=meanper, y=100), position=position_dodge(width=0.85), vjust=0.5, angle=90, hjust=0)+
#   geom_text(aes(label=meancount, y=1000), position=position_dodge(width=0.85), vjust=0.5, angle=90, hjust=0)+
#   coord_cartesian(ylim=c(0.01, 5000))+
#   labs(x="Read Coverage Difference\n(Personal Assembly - GRCh38)", y="Mean % of Reads", fill="Difference in\n# of Mismatches\nand Gaps")+
#   scale_x_discrete(labels=c("More than 1% less"=">1%",
#                           "Up to 1% less"="0-1%",
#                           "No Differences"="0",
#                           "Up to 1% more"="0-1%",
#                           "More than 1% more"=">1%"))



ggplot(mydata, aes(x=newcoverage, y=mean))+
  geom_col(, fill="#4A9F37")+
  geom_errorbar(aes(ymin = mean-sd, ymax=mean+sd), width=0.25)+
  mytheme+
  labs(x="Coverage Difference \n(Personal Assembly - GRCh38)", y="% Reads")+
  geom_text(aes(label=paste0(round(mean, 2), "%")), size=6*0.35, vjust=-0.5)+
  scale_x_discrete(labels=c("More than 1% less"=">1%",
                            "Up to 1% less"="0-1%",
                            "No Differences"="0",
                            "Up to 1% more"="0-1%",
                            "More than 1% more"=">1%"))+
  scale_y_continuous(limits=c(0,100))
ggsave("01_plots/main/fig_05/barplot_diff_in_coverage_reads.pdf", dpi=700, width = 2, height = 2.25,  units = "in")

