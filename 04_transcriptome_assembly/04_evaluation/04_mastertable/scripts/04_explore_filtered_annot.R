## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/"
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

data <- fread("pclavell/04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv")
metadata <- fread("pclavell/00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,2,3,5,4,8,6,7)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
data[, filter:=factor(filter, levels=c("pass", "fail"))]
filt <- data[filter=="pass"]

# PLOT
# data[, EUR:=ifelse(CEU>0 | AJI>0, 1,
#                    ifelse(CEU==0 & AJI==0, 0, "error"))]
# data[, nonEUR:=ifelse(ITU>0 | HAC>0 | PEL>0 | LWK>0 | MPC>0 | YRI>0, 1, 0)]
# data[, eur_assembled := ifelse(EUR==1, 
#                                ifelse(nonEUR==1, "EUR & nonEUR", "EUR"), "nonEUR")]
# ggplot(data[structural_category%in%c("NIC","NNC")], aes(x=eur_assembled, y=exons, fill=eur_assembled))+
#   geom_violin(adjust=3, alpha=0.8)+
#   geom_boxplot(outliers=F, width=0.15)+
#   mytheme+
#   facet_wrap(~structural_category)+
#   labs(x="", y="# Exons")+
#   scale_fill_manual(values=c("#466995", "#76517B","#A53860"))+
#   scale_y_continuous()+
#   stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=-10))+
#   ggpubr::stat_compare_means(comparisons=list(c("EUR", "EUR & nonEUR"),c("EUR", "nonEUR")),
#                              method = "t.test", 
#                              method.args = list(alternative = "two.sided")) +
#   guides(fill="none")

  
# Transcript level
ggplot(data, aes(x=structural_category, fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Transcripts")+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=0, size=6*0.35)+
  facet_wrap(~filter, labeller = labeller(filter = c("pass" = "Passed Filters", "fail" = "Filtered Out")))
ggsave("pclavell/10_figures/01_plots/supp/05_filtering/barplotStack.filtered_trx_perSQANTI.pdf", dpi=700, width = 6.1, height = 2.5,  units = "in")


# ggplot(data, aes(x=structural_category, fill=filter))+
#   geom_bar(position="fill")+
#   scale_fill_manual(values=rev(c("darkred","#356CA1" )))+
#   labs(x="", y="# Transcripts")+
#   mytheme+
#   guides(fill="none")+
#   geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))

# gene level
ggplot(unique(data[,.(associated_geneid.v, structural_category, filter)]), aes(x=structural_category, fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Genes")+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=0, size=6*0.35)+
  facet_wrap(~filter, labeller = labeller(filter = c("pass" = "Passed Filters", "fail" = "Filtered Out")))
ggsave("pclavell/10_figures/01_plots/supp/05_filtering/barplotStack.filtered_genes_perSQANTI.pdf", dpi=700, width = 6.1, height = 2.5,  units = "in")

ggplot(unique(data[,.(associated_geneid.v, structural_category, filter)]), aes(x=structural_category, fill=filter))+
  geom_bar(position="fill")+
  scale_fill_manual(values=rev(c("darkred","#356CA1" )))+
  labs(x="", y="# Genes")+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))

# plot biotypes
data$associated_gene_biotype <- ifelse(data$associated_gene_biotype=="", "novel/ambiguous gene", data$associated_gene_biotype )
ggplot(unique(data[filter=="pass",.(associated_geneid.v, associated_gene_biotype, structural_category, filter)]), aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#F79D5C","darkgrey","#297373" ))+
  labs(x="", y="# Genes")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))
ggplot(unique(data[filter=="pass",.(isoform,associated_geneid.v, associated_gene_biotype, structural_category, filter)]), aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#F79D5C","darkgrey","#297373"))+
  labs(x="", y="# Transcripts")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))

# now exploration just on filtered annotation

ggplot(filt[, .(num_exons = factor(ifelse(exons >= 25, ">=25", as.character(exons)), levels = c(as.character(1:24), ">=25")), 
           structural_category, 
           isoform)],
  aes(x = num_exons, fill = structural_category)) +
  geom_bar() +
  mytheme +
  facet_wrap(~structural_category, scales = "free_y", nrow=2) +
  scale_fill_manual(values = colsqanti)+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  guides(fill="none")+
  labs(x="Exons/transcript", y="# Transcripts")+
  scale_x_discrete(breaks = c(as.character(seq(0, 24, by = 3)), ">=25"))  # Show only even numbers and >=25


ggplot(filt, aes(y=length, x=structural_category, fill=structural_category))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(y="Transcript length", x="")+
  annotation_logticks(sides="l")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=5))
  

filt[, protein_is_nmd:=ifelse(is.na(protein_is_nmd), "No predicted ORF", protein_is_nmd)]
filt[, detected_orf:=ifelse(protein_is_nmd=="No predicted ORF", "No", "Yes")]
filt[, associated_gene_biotype:=factor(associated_gene_biotype, levels=c("protein_coding", "lncRNA", "novel/ambiguous gene"))]
library(ggbreak)
ggplot(filt, aes(fill=detected_orf, x=associated_gene_biotype))+
  geom_bar(position="fill")+
  mytheme+
  labs(y="# Transcript", x="")+
  scale_fill_manual(values=rev(c("#356CA1", "darkred")))+
  labs(fill="Predicted ORF")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))
 # scale_y_break(c(1200, 10000), scales="free") 
ggplot(filt[detected_orf=="Yes"], aes(fill=protein_is_nmd, x=structural_category))+
  geom_bar(position="fill")+
  mytheme+
  labs(y="# Transcript", x="")+
  scale_fill_manual(values=c("#B0DAF1", "#AA6373"))+
  labs(fill="Predicted NMD")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))


myfilters <-data[, .(structural_category,f_sample=ifelse(sample_sharing<2, 1, 0),
            f_tool=ifelse(tool_sharing<2, 1, 0),
            f_length=ifelse(length<300, 1, 0),
            f_expression=ifelse(flair_max_counts<3, 1,0),
            f_recount=ifelse(sj_less_recountsupported_counts<25, 1, 0)
            )]

myfilters[,filter_sharing:=f_sample+f_tool+f_length+f_expression+f_recount]

# plot how many filters are removing each transcript
myfilters_nofsm <-myfilters[!structural_category%in%c("FSM", "ISM")][,filter_sharing:=factor(filter_sharing, levels=c(0,1,2,3,4,5))]
ggplot(myfilters_nofsm, 
       aes(x=filter_sharing, fill=structural_category))+
  geom_bar()+
  mytheme+
  labs(x="# Filters removing a Transcript", y="# Filtered Out Transcripts", fill="Structural\nCategory")+
  scale_fill_manual(values=colsqanti)+
  geom_text(stat = "count",aes(label =after_stat(count)),position=position_stack(vjust = 0.5), size=6*0.35)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7))
ggsave("pclavell/10_figures/01_plots/supp/07_filter_overlap/barplot_number_filters_removing_transcript_SQANTI.pdf", dpi=700, width = 3.3, height = 2.25,  units = "in")

colnames(myfilters) <- c("structural_category", "Sample", "Tool", "Length", "Expression", "Recount3 Support", "filter_sharing")
fwrite(myfilters, "pclavell/04_transcriptome_assembly/04_evaluation/05_mastertable/data/24222024_filtered_transcripts.tsv", quote = F, row.names = F, sep="\t")
pdf("pclavell/10_figures/01_plots/supp/07_filter_overlap/barplot_NIC_filter_overlap.pdf",  width = 3, height = 2.33,onefile=FALSE)
UpSetR::upset(myfilters[structural_category=="NIC", .(Length, Sample, Expression)],
              order.by = "freq", 
              decreasing = T,
              matrix.color= colsqanti["NIC"],
              main.bar.color=colsqanti["NIC"])
dev.off()
UpSetR::upset(myfilters[structural_category=="NNC", .(Length, Sample, Tool, `Recount3 Support`, Expression)],
              order.by = "freq", 
              decreasing = T,
              matrix.color= colsqanti["NNC"],
              main.bar.color=colsqanti["NNC"])
UpSetR::upset(myfilters[structural_category%in%c("Intergenic"), .(Length, Sample, Tool, `Recount3 Support`, Expression)],
              order.by = "freq", 
              decreasing = T,
              matrix.color= colsqanti["Intergenic"],
              main.bar.color=colsqanti["Intergenic"])

# collapse sampela and tools sharing under reproducible
myfilters2 <-data[, .(structural_category,f_reproducible=ifelse(sample_sharing<2 | tool_sharing<2, 1, 0),
                     f_length=ifelse(length<300, 1, 0),
                     f_expression=ifelse(flair_max_counts<3, 1,0),
                     f_recount=ifelse(sj_less_recountsupported_counts<25, 1, 0)
)]

UpSetR::upset(myfilters2,
              order.by = "freq", 
              decreasing = T)

# REPEAT THE POP SPECIFIC PART
popsp <- data[sample_sharing>=2 & population_sharing==1, .(AJI, CEU, YRI, LWK, MPC, ITU, HAC, PEL, isoform, structural_category, associated_geneid.v)]
popsplong <- melt(popsp, value.name = "detected", variable.name = "population", id.vars = c("isoform", "structural_category", "associated_geneid.v"))[detected>=2]
popsplong <- unique(metadata[, .(population, map_reads_assemblymap)][, .(popdepth=sum(map_reads_assemblymap), samples_per_pop=.N), by="population"])[popsplong, on="population"]
popsplong[, trx_per_pop_per_cat:=.N, by=c("population", "structural_category")]
popsplong[, associated_gene:=ifelse(structural_category%in%c("FSM", "ISM", "NIC", "NNC"), "known", "novel/ambiguous" )]
popsplong[, eurpop:= ifelse(population%in%c("CEU", "AJI"), "EUR", "nonEUR")]
ggplot(popsplong, aes(x=reorder(population, popdepth), fill=structural_category))+
  geom_bar(position="fill")+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  labs(x="Population specific transcripts", y="# Transcripts")
ggplot(unique(popsplong[, .(population, eurpop, structural_category, popdepth, samples_per_pop,trx_per_pop_per_cat)]), 
       aes(x=popdepth, y=trx_per_pop_per_cat))+
  # geom_rect(aes( fill = factor(samples_per_pop), xmin = popdepth-1000000, xmax=popdepth+1000000, ymin=0, ymax=300), 
  #           alpha = 0.8) +  # Background tiles
  # scale_fill_manual(values = c("4" = "lightblue", 
  #                              "5" = "lightgreen", 
  #                              "6" = "lightcoral"),
  #                   na.value = "white")+
  # ggnewscale::new_scale_color()+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total population mapped reads", y="# Transcripts")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(5,3))+
  scale_alpha_manual(values=c(0.6, 0.9))+
  ylim(c(0, 265))
ggplot(popsplong[!population%in%c("CEU", "AJI")], aes(x=popdepth, y=trx_per_pop_per_cat,col=structural_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total population mapped reads", y="# Transcripts")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c( "grey"))+
  scale_size_manual(values=c(3))+
  scale_alpha_manual(values=c( 0.9))+
  ylim(c(0, 265))


popspgenes <-unique(popsplong[, .(population, popdepth, associated_geneid.v, eurpop)])[, howManyPopHaveTheGene:=.N, by="associated_geneid.v"][howManyPopHaveTheGene==1]
popspgenes[, countGenes :=.N, by="population"]
ggplot(unique(popspgenes[eurpop=="nonEUR", .(popdepth, countGenes, eurpop)]), aes(x=popdepth, y=countGenes))+
  geom_line(linewidth=1.5)+
  mytheme+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c( "grey"))+
  scale_size_manual(values=c(3))+
  scale_alpha_manual(values=c(0.9))+
  ylim(c(0, 375))+
  labs(x="Total population mapped reads", y="# Population Specific Genes")

ggplot(unique(popspgenes[, .(popdepth, countGenes, eurpop)]), aes(x=popdepth, y=countGenes))+
  geom_line(linewidth=1.5)+
  mytheme+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("red", "grey"))+
  scale_size_manual(values=c(5,3))+
  scale_alpha_manual(values=c(0.6, 0.9))+
  ylim(c(0, 375))+
  labs(x="Total population mapped reads", y="# Population Specific Genes")


# replot blasttscore
ggplot(filt[, associated_gene_biotype:=ifelse(associated_gene_biotype=="", "novel/ambiguous", associated_gene_biotype)], aes(y=blastp_identity, x=associated_gene_biotype, fill=associated_gene_biotype))+
  geom_violin(alpha=0.7, scale="width")+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=4))+
  scale_fill_manual(values=c( "#F79D5C","darkgrey","#297373" ))+
  facet_wrap(~structural_category)


### EXPLORE THOSE 800 NNC TRANSCRIPTS THAT ARE ONLY FILTERED OUT BY RECOUNT FILTER
data[, `:=`(f_sample=ifelse(sample_sharing<2, 1, 0),
         f_tool=ifelse(tool_sharing<2, 1, 0),
         f_length=ifelse(length<300, 1, 0),
         f_expression=ifelse(flair_max_counts<3, 1,0),
         f_recount=ifelse(sj_less_recountsupported_counts<25, 1, 0)
)]

# identify said transcripts
onlyrecountfilt <-data[f_recount==1 & f_sample==0 & f_tool==0 & f_expression==0 & f_length==0]

# identify in which populations are assembled
onlyrecountfiltlong <-melt(onlyrecountfilt, measure.vars = c("YRI", "LWK", "MPC", "CEU", "AJI", "ITU", "HAC", "PEL"), variable.name = "population", value.name = "detected")[detected>0]
onlyrecountfiltlong <- unique(metadata[, popdepth:=sum(map_reads_assemblymap), by="population"][, .(popdepth, population)])[onlyrecountfiltlong, on="population"]
onlyrecountfiltlong2 <- unique(onlyrecountfiltlong[filter=="fail" & !structural_category%in%c("FSM", "ISM")][, countperpopcat:=.N, by=c("population", "structural_category")][, eurpop:=ifelse(population%in%c("CEU", "AJI"), "EUR", "nonEUR")][, .(population, eurpop, popdepth, countperpopcat, structural_category)])

# intersection between ancestries removing these transripts based on just recount
filteredNNC <-onlyrecountfilt[structural_category=="NNC", .(AJI, CEU, YRI, LWK, MPC, ITU, HAC, PEL)]

filteredNNC[, c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL") := 
              lapply(.SD, function(x) ifelse(x != 0, 1, 0)), 
            .SDcols = c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL")]
UpSetR::upset(filteredNNC,
              order.by = "freq", 
              decreasing = T,
              nsets=8,
              matrix.color= colsqanti["NNC"],
              main.bar.color=colsqanti["NNC"])
# now with intergenics
filteredNNC <-onlyrecountfilt[structural_category=="Intergenic", .(AJI, CEU, YRI, LWK, MPC, ITU, HAC, PEL)]

filteredNNC[, c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL") := 
              lapply(.SD, function(x) ifelse(x != 0, 1, 0)), 
            .SDcols = c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL")]
UpSetR::upset(filteredNNC,
              order.by = "freq", 
              decreasing = T,
              nsets=8,
              matrix.color= colsqanti["Intergenic"],
              main.bar.color=colsqanti["Intergenic"])
# now with the rest
filteredNNC <-onlyrecountfilt[structural_category%in%c("Fusion","Intergenic", "Antisense", "Genic"), .(AJI, CEU, YRI, LWK, MPC, ITU, HAC, PEL)]

filteredNNC[, c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL") := 
              lapply(.SD, function(x) ifelse(x != 0, 1, 0)), 
            .SDcols = c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL")]
UpSetR::upset(filteredNNC,
              order.by = "freq", 
              decreasing = T,
              nsets=8,
              matrix.color="black",
              main.bar.color="black")

ggplot(onlyrecountfiltlong2, 
       aes(x=popdepth, y=countperpopcat,col=structural_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("red", "grey"))+
  scale_size_manual(values=c(5,3))+
  scale_alpha_manual(values=c(0.6, 0.9))+
  labs(x="Total population throughput", y="# Detected (dropped) trancripts\nby recount3 filter")


#### Compare recounts between ancestries to discoveer if there is any bias
datalong <-melt(data, measure.vars = c("YRI", "LWK", "MPC", "CEU", "AJI", "ITU", "HAC", "PEL"), variable.name = "population", value.name = "detected")[detected>0]
datalong <- unique(metadata[, popdepth:=sum(map_reads_assemblymap), by="population"][, .(popdepth, population)])[datalong, on="population"]
datalong2 <- unique(datalong[ !structural_category%in%c("FSM", "ISM")][, countperpopcat_allunfilteredannot:=.N, by=c("population", "structural_category")][, eurpop:=ifelse(population%in%c("CEU", "AJI"), "EUR", "nonEUR")][, .(population, eurpop, popdepth, countperpopcat_allunfilteredannot, structural_category)])

merged <- onlyrecountfiltlong2[datalong2, on=c("population", "eurpop", "popdepth", "structural_category")]
mergedlong <- melt(merged, measure.vars = c("countperpopcat", "countperpopcat_allunfilteredannot"), value.name = "count", variable.name = "level")
mergedlong[, level:=ifelse(level=="countperpopcat", "Filtered by recount", "Whole unfiltered annotation")]

ggplot(mergedlong,
       aes(x=countperpopcat_allunfilteredannot, y=countperpopcat, col=population))+
  geom_point()+
  scale_color_manual(values=popcol)+
  mytheme+
  facet_wrap(~structural_category, scales="free_y")+
  geom_abline(intercept = 0, slope=0.01)


ggplot(mergedlong[structural_category!="NIC"], 
       aes(x=popdepth, y=count,col=structural_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("red", "grey"))+
  scale_size_manual(values=c(5,3))+
  scale_alpha_manual(values=c(0.6, 0.9))+
  labs(x="Total population throughput", y="# Transcripts")+
  facet_wrap(~level, scales="free_y")

# do density plots on geom_ridge to visualize the recount support

ggplot(unique(datalong[!structural_category%in%c("FSM", "ISM", "NIC"), .(isoform,population, structural_category, sj_less_recountsupported_counts, popdepth)]),
       aes(x=log10(sj_less_recountsupported_counts+1),y=reorder(population, popdepth), fill=population))+
  ggridges::geom_density_ridges(alpha=0.5)+
  mytheme+
  scale_fill_manual(values=popcol)+
  facet_wrap(~structural_category)+
  annotation_logticks(sides="b")+
  ggtitle("Least supported SJ")+
  labs(x="Recount3 support of Least supported SJ/Transcript", y="")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=9), vjust=-0.5)+
  geom_vline(xintercept=log10(25+1), linetype="dashed", color="darkgrey")


ggplot(unique(datalong[!structural_category%in%c("FSM", "ISM", "NIC"), .(isoform,population, structural_category, sj_less_recountsupported_novel_counts, popdepth)]),
       aes(x=log10(sj_less_recountsupported_novel_counts+1),y=reorder(population, popdepth), fill=population))+
  ggridges::geom_density_ridges(alpha=0.5)+
  mytheme+
  scale_fill_manual(values=popcol)+
  facet_wrap(~structural_category)+
  annotation_logticks(sides="b")+
  ggtitle("Least supported novel SJ")+
  labs(x="Recount3 support of Least supported NOVEL SJ/Transcript", y="")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=9), vjust=-0.5)+
  geom_vline(xintercept=log10(25+1), linetype="dashed", color="darkgrey")



ggplot(unique(datalong[sample_sharing>1 & population_sharing==1][!structural_category%in%c("FSM", "ISM", "NIC"), .(isoform, population, structural_category, popdepth,sj_less_recountsupported_counts)]),
       aes(x=log10(sj_less_recountsupported_counts+1),y=reorder(population, popdepth), fill=population))+
  ggridges::geom_density_ridges(alpha=0.5)+
  mytheme+
  scale_fill_manual(values=popcol)+
  facet_wrap(~structural_category)+
  annotation_logticks(sides="b")+
  ggtitle("Least supported SJ (in Population Specific Transcripts)")+
  labs(x="Recount3 support of Least supported SJ/pop specific Transcript", y="")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=9), vjust=-0.5)+
  geom_vline(xintercept=log10(25+1), linetype="dashed", color="darkgrey")



ggplot(unique(datalong[sample_sharing>1 & population_sharing==1][!structural_category%in%c("FSM", "ISM", "NIC"), .(isoform,population, structural_category, sj_less_recountsupported_novel_counts, popdepth)]),
       aes(x=log10(sj_less_recountsupported_novel_counts+1),y=reorder(population, popdepth), fill=population))+
  ggridges::geom_density_ridges(alpha=0.5)+
  mytheme+
  scale_fill_manual(values=popcol)+
  facet_wrap(~structural_category)+
  annotation_logticks(sides="b")+
  ggtitle("Least supported novel SJ (in Population Specific Transcripts)")+
  labs(x="Recount3 support of Least supported NOVEL SJ/pop specific Transcript", y="")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=9), vjust=-0.5)+
  geom_vline(xintercept=log10(25+1), linetype="dashed", color="darkgrey")


# check if intra sharing is higher than across population sharing
dataupset <-data[sample_sharing>1, .(AJI,CEU, YRI, MPC, LWK, ITU, HAC, PEL)]
dataupset[, c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL") := 
              lapply(.SD, function(x) ifelse(x != 0, 1, 0)), 
            .SDcols = c("AJI", "CEU", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL")]
UpSetR::upset(dataupset,
              order.by = "freq", 
              decreasing = T,
              nsets=8,
              matrix.color="black",
              main.bar.color="black")

### Structural categories by tool

ggplot(melt(filt, measure.vars = c("flair", "isoquant", "lyric", "espresso"), variable.name = "tool", value.name = "detected")[detected==1][, known_associated:=ifelse(structural_category%in%c("FSM", "ISM", "NIC", "NNC"), "known gene", "novel/ambiguous gene")], aes(x=tool, fill=structural_category))+
  mytheme+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  labs(x="", y="# Transcripts")


##### POPULATION SPECIFIC TRANSCRIPTS
# LINE PLOT FOR 2 SAMPLES SHARING FACETED
popsplong[, eur:=fifelse(eurpop=="EUR", "European", "Non-European")]
popsplong[, total_throughput:=popdepth]
popsplong[, trx_per_cat_per_pop:=trx_per_pop_per_cat]
subpopsplongmeta<- popsplong

ggplot(unique(popsplong[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), .(population, eur, structural_category, total_throughput, samples_per_pop,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(data=unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC")&eur=="Non-European", 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
            linewidth=1, lty="11", color="darkgrey")+
  geom_line(aes(col=structural_category),linewidth=1)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific UMA Transcripts")+
  guides(color="none")+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  # geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("FSM"), 
  #                                           .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
  #              aes(xend = total_throughput/10^6, # Adjust for arrow length
  #                  yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
  #              arrow = arrow(length = unit(0.2, "cm")), # Arrow size
  #              color = "darkgrey") +
  # geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"), 
  #                                           .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
  #              aes(xend = total_throughput/10^6, # Adjust for arrow length
  #                  yend = ifelse(eur == "European", trx_per_cat_per_pop - 25, trx_per_cat_per_pop + 25)),
  #              arrow = arrow(length = unit(0.2, "cm")), # Arrow size
  #              color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop + 45, trx_per_cat_per_pop - 45),
                label = population,
                color=eur),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=c("#466995", "#A53860"))+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop - 45, trx_per_cat_per_pop + 45),
                label = population,
                color=eur),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=c("#466995", "#A53860"))+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.90))+
  xlim(c(52, 105))+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(4,1.5))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(alpha="", size="", color="")
ggsave("pclavell/10_figures/01_plots/supp/15_uma_popsp/line_UMA_popSpecific_transcripts.2samplesSharing.faceted.pdf", dpi=700, width = 6, height = 6,  units = "in")
ggsave("pclavell/10_figures/02_panels/svg/supp/15_uma_popsp.svg", dpi=700, width = 6, height = 6,  units = "in")
ggsave("pclavell/10_figures/02_panels/png/supp/15_uma_popsp.png", dpi=700, width = 6, height = 6,  units = "in")


# LINE PLOT FOR 2 SAMPLES SHARING FACETED with ----- REFSEQ ----------
refseq <- fread("pclavell/04_transcriptome_assembly/04_evaluation/02_sqanti/data/uma_evaluatedBy_RefSeq/uma_evaluatedBy_RefSeq_classification.txt")
refseq <- refseq[, .(isoform, structural_category)]
categories <-c("Antisense","Intergenic", "NNC", "FSM", "NIC", "ISM", "Genic", "Fusion" )
names(categories) <- unique(refseq$structural_category)
refseq[, structural_category_refseq:=categories[structural_category]]
subpopsplongmeta_refseq <- refseq[, .(isoform, structural_category_refseq)][subpopsplongmeta, on="isoform"]
subpopsplongmeta_refseq[, trx_per_cat_per_pop:=uniqueN(isoform), by=c("structural_category_refseq", "population")]


ggplot(unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM", "NIC", "NNC"),][,.(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM", "NIC", "NNC")&eur=="Non-European", 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]), 
            linewidth=1, lty="11", color="darkgrey")+
  geom_line(aes(col=structural_category_refseq),linewidth=1)+
  mytheme+
  scale_color_manual(values=c(colsqanti, "ISM"="#8EDE95"))+
  labs(x="Total Mapped Reads per Population (M)", y="# Population-Specific UMA Transcripts")+
  guides(color="none")+
  facet_wrap(~structural_category_refseq)+
  labs(col="", alpha="", size="")+
  geom_segment(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM"),
                                                   .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
               color = "darkgrey") +
  geom_segment(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("NIC", "NNC"),
                                                   .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop - 25, trx_per_cat_per_pop + 25)),
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM"), 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop + 35, trx_per_cat_per_pop - 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  geom_text(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("NIC", "NNC"), 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop - 35, trx_per_cat_per_pop + 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.90))+
  xlim(c(52, 105))+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(4,1.5))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(alpha="", size="", color="")
ggsave("pclavell/10_figures/01_plots/supp/18_popsp_val_refseq_uma/line_UMA_popSpecific_transcripts.2samplesSharing.faceted_refseq_UMA.pdf", dpi=700, width = 7, height = 6,  units = "in")
ggsave("pclavell/10_figures/02_panels/svg/supp/18_popsp_val_refseq_uma.svg", dpi=700, width = 7, height = 6,  units = "in")
ggsave("pclavell/10_figures/02_panels/png/supp/18_popsp_val_refseq_uma.png", dpi=700, width = 7, height = 6,  units = "in")
