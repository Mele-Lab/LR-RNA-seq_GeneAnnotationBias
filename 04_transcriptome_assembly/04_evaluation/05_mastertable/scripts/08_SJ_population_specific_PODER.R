## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation"
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

# load data
data <- fread("05_mastertable/data/29102024_PODER_mastertable.tsv")
data <- data[structural_category%in%c("FSM","NIC" ,"NNC")]
sj <- fread("02_sqanti/data/poder/poder_junctions.txt")
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- unique(metadata$color_pop)
names(popcol) <- unique(metadata$population)

# merge
sj[, junction:=paste(chrom, strand, genomic_start_coord,genomic_end_coord, sep="_")]
sj[, sj_category := fifelse(junction_category == "known", "Known SJ",
                            fifelse(start_site_category == "known" & end_site_category == "known", "Novel SJ: Known SS",
                                    fifelse(xor(start_site_category == "known", end_site_category == "known"), 
                                            "Novel SJ: 1 Novel SS", 
                                            "Novel SJ: 2 Novel SS")))]

sj1 <- data[,.(isoform, structural_category, AJI, CEU, MPC, YRI, LWK, HAC, ITU, PEL, population_sharing, sample_sharing)][sj, on="isoform"]

# keep only sj belonging to population specific transcripts
sjone <- sj1[population_sharing==1 & sample_sharing>=2,]
sjonelong <- melt(sjone, measure.vars = c("AJI", "CEU", "MPC", "YRI", "LWK", "HAC", "ITU", "PEL"), value.name = "detected", variable.name = "population")[detected>=2]







#### PARENTHESIS TO TRY TO IDENTIFY POP SPECIFIC JUNCTIONS------------------

# SJ from pop specific trx
sjonelong[, sj_total_ocurrences:=.N, by="junction"]
sjonelong[, sj_pop_ocurrences:=.N, by=c("junction", "population")]
sjonelong[, sj_pop_specific:=ifelse(sj_total_ocurrences==sj_pop_ocurrences, "Pop Specific", "Pop Shared")]
sjonelong <- unique(metadata[, .(population, map_reads_assemblymap)][, .(popdepth=sum(map_reads_assemblymap)), by="population"])[sjonelong, on="population"]
sjonelong[, sjcount:=uniqueN(junction), by=c("sj_category", "population")]
sjonelong[,  eurpop:= ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]
sjonelong[, sjcategory:=factor(sj_category, levels = c("Known SJ", "Novel SJ: Known SS", "Novel SJ: 1 Novel SS", "Novel SJ: 2 Novel SS"))]
ggplot(unique(sjonelong[, .(popdepth, sjcount,sj_category, eurpop, population)]), aes(x=popdepth/10^6, y=sjcount, col=sj_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=c("#8DA04E","#F3D9B1", "#C29979", "#A22522"))+
  labs(x="Total population mapped reads", y="# Population Specific Splice Junctions", col="")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("#466995", "#A53860"), name="")+
  scale_size_manual(values=c(5,3), name="")+
  scale_alpha_manual(values=c(0.6, 0.9), name="")+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts' Splice Junctions")+
  geom_segment(data=unique(sjonelong[sj_category=="Known SJ", .(popdepth, sjcount,sj_category, eurpop, population)]),
               aes(xend = popdepth/10^6, # Adjust for arrow length
                   yend = ifelse(eurpop == "European", sjcount + 125, sjcount - 125)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(sjonelong[sj_category=="Known SJ", .(popdepth, sjcount,sj_category, eurpop, population)]),
            aes(x = popdepth/10^6,  # Position text away from point
                y = ifelse(eurpop == "European", sjcount + 175, sjcount - 175),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 4.5,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.85, 0.45))
ggsave("../../10_figures/suppfig/line_PODER_popSpecific_trxSJ.2samplesSharing.pdf", dpi=700, width = 18, height = 17,  units = "cm")

# Calculate the proportion of "Known SJ" per population
known_proportion <- sjonelong[, .(known_percentage = mean(sj_category == "Known SJ")), by = population]
# Reorder population by known_percentage
sjonelong[, population := factor(population, levels = known_proportion[order(-known_percentage)]$population)]
sjonelong[, sj_category:=factor(sj_category, levels = c("Known SJ", "Novel SJ: Known SS", "Novel SJ: 1 Novel SS", "Novel SJ: 2 Novel SS"))]
toplot <- unique(sjonelong[, .(junction,sj_category, population)])
toplot[, sjcategory:=factor(sjcategory, levels = c("Known SJ", "Novel SJ: Known SS", "Novel SJ: 1 Novel SS", "Novel SJ: 2 Novel SS"))]

ggplot(unique(sjonelong[, .(junction,sj_category, population)]), aes(x=population, fill=sj_category))+
  geom_bar(position="fill")+
  labs(x="", y="Proportion of Population-Specific\nSplice Junctions", fill="")+
  scale_fill_manual(values=c("#8DA04E","#F3D9B1", "#C29979", "#A22522"),
                    labels = c("Known SJ"="Known SJ", "Novel SJ: Known SS"="Novel SJ:\nKnown SS", "Novel SJ: 1 Novel SS"="Novel SJ:\n1 Novel SS", "Novel SJ: 2 Novel SS"="Novel SJ:\n2 Novel SS"))+
  geom_text(
    aes(label = paste0(round(after_stat(count) / tapply(after_stat(count), after_stat(x), sum)[after_stat(x)], digits = 3)*100, "%"),
        alpha=sj_category),
    stat = "count",
    position = position_fill(vjust = 0.5),
    size=6*0.35
  )+
  scale_alpha_manual(values=c(1, 0, 1, 0))+
  guides(alpha="none")+
  mytheme+
  theme(legend.position = "top",
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.x = unit(1, "cm"),           # Reduce horizontal spacing between keys
        legend.text = element_text(margin = margin(r = -1, unit = "pt")),
        legend.box.margin=margin(-10,-10,-10,-10)
  )
ggsave("../../10_figures/fig_03/barplot_PODER_popSpecificTrx_SJ.bycategory.pdf", dpi=700, width = 3, height = 2.25,  units = "in")

popspsj <-sjonelong[sj_pop_specific=="Pop Specific"]
popspsj <- unique(metadata[, .(population, map_reads_assemblymap)][, .(popdepth=sum(map_reads_assemblymap)), by="population"])[popspsj, on="population"]
popspsj[, popspsj_x_pop_x_cat:=.N, by=c("population", "sj_category")]
popspsj[,  eurpop:= ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]

ggplot(unique(popspsj[, .(popdepth, popspsj_x_pop_x_cat,sj_category, eurpop, population)]), aes(x=popdepth/10^6, y=popspsj_x_pop_x_cat, col=sj_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=c("#8DA04E","#F3D9B1", "#C29979", "#A22522"))+
  labs(x="Total population mapped reads", y="# Population Specific Splice Junctions", col="")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("#466995", "#A53860"), name="")+
  scale_size_manual(values=c(5,3), name="")+
  scale_alpha_manual(values=c(0.6, 0.9), name="")+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Splice Junctions")+
  geom_segment(data=unique(popspsj[sj_category=="Known SJ", .(popdepth, popspsj_x_pop_x_cat,sj_category, eurpop, population)]),
    aes(xend = popdepth/10^6, # Adjust for arrow length
        yend = ifelse(eurpop == "European", popspsj_x_pop_x_cat + 125, popspsj_x_pop_x_cat - 125)),
    arrow = arrow(length = unit(0.2, "cm")), # Arrow size
    color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(popspsj[sj_category=="Known SJ", .(popdepth, popspsj_x_pop_x_cat,sj_category, eurpop, population)]),
    aes(x = popdepth/10^6,  # Position text away from point
        y = ifelse(eurpop == "European", popspsj_x_pop_x_cat + 175, popspsj_x_pop_x_cat - 175),
        label = population,
        color=population),
    fontface="bold",
    family="Helvetica",
    size = 4.5,
    alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.85, 0.45))
ggsave("../../10_figures/suppfig/line_PODER_popSpecific_SJ.2samplesSharing.pdf", dpi=700, width = 18, height = 17,  units = "cm")
# Calculate the proportion of "Known SJ" per population
known_proportion <- popspsj[, .(known_percentage = mean(sj_category == "Known SJ")), by = population]
# Reorder population by known_percentage
popspsj[, population := factor(population, levels = known_proportion[order(-known_percentage)]$population)]
popspsj[, sj_category:=factor(sj_category, levels=c("Known SJ", "Novel SJ: Known SS", "Novel SJ: 1 Novel SS", "Novel SJ: 2 Novel SS"))]
ggplot(unique(popspsj[, .(junction,sj_category, population)]), aes(x=population, fill=sj_category))+
  geom_bar(position="fill")+
  mytheme+
  labs(x="", y="Proportion of Population-Specific\nSplice Junctions by Category", fill="")+
  scale_fill_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  geom_text(
    aes(label = round(after_stat(count) / tapply(after_stat(count), after_stat(x), sum)[after_stat(x)], digits = 3),
        alpha=sj_category),
    stat = "count",
    position = position_fill(vjust = 0.5)
  )+
  scale_alpha_manual(values=c(1, 0, 0, 0))+
  guides(alpha="none")
ggsave("../../10_figures/suppfig/barplot_PODER_popSpecific_SJ.bycategory.pdf", dpi=700, width = 18, height = 10,  units = "cm")

# PERMUTATION TEST -------I'll test proportions, are CEU having a largest proportion of knownSJ within their population specific SJ ??
pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(data)[grepl(pattern, colnames(data))]

# keep only transcripts that are shared at least across 2 samples because (if not, they won't pass the popsp filter)
sj2 <- data[, .SD, .SDcols=c(samplecols, "isoform", "sample_sharing", "population_sharing")][sj, on="isoform"]
sj2 <- sj2[sample_sharing>=2,]

# sj2 <- sj[sample_sharing>1][, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]
sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]

samplingvec <- gsub(".$", "", unique(sjlong$sample))

resvecpos <- c()
resvecprop <- c()
for(iteration in 1:1000){
  print(paste0("Iteration ", iteration))
  randompop <-sample(samplingvec, 43)
  sjlong3 <- copy(sjlong)
  names(randompop) <- unique(sjlong3$sample)
  
  sjlong3[, `:=`(population=randompop[sample])]
  sjlong3 <- unique(metadata[, .(sample, map_reads_assemblymap)])[sjlong3, on="sample"]
  sjlong3[, popreads := sum(unique(map_reads_assemblymap)), by="population"]
  
  # detect pseupop sp transcripts
  # count ocurrences of isoform by population (how many samples in a pop detect the transcript)
  pseudopopsp <- unique(unique(sjlong3[, .(isoform, population, sample)])[, samples_x_pseudopop:=.N, by=c("isoform", "population")][samples_x_pseudopop>1][, .(isoform, population, samples_x_pseudopop)])
  # keep pop specific transcripts
  pseudopopsp <- pseudopopsp[, popsharing := .N, by="isoform"][popsharing==1]
  pseudopopsp <- sjlong3[pseudopopsp, on=c("isoform","population" )]
  # keep population specific junctions
  pseudopopsp[, sj_total_ocurrences:=.N, by="junction"]
  pseudopopsp[, sj_pop_ocurrences:=.N, by=c("junction", "population")]
  pseudopopsp[, sj_pop_specific:=ifelse(sj_total_ocurrences==sj_pop_ocurrences, "Pop Specific", "Pop Shared")]
  pseudopopsp[sj_pop_specific=="Pop Specific"]
  
  # keep population specific junctions
  sjlong3 <- sjlong3[population_sharing==1]
  pseudodt <-as.data.table(as.data.frame(table(unique(sjlong3[, .(population, popreads, sj_category, junction)])[, .(population, sj_category)])))
  finaldt <-pseudodt[, propsj:=Freq/sum(Freq), by="population"][sj_category=="Known SJ", .(population, propsj)]
  setorder(finaldt, -propsj)
  resvecprop <- c(resvecprop, which(finaldt$population=="AJI"))
  resvecpos <- c(resvecpos, finaldt[population=="AJI", .(propsj)])
  rm(sjlong3, pseudodt,pseudopopsp)
  gc()}

fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/PODER_sj_permutation.propKnownSJ_in_popSpSJ.AJI.tsv")

fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/PODER_sj_permutation.propKnownSJ_in_popSpSJ.tsv")
emppval <- sum(as.numeric(unlist(resvecpos))>0.949)/1000
emppval <- sum(as.numeric(unlist(resvecpos))>0.963)/1000
ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram(bins=30, color="black", fill="darkgrey")+
  geom_vline(xintercept=0.949, linetype="dashed", color="darkred")+
  mytheme+
  annotate(geom="text",label=paste0("Empirical      \nP-value= ", emppval), x=0.949, y=75, hjust=-0.1)+
  labs(x="Proportion of Known SJ in pseudoPopulation-Specific SJ", y="Permutations Count")
ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram(bins=30, color="black", fill="darkgrey")+
  geom_vline(xintercept=0.963, linetype="dashed", color="darkred")+
  mytheme+
  annotate(geom="text",label=paste0("Empirical      \nP-value= ", emppval), x=0.963, y=75, hjust=1.1)+
  labs(x="Proportion of Known SJ in pseudoPopulation-Specific SJ", y="Permutations Count")
ggsave("../../10_figures/suppfig/barplot.permutation_knownSJpopspecific.pdf", dpi=700, width = 14, height = 11,  units = "cm")





##############################################################--------------



pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(sjone)[grepl(pattern, colnames(sjone))][]
sjonlong <- melt(sjone, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]

sjonlong <- unique(metadata[,.(sample, map_reads_assemblymap)])[sjonlong, on="sample"]
sjonlong[, population:=gsub(".$", "", sample)]
sjonlong[, population:=factor(population, levels=names(popcol))]
popdata <-unique(unique(sjonlong[,.(population, map_reads_assemblymap)])[, popreads :=sum(map_reads_assemblymap), by="population"][, individuals_per_pop:=.N, by="population"][, map_reads_assemblymap:=NULL])
sjonlong <- popdata[sjonlong, on="population"]
sjonlong[, SJ_count:= .N, by=c("population", "sj_category")]
sjonlong[, norm_SJ_count:= .N/popreads/individuals_per_pop, by=c("population", "sj_category")]

sjonlong[, sj_category:=factor(sj_category, level=c("knownSJ", "novelSJ_knownSS_2", "novelSJ_knownSS_1", "novelSJ_knownSS_0"))]


pseudodt1 <-as.data.table(as.data.frame(table(unique(sjonlong[, .(population, popreads, sj_category, junction)])[, .(population, sj_category)])))
finaldt1 <-pseudodt1[, propsj:=Freq/sum(Freq), by="population"][sj_category=="knownSJ", .(population, propsj)]


ggplot(unique(sjonlong[, .(population, SJ_count, sj_category, popreads)]), aes(x=reorder(population, popreads), y=SJ_count, fill=sj_category))+
  geom_col(position="dodge")+
  scale_fill_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  mytheme+
  geom_text(aes(label = SJ_count), stat="identity", position=position_dodge(width=0.9),vjust=-0.5)+
  labs(x="", y="SJ count")
ggplot(unique(sjonlong[, .(population, SJ_count, sj_category, popreads)]), aes(x=reorder(population, popreads), y=SJ_count, fill=sj_category))+
  geom_col(position="fill")+
  scale_fill_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  mytheme+
  labs(x="", y="SJ count")
ggplot(unique(sjonlong[, .(population, norm_SJ_count, sj_category, popreads)]), aes(x=reorder(population, popreads), y=norm_SJ_count, fill=sj_category))+
  geom_col(position="dodge")+
  scale_fill_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  mytheme+
  labs(x="", y="SJ count")


# # PERMUTATION TEST -------I'll test proportions, are CEU having a largest proportion of knownSJ in population specific transcripts??
# pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
# pattern <- paste(pops, collapse = "|")
# samplecols <- colnames(sj)[grepl(pattern, colnames(sj))]
# 
# sj2 <- sj[sample_sharing>1][, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]
# 
# sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]
# samplingvec <- gsub(".$", "", unique(sjlong$sample))
# 
# resvecpos <- c()
# resvecprop <- c()
# for(iteration in 1:1000){
#   randompop <-sample(samplingvec, 43)
#   sjlong3 <- copy(sjlong)
#   names(randompop) <- unique(sjlong3$sample)
#   
#   sjlong3[, `:=`(population=randompop[sample])]
#   sjlong3 <- unique(metadata[, .(sample, map_reads_assemblymap)])[sjlong3, on="sample"]
#   sjlong3[, popreads := sum(unique(map_reads_assemblymap)), by="population"]
#   
#   # detect pseupop sp transcripts
#   pseudopopsp <- unique(unique(sjlong3[, .(isoform, population, sample)])[, samples_x_pseudopop:=.N, by=c("isoform", "population")][samples_x_pseudopop>1][, .(isoform, population, samples_x_pseudopop)])
#   pseudopopsp <- pseudopopsp[, popsharing := .N, by="isoform"][popsharing==1]
#   # keep population specific junctions
#   sjlong3 <- pseudopopsp[sjlong3, on=c("isoform","population" )]
#   sjlong3 <- sjlong3[popsharing==1]
#   pseudodt <-as.data.table(as.data.frame(table(unique(sjlong3[, .(population, popreads, sj_category, junction)])[, .(population, sj_category)])))
#   finaldt <-pseudodt[, propsj:=Freq/sum(Freq), by="population"][sj_category=="knownSJ", .(population, propsj)]
#   setorder(finaldt, -propsj)
#   resvecprop <- c(resvecprop, which(finaldt$population=="CEU"))
#   resvecpos <- c(resvecpos, finaldt[population=="CEU", .(propsj)])
#   rm(sjlong3, pseudodt,pseudopopsp)
#   gc()}
# 
# fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_permutationSJinpopspTranscripts.tsv")
# 
# 
# 


