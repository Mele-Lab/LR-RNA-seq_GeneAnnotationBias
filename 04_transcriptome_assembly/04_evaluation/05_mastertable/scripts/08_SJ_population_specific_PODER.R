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
data <- fread("05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv")
sj <- fread("05_mastertable/data/240909merge_sj_table.tsv")
#sj <- fread("02_sqanti/data/240909merge_withoutchrEBV_sqanti_output/240909merge_junctions.txt")
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
popcol <- unique(metadata$color_pop)
names(popcol) <- unique(metadata$population)

# merge
sj <- data[filter=="pass",sj_category_full_trx:=sj_category_full][,sj_category_full:=NULL][sj, on="isoform"]

# keep only sj belonging to population specific transcripts
sjone <- sj[population_sharing==1 & sample_sharing>=2,]
sjonelong <- melt(, measure.vars = c("AJI", "CEU", "MPC", "YRI", "LWK", "HAC", "ITU", "PEL"), value.name = "detected", variable.name = "population")[detected>=2]
#### PARENTHESIS TO TRY TO IDENTIFY POP SPECIFIC JUNCTIONS------------------


sjonelong[, sj_total_ocurrences:=.N, by="junction"]
sjonelong[, sj_pop_ocurrences:=.N, by=c("junction", "population")]
sjonelong[, sj_pop_specific:=ifelse(sj_total_ocurrences==sj_pop_ocurrences, "Pop Specific", "Pop Shared")]

popspsj <-sjonelong[sj_pop_specific=="Pop Specific"]
popspsj <- unique(metadata[, .(population, map_reads_assemblymap)][, .(popdepth=sum(map_reads_assemblymap)), by="population"])[popspsj, on="population"]
popspsj[, popspsj_x_pop_x_cat:=.N, by=c("population", "sj_category")]
popspsj[,  eurpop:= ifelse(population%in%c("CEU", "AJI"), "EUR", "nonEUR")]

ggplot(unique(popspsj[, .(popdepth, popspsj_x_pop_x_cat,sj_category, eurpop)]), aes(x=popdepth, y=popspsj_x_pop_x_cat, col=sj_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  labs(x="Total population mapped reads", y="# Population Specific Splice Junctions")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c("red", "grey"))+
  scale_size_manual(values=c(5,3))+
  scale_alpha_manual(values=c(0.6, 0.9))

ggplot(unique(popspsj[eurpop=="nonEUR", .(popdepth, popspsj_x_pop_x_cat,sj_category, eurpop)]), aes(x=popdepth, y=popspsj_x_pop_x_cat, col=sj_category))+
  geom_line(linewidth=1.5)+
  mytheme+
  scale_color_manual(values=c("darkgreen","deepskyblue3", "darkgoldenrod3", "darkred"))+
  labs(x="Total population mapped reads", y="# Population Specific Splice Junctions")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eurpop, size=eurpop, alpha=eurpop))+
  scale_color_manual(values=c( "grey"))+
  scale_size_manual(values=c(3))+
  scale_alpha_manual(values=c( 0.9))

# PERMUTATION TEST -------I'll test proportions, are CEU having a largest proportion of knownSJ within their population specific SJ ??
pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(sj)[grepl(pattern, colnames(sj))]

# keep only transcripts that are shared at least across 2 samples because (if not, they won't pass the popsp filter)
sj2 <- sj[sample_sharing>1 & filter=="pass"][, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]

sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]
samplingvec <- gsub(".$", "", unique(sjlong$sample))

resvecpos <- c()
resvecprop <- c()
for(iteration in 1:1000){
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
  
  pseudopopsp[, sj_total_ocurrences:=.N, by="junction"]
  pseudopopsp[, sj_pop_ocurrences:=.N, by=c("junction", "population")]
  pseudopopsp[, sj_pop_specific:=ifelse(sj_total_ocurrences==sj_pop_ocurrences, "Pop Specific", "Pop Shared")]
  pseudopopsp[sj_pop_specific=="Pop Specific"]
  
  # keep population specific junctions
  sjlong3 <- sjlong3[popsharing==1]
  pseudodt <-as.data.table(as.data.frame(table(unique(sjlong3[, .(population, popreads, sj_category, junction)])[, .(population, sj_category)])))
  finaldt <-pseudodt[, propsj:=Freq/sum(Freq), by="population"][sj_category=="knownSJ", .(population, propsj)]
  setorder(finaldt, -propsj)
  resvecprop <- c(resvecprop, which(finaldt$population=="CEU"))
  resvecpos <- c(resvecpos, finaldt[population=="CEU", .(propsj)])
  rm(sjlong3, pseudodt,pseudopopsp)
  gc()}





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


# PERMUTATION TEST -------I'll test proportions, are CEU having a largest proportion of knownSJ in population specific transcripts??
pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(sj)[grepl(pattern, colnames(sj))]

sj2 <- sj[sample_sharing>1][, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]

sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]
samplingvec <- gsub(".$", "", unique(sjlong$sample))

resvecpos <- c()
resvecprop <- c()
for(iteration in 1:1000){
  randompop <-sample(samplingvec, 43)
  sjlong3 <- copy(sjlong)
  names(randompop) <- unique(sjlong3$sample)
  
  sjlong3[, `:=`(population=randompop[sample])]
  sjlong3 <- unique(metadata[, .(sample, map_reads_assemblymap)])[sjlong3, on="sample"]
  sjlong3[, popreads := sum(unique(map_reads_assemblymap)), by="population"]
  
  # detect pseupop sp transcripts
  pseudopopsp <- unique(unique(sjlong3[, .(isoform, population, sample)])[, samples_x_pseudopop:=.N, by=c("isoform", "population")][samples_x_pseudopop>1][, .(isoform, population, samples_x_pseudopop)])
  pseudopopsp <- pseudopopsp[, popsharing := .N, by="isoform"][popsharing==1]
  # keep population specific junctions
  sjlong3 <- pseudopopsp[sjlong3, on=c("isoform","population" )]
  sjlong3 <- sjlong3[popsharing==1]
  pseudodt <-as.data.table(as.data.frame(table(unique(sjlong3[, .(population, popreads, sj_category, junction)])[, .(population, sj_category)])))
  finaldt <-pseudodt[, propsj:=Freq/sum(Freq), by="population"][sj_category=="knownSJ", .(population, propsj)]
  setorder(finaldt, -propsj)
  resvecprop <- c(resvecprop, which(finaldt$population=="CEU"))
  resvecpos <- c(resvecpos, finaldt[population=="CEU", .(propsj)])
  rm(sjlong3, pseudodt,pseudopopsp)
  gc()}

fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_permutationSJinpopspTranscripts.tsv")


ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram()+
  geom_vline(xintercept=0.902)+
  mytheme+
  annotate(geom="text",label="p=0.04", x=0.902, y=75)+
  labs(x="Proportion of known SJ", y="Permutations count")

ggplot(data.frame("permutations"=as.numeric(unlist(resvecprop))), aes(x=permutations))+
  geom_histogram()+
  geom_vline(xintercept=0.902)+
  mytheme+
  annotate(geom="text",label="p=0.04", x=0.902, y=75)+
  labs(x="Proportion of known SJ", y="Permutations count")
#############


