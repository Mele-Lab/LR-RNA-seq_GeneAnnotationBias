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



# PERMUTATION TEST -------I'll test proportions, are CEU having a largest proportion of knownSJ in population specific transcripts??

# Identify columns showing in which samples a transcript is found
pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(sj)[grepl(pattern, colnames(sj))]

# Keep trx found in at least 2 samples
sj2 <- sj[sample_sharing>1][, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]

sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")[detected!=0]
samplingvec <- gsub(".$", "", unique(sjlong$sample))



## Find proportions of SJ
sjlong[, realpop := gsub(".$", "", sample)]

# Deduplicate and filter directly
sjlong3 <- unique(sjlong[, .(sample, sj_category, junction, realpop)])
pseudopopsp <- sjlong3[, .N, by = .(junction, realpop)][N > 1, .(junction)]












# 
# start_time <- Sys.time()
# 
# resvecprop <- numeric(200)  # Pre-allocate vector for efficiency
# 
# setDT(sjlong)  # Ensure sjlong is a data.table
# samplingvec_len <- length(samplingvec)
# 
# for(iteration in 1:200){
#   # Efficient sampling and mapping population names using vectorized operations
#   randompop <- setNames(sample(samplingvec, 43, replace = FALSE), unique(sjlong$sample))
#   
#   # Modify population directly
#   sjlong[, population := randompop[sample]]
#   
#   # Deduplicate and filter directly
#   sjlong3 <- unique(sjlong[, .(sample, sj_category, junction, population)])
#   pseudopopsp <- sjlong3[, .N, by = .(junction, population)][N > 1, .(junction)]
#   
#   sjlong3 <- sjlong3[junction %chin% pseudopopsp$junction]  # Use %chin% for faster lookup
#   
#   # Compute frequency directly
#   finaldt <- sjlong3[, .(propsj = sum(sj_category == "knownSJ") / .N), by = population]
#   
#   # Direct assignment
#   resvecprop[iteration] <- finaldt[population == "CEU", propsj]
#   
#   # Memory cleanup
#   rm(randompop, sjlong3, pseudopopsp, finaldt)
#   gc()
# }
# 
# end_time <- Sys.time()
# time_taken <- end_time - start_time
# print(time_taken)

# fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_permutationSJinpopspTranscripts.tsv")


ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram()+
  geom_vline(xintercept=0.902)+
  mytheme+
  annotate(geom="text",label="p=0.04", x=0.902, y=75)+
  labs(x="Proportion of known SJ", y="Permutations count")



