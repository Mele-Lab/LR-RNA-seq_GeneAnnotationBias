## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES: Script to create a master table of a gff
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly"
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

################THRESHOLDS FOR RECOUNT




##
## 0----------------------------END OF HEADER----------------------------------0

#### LOAD DATA
metadata <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv") # metadata
annot <- fread("../../novelannotations/merged/240909_merged_withoutchrEBV.gtf") # merged annotation

#### PARSE ANNOTATION
colnames(annot) <- c("contig", "platform", "feature", "start", "end", "unk1", "strand", "unk2", "info")

annot <- annot[feature=="transcript",]
gc()
annot[, tool_sample_pairs := tstrsplit(info, "\"")[[4]]]
annot[, transcript := tstrsplit(info, "\"")[[2]]]

# Split the tool_sample_pairs into individual pairs
annot[, tool_sample_pairs_split := strsplit(tool_sample_pairs, ",")]
# Transform the annotation to a long format to be able to compute sharing across samples, tools and populations
all_tool_sample_pairs <- unique(unlist(annot$tool_sample_pairs_split))
gc()
# Create a data.table with the transcripts and their corresponding tool_sample_pairs
annot_expanded <- annot[, .(transcript = transcript, tool_sample_pair = unlist(tool_sample_pairs_split)), by = 1:nrow(annot)]
gc()
# Create the result matrix directly using sparse matrix construction
rows<-nrow(annot)
cols <-length(all_tool_sample_pairs)
mylist<-list(annot$transcript, all_tool_sample_pairs)
gc()
result_matrix <- Matrix::sparseMatrix(
  i = match(annot_expanded$transcript, annot$transcript),
  j = match(annot_expanded$tool_sample_pair, all_tool_sample_pairs),
  x = 1,
  dims = c(rows, cols),
  dimnames = mylist
)

# Convert the sparse matrix to a data frame
result_matrix <- as.data.frame(as.matrix(result_matrix))

# melt
reslong <- melt(rownames_to_column(result_matrix), variable.name="sample_tool_pair", value.name="detected")
rm(result_matrix, mylist, rows, cols, all_tool_sample_pairs, annot_expanded)
gc()
setDT(reslong)
reslong[, transcript:=rowname]
reslong[, tool:=gsub("_.*","",sample_tool_pair)]
reslong[, cell_line_id:=gsub(".*_","",gsub("_1", "", sample_tool_pair))]
gc()
# add metadata
reslongmeta <- metadata[, .(sample, cell_line_id, population)][reslong, on="cell_line_id"]
rm(reslong)
gc()
reslongmeta <- reslongmeta[detected>0,][, `:=`(rowname=NULL, detected=NULL, cell_line_id=NULL)]
gc()
reslongmeta$tool <- ifelse(reslongmeta$tool=="iq", "isoquant", reslongmeta$tool)
#save checkpoint
fwrite(reslongmeta, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_withoutchrEBV_reslongmeta.tsv")
reslongmeta <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta.tsv")
# compute sharings
reslongmeta[, detected:=1]
tool_sharing_temp <- unique(reslongmeta[, .(transcript, tool,detected)])
tool_sharing <- dcast(tool_sharing_temp, transcript~tool, fill=0)[, tool_sharing:=espresso+flair+isoquant+lyric]
population_sharing <-unique(unique(reslongmeta[, .(transcript, population)])[ ,population_sharing:=.N, by=transcript][, population:=NULL])
sample_sharing_temp <- unique(reslongmeta[, .(transcript, sample, detected)])
sample_sharing <- dcast(sample_sharing_temp, transcript~sample, fill=0)
sample_sharing[,sample_sharing := rowSums(.SD), .SDcols = colnames(sample_sharing)[2:ncol(sample_sharing)]]
samples_x_pop_sharing <-unique(unique(reslongmeta[, .(transcript, sample, population, detected)])[ ,samples_x_pop_sharing:=.N, by=c("transcript", "population")][, sample:=NULL])
samples_x_pop_sharing_wide <- dcast(samples_x_pop_sharing, transcript~population, fill=0)
gc()


# merge  to annot
annot <- tool_sharing[annot, on="transcript"]
annot <- population_sharing[annot, on="transcript"]
annot <- sample_sharing[annot, on="transcript"]
annot <- samples_x_pop_sharing_wide[annot, on="transcript"]
rm(population_sharing, sample_sharing, sample_sharing_temp, samples_x_pop_sharing, samples_x_pop_sharing_wide, tool_sharing, tool_sharing_temp)
gc()




#save checkpoint
fwrite(annot, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot.tsv")

##### SQANTI -------------------------------------------------------------------------------------------------
# LOAD MORE DATA
sqanti <- fread("04_evaluation/02_sqanti/data/240909merge/240909merge_classification.txt") # sqanti table
sqanti <- sqanti[, .(isoform, length, exons, structural_category, associated_gene, associated_transcript, ref_length, ref_exons,  subcategory, all_canonical)]
annot <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot.tsv")
annot <- sqanti[annot, on=c("isoform"="transcript")]
annot <- annot[contig!="chrEBV"]
#save checkpoint
fwrite(annot, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti.tsv")

##### JUNCTIONS ------------------------------------------------------------------------------------------------
# load data
myjunctions <- fread("04_evaluation/02_sqanti/data/240909merge/240909merge_junctions.txt") # sqanti sj table
recount <- fread("04_evaluation/03_recount3/data/recount3_srv3h_morethan10counts.tsv") # recount 3 with more than 10 counts per SJ
colnames(recount) <- c("contig", "start", "end", "strand", "novel", "donor", "acceptor", "recount_samples", "recount_counts", "junction")

# parse my junctions
myjunctions <-myjunctions[, .(isoform, chrom, strand, junction_number, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, canonical)]
myjunctions[, junction := paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand)]

mergedsj <-recount[, .(junction, recount_samples, recount_counts)][myjunctions, on="junction"]
rm(recount, myjunctions)
gc()
mergedsj$recount_counts <- replace_na(mergedsj$recount_counts, 0)
mergedsj$recount_samples <- replace_na(mergedsj$recount_samples, 0)

mergedsj[, .(recount_samples, recount_counts, isoform, junction, junction_category, start_site_category, end_site_category, canonical)]


mergedsj[, sj_category:=ifelse(junction_category=="novel", 
                               ifelse(rowSums(cbind(start_site_category=="known", end_site_category=="known"))>0, 
                                      ifelse(rowSums(cbind(start_site_category=="known", end_site_category=="known"))==1,"novelSJ_knownSS_1", 
                                             "novelSJ_knownSS_2"), 
                                      "novelSJ_knownSS_0"),
                               "knownSJ")]


mergedsj[, sj_category_full := paste(sj_category,canonical, sep="_")]


# Convert SJ data.table to transcript-wise SJ info
mergedsj[, sj_recount_counts := paste(recount_counts, collapse=","), by=isoform]
mergedsj[, sj_category_full := paste(sj_category_full, collapse=","), by=isoform]
mergedsj[, sj_less_recountsupported_counts:=min(recount_counts), by=isoform]
mergedsj[, sj_less_recount_counts_info:=sj_category[which.min(recount_counts)], by=isoform]
mergedsj[, sj_less_recountsupported_samples:=min(recount_samples), by=isoform]
mergedsj[, sj_less_recount_samples_info:=sj_category[which.min(recount_samples)], by=isoform]


novelsj <-mergedsj[junction_category=="novel",][, `:=`(sj_less_recountsupported_novel_counts= min(recount_counts),
                                                       sj_less_recount_counts_novel_info=sj_category[which.min(recount_counts)],
                                                       sj_less_recountsupported_novel_samples= min(recount_samples),
                                                       sj_less_recount_samples_novel_info=sj_category[which.min(recount_samples)]), by="isoform"][, .(isoform, sj_less_recountsupported_novel_counts, sj_less_recount_counts_novel_info,sj_less_recountsupported_novel_samples,sj_less_recount_samples_novel_info)]
newmerged <- novelsj[mergedsj, on="isoform" ]
# create the datatable to merge with the sqanti output
mylogical <- c(colnames(newmerged)=="isoform" | grepl("sj_", colnames(newmerged)))
sj_trxwise <- unique(newmerged[,..mylogical][, sj_category:=NULL])

fwrite(sj_trxwise, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_sj_trxwise.tsv")
# merge with annotation
sj_trxwise <-fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_sj_trxwise.tsv")
annot <-fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti.tsv")
annot <- sj_trxwise[annot, on=c("isoform")]
annot[, `:=`(platform=NULL, feature=NULL, unk1=NULL, unk2=NULL, info=NULL, tool_sample_pairs=NULL, tool_sample_pairs_split=NULL)]
fwrite(annot, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti_sj.tsv")

# add information about the associated genes
annot <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti_sj.tsv")
gencode <- fread("/home/pclavell/mounts/mn5/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
gencode <- gencode[, .(geneid.v, transcriptid.v, gene_biotype, gene_name, transcript_biotype)]

# count the number of transcripts per annotated gene
gencode[, transcripts_per_gene:=.N, by=geneid.v]

annot <- unique(gencode[, .(geneid.v, gene_biotype, gene_name, transcripts_per_gene)])[annot, on=c("geneid.v"="associated_gene")]
annot <- unique(gencode[, .(transcriptid.v, transcript_biotype)])[annot, on=c("transcriptid.v"="associated_transcript")]

replacenames <- c("transcriptid.v", "transcript_biotype", "geneid.v", "gene_biotype", "gene_name", "transcripts_per_gene")
newnames <- paste0("associated_", replacenames)
names(newnames) <- replacenames
colnames(annot)[colnames(annot)%in%replacenames] <- newnames[colnames(annot)[colnames(annot)%in%replacenames]]
# count number of discovered transcripts per gene
annot[, discovered_transcripts_per_gene:=.N, by=associated_geneid.v]
# save
fwrite(annot, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti_sj_gencodev47.tsv")
#annot <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti_sj_gencodev47.tsv")

# # create table relating transcripts and geneid
# fwrite(unique(annot[, .(associated_geneid.v, isoform)]), "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_transcript_associatedgene_correspondence.tsv", col.names = F, quote = F, sep="\t")
# fwrite(unique(annot[structural_category=="intergenic", .(associated_geneid.v, isoform)]), "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_novelgene_transcripts.tsv", col.names = F, quote = F, sep="\t")
# test <-unique(annot[structural_category=="intergenic", .(associated_geneid.v, isoform, contig)])
# table(test$contig)              


test <- fread("../02_ONT_preprocessing/arrayq10_fastqgz", header = F)
test[, lab_sampleid:=tstrsplit(V1,"_")[[4]]]
myfile <-metadata[, .(lab_sampleid, sample)][test, on="lab_sampleid"]
fwrite(myfile[, `:=`(condition="condition1", batch="batch1")][, .(sample, condition, batch, V1)],
       "../06_quantification/01_isoquantify/array_flair", row.names = F, col.names = F, quote = F, sep="\t")
