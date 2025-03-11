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
# fwrite(reslongmeta, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_withoutchrEBV_reslongmeta.tsv")
# reslongmeta <- fread("04_evaluation/05_mastertable/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta.tsv")
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
sqanti <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_sqanti_output/240909merge_classification.txt") # sqanti table
sqanti <- sqanti[, .(isoform, length, exons, structural_category, associated_gene, associated_transcript, ref_length, ref_exons,  subcategory, all_canonical)]
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot.tsv")
annot <- sqanti[annot, on=c("isoform"="transcript")]
annot <- annot[contig!="chrEBV"]
#save checkpoint
fwrite(annot, "04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti.tsv")

##### JUNCTIONS ------------------------------------------------------------------------------------------------
# load data
myjunctions <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_sqanti_output/240909merge_junctions.txt") # sqanti sj table
recount <- fread("04_evaluation/03_recount3/data/recount3_srv3h_morethan10counts.tsv") # recount 3 with more than 10 counts per SJ
colnames(recount) <- c("contig", "start", "end", "strand", "novel", "donor", "acceptor", "recount_samples", "recount_counts", "junction")
correspondence <- fread("04_evaluation/01_adapt_lyric_output/data/ref_contigs_correspondence.tsv", header=F)
# parse my junctions
myjunctions <-myjunctions[, .(isoform, chrom, strand, junction_number, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, canonical)]
myjunctions[, junction := paste0(chrom, ":", genomic_start_coord, ":", genomic_end_coord, ":", strand)]

correspondencevec <- correspondence$V1
names(correspondencevec) <- correspondence$V2
recount[, contig:=correspondencevec[contig]][, junction:=paste(contig, start,end,strand, sep=":")]
mergedsj <-recount[, .(junction, recount_samples, recount_counts)][myjunctions, on="junction"]
rm(recount, myjunctions, correspondence)
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
fwrite(mergedsj, "04_evaluation/05_mastertable/data/240909merge_sj_table_updatedrecount.tsv", quote = F, row.names = F, sep="\t")


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

fwrite(sj_trxwise, "04_evaluation/05_mastertable/data/240909merge_sj_trxwise_updatedrecount.tsv")
rm(mergedsj, newmerged, novelsj)
gc()
# merge with annotation
sj_trxwise <-fread("04_evaluation/05_mastertable/data/240909merge_sj_trxwise_updatedrecount.tsv")
annot <-fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti.tsv")
annot <- sj_trxwise[annot, on=c("isoform")]
annot[, `:=`(platform=NULL, feature=NULL, unk1=NULL, unk2=NULL, info=NULL, tool_sample_pairs=NULL, tool_sample_pairs_split=NULL)]
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_updatedrecount.tsv")

# add information about the associated genes
annot <- fread("04_evaluation/02_sqanti/data/240909merge_withoutchrEBV_mastertable/240909merge_reslongmeta_annot_sqanti_sj_updatedrecount.tsv")
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
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_updatedrecount.tsv")

# read
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_updatedrecount.tsv")
quanti <- fread("../06_quantification/02_flairquantify/data/240917_merge_geneEntry_correctedScaffolds_nochrEBV/quantification_stats_stringent_for_master_table.tsv")

# add quantification data to master table
annot <- quanti[annot, on="isoform"]
# save
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_updatedrecount.tsv")


# add info about novel loci
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_updatedrecount.tsv")
transcriptlvl <- fread("../../novelannotations/merged/240917_merge_geneEntry_correctedScaffolds_nochrEBV_onlyTranscriptFeatures.gtf")
transcriptlvl[,geneid_afterbuildloci := tstrsplit(V9, "\"")[[6]]]
transcriptlvl[,isoform := tstrsplit(V9, "\"")[[2]]]
annot <- unique(transcriptlvl[, .(isoform, geneid_afterbuildloci)])[annot, on="isoform"]
annot[, associated_geneid.v := ifelse(structural_category=="intergenic", geneid_afterbuildloci, associated_geneid.v)]
annot[, geneid_afterbuildloci:=NULL]
annot[, discovered_transcripts_per_gene:=.N, by=associated_geneid.v]
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_updatedrecount.tsv")

#### ADD protein prediciton info

protein <- fread("../../novelannotations/protein/240919_protein_annotation.tsv")
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_updatedrecount.tsv")

protein <- protein[, .(tid, blastp_identity, blastp_bitscore, transcript_length_nt, orf_length_nt, protein_length_cd, protein_is_nmd, protein_splice_category, protein_splice_subcategory, protein_sequence)]
protein[, `:=`(isoform=tid, tid=NULL)]

annot <- protein[annot, on="isoform"]
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount.tsv")

# disambiguoate geneids
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount.tsv")
gencode <- fread("/home/pclavell/mounts/mn5/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
gencode <- gencode[, .(geneid.v, transcriptid.v, gene_biotype, gene_name, transcript_biotype)]
gencodevec <- gencode$geneid.v
names(gencodevec) <- gencode$transcriptid.v

biotypevec <- gencode$gene_biotype
names(biotypevec) <- gencode$geneid.v
annot[associated_gene_biotype=="" & structural_category%in%c("full-splice_match", "incomplete-splice_match")]
annot[, old_associated_geneid.v:=associated_geneid.v]
annot[, old_associated_gene_biotype:=associated_gene_biotype]
annot[, associated_geneid.v:=ifelse(associated_gene_biotype=="" & structural_category%in%c("full-splice_match", "incomplete-splice_match"), gencodevec[associated_transcriptid.v], associated_geneid.v)]
annot[, associated_gene_biotype:=ifelse(associated_gene_biotype=="" & structural_category%in%c("full-splice_match", "incomplete-splice_match"), biotypevec[associated_geneid.v], associated_gene_biotype)]
fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes.tsv")


# Compute some extra metrics and remove some NA
annot <- fread("04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes.tsv")

# How many transcripts are there associated to each gene by sqanti category
geneperstrcount <-unique(annot[, .(trxPerGeneAndCategory = .N), by=c("associated_geneid.v", "structural_category")])
transcriptperstrcount <-unique(annot[, .(trxPerAssTrxAndCategory = .N), by=c("associated_transcriptid.v", "structural_category")])


geneperstrcount_wide <- dcast(unique(geneperstrcount[, .(associated_geneid.v, structural_category, trxPerGeneAndCategory)]), associated_geneid.v~structural_category, fill=0)
transcriptperstrcount_wide <- dcast(unique(transcriptperstrcount[, .(associated_transcriptid.v, structural_category, trxPerAssTrxAndCategory)]), associated_transcriptid.v~structural_category, fill=0)
colnames(geneperstrcount_wide)[2:ncol(geneperstrcount_wide)] <- paste0("trx_pergene_count_",colnames(geneperstrcount_wide)[2:ncol(geneperstrcount_wide)] )
colnames(transcriptperstrcount_wide)[2:ncol(transcriptperstrcount_wide)] <- paste0("trx_per_asstrx_count_",colnames(transcriptperstrcount_wide)[2:ncol(transcriptperstrcount_wide)] )

annot <- geneperstrcount_wide[annot, on = .(associated_geneid.v)]
annot <- transcriptperstrcount_wide[annot, on = .(associated_transcriptid.v)]

annot[, existsFSMinGene:=ifelse(`trx_pergene_count_full-splice_match`>0, TRUE, FALSE)]
annot[, existsFSMinTranscript:=ifelse(`trx_per_asstrx_count_full-splice_match`>0, TRUE, FALSE)]

# replace NA by 0
annot$flair_max_counts <- replace_na(annot$flair_max_counts, 0)
annot$flair_min_counts <- replace_na(annot$flair_min_counts, 0)
annot$flair_mean_counts <- replace_na(annot$flair_mean_counts, 0)
annot$flair_expressed_samples <- replace_na(annot$flair_expressed_samples, 0)
annot$flair_total_counts <- replace_na(annot$flair_total_counts, 0)

# relabel empty gene biotypes
annot[, associated_gene_biotype:=ifelse(associated_gene_biotype=="", "novel/ambiguous gene")]

# change categories names
categories <- c("FSM", "ISM", "NIC", "Intergenic", "NNC", "Fusion", "Antisense", "Genic")
names(categories) <- unique(annot$structural_category)
annot$structural_category <- categories[annot$structural_category]

fwrite(annot, "04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo.tsv")


# 
# # is flair quantification biasing the results?
# 
# ggplot(annot, aes(x=factor(flair), y=log10(flair_total_counts+1)))+
#   geom_violin()+
#   geom_boxplot(outliers=F, width=0.15)+
#   mytheme
# # what percentage of transcripts discovered by each tool are found in more tools
# toollong <-melt(annot[, .(flair, espresso, isoquant, lyric, tool_sharing, flair_mean_counts, flair_total_counts)], measure.vars=c("flair", "espresso", "isoquant", "lyric"), variable.name = "tool", value.name="detected")
# toollong <- toollong[detected==1]
# ggplot(toollong, aes(x=tool, fill=factor(tool_sharing)))+
#   geom_bar()+
#   mytheme+
#   labs(x="", y="# Detected transcripts", fill="Tool sharing")
# 
# 
# # which tool does have a higher percentage of stuff passing expresion filter
# toollong[, flair_mean_counts:=ifelse(is.na(flair_mean_counts), 0, flair_mean_counts)]
# toollong[, overfilt:=ifelse(flair_mean_counts>=1.3, "kept transcript", "dropped transcript")]
# ggplot(toollong, aes(x=tool, fill=factor(overfilt)))+
#   geom_bar()+
#   mytheme+
#   labs(x="", y="# Detected transcripts", fill="Mean counts>=1.3")+
#   facet_wrap(~factor(tool_sharing))+
#   scale_fill_manual(values=c("darkred", "#2D93AD"))
# ggplot(toollong, aes(x=tool, fill=factor(overfilt)))+
#   geom_bar(position="fill")+
#   mytheme+
#   labs(x="", y="# Detected transcripts", fill="Mean counts>=1.3")+
#   facet_wrap(~factor(tool_sharing))+
#   scale_fill_manual(values=c("darkred", "#2D93AD"))
# 
# ggplot(toollong, aes(x=log10(flair_total_counts+1), fill=factor(tool_sharing)))+
#   geom_bar()+
#   facet_wrap(~factor(tool_sharing))+
#   mytheme+
#   ylab("Tool sharing")+
#   guides(fill="none")
