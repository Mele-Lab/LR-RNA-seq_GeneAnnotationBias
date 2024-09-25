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

MINTHRESHOLD=50

##
## 0----------------------------END OF HEADER----------------------------------0

# load metadata
metadata <- fread("../00_metadata/pantranscriptome_samples_metadata.tsv")


# # load gtf
# rawgtf <- fread("../../novelannotations/espresso/1_PY1_GM10492_espresso.gtf", skip = 1)
# colnames(rawgtf) <- c("contig", "novelty", "feature", "start", "end", "v6", "strand", "v8", "info")
# 
# # keep only transcripts and parse info column
# gtf <- rawgtf[feature=="transcript",]
# 
# mysplit <- strsplit(gtf$info, ";")
# transcriptidvs <- lapply(mysplit, \(x) gsub("\"", "", gsub("transcript_id \"", "", x[[1]])))
# geneidvs <- lapply(mysplit, \(x) gsub("\"", "", gsub("gene_id \"", "", x[[2]])))
# # rename columns
# gtf <- gtf[, c("transcriptid.v", "geneid.v", "info", "v6", "v8") := list(transcriptidvs, geneidvs, NULL, NULL, NULL)][]
# 
# # load sj info and recount data
# sjinfo <- fread("04_evaluation/02_sqanti/data/pseudomasked_genomic_espresso/pseudomasked_genomic_espresso_junctions.txt")
# recount <- fread("04_evaluation/03_recount3/data/recount3_srv3h_morethan10counts.tsv")
# colnames(recount) <- c("contig", "start", "end", "strand", "novel", "donor", "acceptor", "recount_samples", "recount_counts", "junction")
# 
# # PARSE BOTH DATAFRAMEs
# sjinfo <-sjinfo[, .(isoform, chrom, strand, junction_number, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, canonical)]
# sjinfo[, junction := paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand)]
# 
# mergedsj <-recount[, .(junction, recount_samples, recount_counts)][sjinfo, on="junction"]
# mergedsj$recount_counts <- replace_na(mergedsj$recount_counts, 0)
# mergedsj$recount_samples <- replace_na(mergedsj$recount_samples, 0)
# 
# 
# mergedsj[, recountsupported:=ifelse(recount_counts>=10, TRUE, FALSE)]
# 
# # Convert SJ data.table to transcript-wise SJ info
# mergedsj[, sj_support := paste(recount_counts, collapse=","), by=isoform]
# mergedsj[, sj_supported := paste(recountsupported, collapse=","), by=isoform]
# mergedsj[, sj_novelty := paste(junction_category, collapse=","), by=isoform]
# mergedsj[, sj_start_site_novelty := paste(start_site_category, collapse=","), by=isoform]
# mergedsj[, sj_end_site_novelty := paste(end_site_category, collapse=","), by=isoform]
# mergedsj[, sj_canonical := paste(canonical, collapse=","), by=isoform]
# mergedsj[, exons := sapply(strsplit(sj_canonical, ","), length)+1] # count SJ and add 1
# 
# # create the datatable to merge with the sqanti output
# mylogical <- c(colnames(mergedsj)%in%c("isoform", "chrom", "strand") | grepl("sj_", colnames(mergedsj)))
# sj_trxwise <- unique(mergedsj[,..mylogical])

## PLOT RECOUNT3 SUPPORT FOR MY JUNCTIONS
# ggplot(mergedsj, aes(x=junction_category, fill=recountsupported))+
#   geom_histogram(stat="count")
# 
# ggplot(mergedsj, aes(x=junction_category, y=log10(counts), fill=junction_category))+
#   geom_violin(alpha=0.5)+
#   geom_boxplot(width=0.15, outliers=F)+
#   mytheme
# ggplot(mergedsj, aes(x=junction_category, y=counts/samples, fill=junction_category))+
#   geom_boxplot(outliers=F)+
#   mytheme


# Load annotation
annot <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/novelannotations/merged/merged.gtf")
colnames(annot) <- c("contig", "platform", "feature", "start", "end", "unk1", "strand", "unk2", "info")

# parse
annot <- annot[feature=="transcript",]

annot[, tool_sample_pairs := tstrsplit(info, "\"")[[4]]]
annot[, transcript_id := tstrsplit(info, "\"")[[2]]]

# Split the tool_sample_pairs into individual pairs
annot[, tool_sample_pairs_split := strsplit(tool_sample_pairs, ",")]

# Create a vector of all unique tool_sample_pairs
all_tool_sample_pairs <- unique(unlist(annot$tool_sample_pairs_split))

# Initialize the matrix with 0s
result_matrix <- matrix(0, nrow = nrow(annot), ncol = length(all_tool_sample_pairs))
rownames(result_matrix) <- annot$transcript_id
colnames(result_matrix) <- all_tool_sample_pairs

# Fill in the matrix with 1s where appropriate
for (i in 1:nrow(annot)) {
  result_matrix[i, annot$tool_sample_pairs_split[[i]]] <- 1
}

result_matrix <- as.data.frame(result_matrix)
# # Plot the UpSet plot with sample-tools pairs overlaps
# UpSetR::upset(result_matrix,
#       sets = colnames(result_matrix),
#       order.by = "freq")

# melt

reslong <- melt(rownames_to_column(result_matrix), variable.name="sample_tool_pair", value.name="detected")
rm(result_matrix)
gc()
setDT(reslong)
reslong[, transcript:=rowname]
reslong[, tool:=tstrsplit(sample_tool_pair, "_")[[1]]]
reslong[, cell_line_id:=tstrsplit(sample_tool_pair, "_")[[2]]]

# add metadata
reslongmeta <- metadata[, .(sample, cell_line_id, sex, population)][reslong, on="cell_line_id"]

reslongmeta$tool <- ifelse(reslongmeta$tool=="iq", "isoquant", reslongmeta$tool)
rm(reslong)
gc()
# TOOL SPECIFICITY
# Step 1: Calculate counts
reslongmeta_filtered <- reslongmeta[detected == 1, ]
rm(reslongmeta)
gc()
reslongmeta_filtered[, tool := factor(tool, levels = names(sort(table(tool), decreasing = TRUE)))]

# Step 2: Plot using the reordered tool factor
ggplot(unique(reslongmeta_filtered[, c("transcript", "tool"), with=F]), aes(x = tool, fill = tool)) +
  geom_bar(stat = "count")+
  mytheme+
  scale_fill_manual(values=c("#F4AC45", "#004643","#6A0F49" ))+
  labs(y="Transcripts", x="")+
  guides(fill="none")+
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5)+
  ylim(0, max(table(unique(reslongmeta_filtered[, c("transcript", "tool"), with=F])$tool)) * 1.1)

toolupset <- dcast(unique(reslongmeta_filtered[, .(transcript, tool, detected)]), transcript~tool,fill=0)
# which percentage is tool specific
sum(rowSums(toolupset[, 2:4])==1)/nrow(toolupset)

toolupset_binary<-column_to_rownames(as.data.frame(toolupset), var="transcript")
UpSetR::upset(toolupset_binary,
              sets = colnames(toolupset_binary),
              order.by = "freq")
rm(toolupset, toolupset_binary)
gc()

#### Tool sharing
mytable <- data.frame(table(unlist(lapply(annot$tool_sample_pairs_split, 
                                          \(x) sum(c(sum(grepl("flair", x)), sum(grepl("espresso", x)), sum(grepl("iq", x)))>=1)
))))

ggplot(mytable, aes(fill=Var1, y=Freq, x=""))+
  geom_col(position="stack")+
  mytheme+
  scale_fill_manual(values=c("#FF8552", "#297373", "#39393A"))+
  labs(x="", y="Transcript Count", fill="Tool\nsharing")+
  geom_text(position = "stack", aes(label = Freq), vjust=1.5)



### POPULATION SPECIFICITY
# prepare upsetplot between tools
popupset <- dcast(unique(reslongmeta_filtered[, .(population, transcript, sample, detected)]), transcript~population)
popupset_binary <- as.data.frame(ifelse(popupset>1, 1, 0))# depending on the threshold the results change
popupset_binary$transcript <- popupset$transcript
popupset_binary<-column_to_rownames(popupset_binary, var="transcript")
popupset_binary <-popupset_binary[rowSums(popupset_binary)>0,]

library(UpSetR)
UpSetR::upset(popupset_binary,
      sets = colnames(popupset_binary),
      order.by = "freq",
      queries = list(
                     list(query = intersects, 
                          params = list("LWK", "YRI"), 
                          color = "orange", active = T),
                     list(query = intersects, 
                          params = list("LWK"), 
                          color = "orange", active = T),
                     list(query = intersects, 
                          params = list("YRI"), 
                          color = "orange", active = T),
                     list(query = intersects, 
                          params = list("MPC"), 
                          color = "orange", active = T)))
rm(popupset, popupset_binary)
gc()


ggplot(data.frame(table(rowSums(popupset_binary))), aes(x=as.numeric(Var1), y=Freq))+
  mytheme+
  geom_line() +
  geom_point()+
  labs(x="# populations sharing", y="# transcripts")

# SAMPLE SHARING
# Step 1: Calculate counts
reslongmeta_filtered[, cell_line_id := factor(cell_line_id, levels = names(sort(table(cell_line_id), decreasing = TRUE)))]

colorpop <-unique(metadata$color_pop)
names(colorpop) <- unique(metadata$population)

# Step 2: Plot using the reordered  factor
ggplot(unique(reslongmeta_filtered[, .(cell_line_id, transcript, population)]), aes(x = cell_line_id, fill = population)) +
  geom_bar(stat = "count")+
  mytheme+
  labs(y="Transcripts", x="")+
  guides(fill="none")+
  scale_fill_manual(values=colorpop)

cell_line_idupset <- dcast(unique(reslongmeta_filtered[, .(transcript, cell_line_id, detected)]), transcript~cell_line_id,fill=0)
# which percentage is cell_line_id specific
sum(rowSums(cell_line_idupset[, 2:4])==1)/nrow(cell_line_idupset)

cell_line_idupset_binary<-column_to_rownames(as.data.frame(cell_line_idupset), var="transcript")
UpSetR::upset(cell_line_idupset_binary,
              sets = colnames(cell_line_idupset_binary),
              order.by = "freq")
##### Sample sharing
mytable <- data.frame(table(rowSums(cell_line_idupset[,2:ncol(cell_line_idupset)])))

ggplot(mytable, aes(x=as.numeric(Var1), y=Freq))+
  labs(x="Sample Sharing", y="Transcript count (log10)")+
  mytheme+
  scale_y_log10()+
  geom_line() +
  geom_point()+
  annotation_logticks(sides="l")
ggplot(mytable, aes(x=as.numeric(Var1), y=Freq))+
  labs(x="Sample Sharing", y="Transcript count")+
  mytheme+
  geom_line() +
  geom_point()


# in how many ancestries are transcripts detected?
toolsdetection <- unique(toolupset[, .(transcript, detected, population)])[, .(populations_detected=sum(detected)), by=c( "transcript")]
ggplot(toolsdetection[, .N, by=.(populations_detected)], aes(x=populations_detected, y=N))+
  geom_point()+
  geom_line()+
  mytheme+
  guides(color="none")+
  labs(y="transcripts", x="Population sharing")
