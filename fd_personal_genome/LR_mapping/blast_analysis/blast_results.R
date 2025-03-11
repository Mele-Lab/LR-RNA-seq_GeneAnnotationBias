setwd("[PATH]/mapping_minimap/blast_analysis/results")


files.list <- "../blastAnalysis_cluster.list"
files.list <- read.table(files.list, header = F, stringsAsFactors = F)
source("~/template/pau_theme.R")
line <- files.list[1,]
all <- NULL
i <- 1
for (i in 1:nrow(files.list)){
  genome <- as.character(files.list[i,1])
  ref <- as.character(files.list[i,2])
  cat(genome, "-", ref,"\n")
  
  blast <- read.delim(paste0(genome,"_", ref, "/alignment.tsv"), header = F, stringsAsFactors = F)
  blast_bestHit <- read.delim(paste0(genome,"_", ref, "/alignment_bestHit.tsv"), header = F, stringsAsFactors = F)
  
  colnames(blast) <- c("query_id", "subject_id", "prct_identity", "alignement_length", "nb_mismatches",
                       "nb_gapOpens", "query_start", "query_end", "subject_start", "subject_end",
                       "E_value", "bit_score")
  colnames(blast_bestHit) <- c("query_id", "subject_id", "prct_identity", "alignement_length", "nb_mismatches",
                               "nb_gapOpens", "query_start", "query_end", "subject_start", "subject_end",
                               "E_value", "bit_score")
  
  #Import transcript ID
  trID_query <- scan(paste0(genome,"_", ref, "/transcript_id_PG.txt"), what = "character")
  trID_target <- scan(paste0(genome,"_", ref, "/transcript_id_",ref,".txt"), what = "character")
  
  ## TODO: Check if all the transcipts provided as inputs (target/query) are found in the alignemnt
  
  cat( "tr query: ", length(trID_query),"\n")
  cat( "tr target: ", length(trID_target),"\n")
  trID_query.nb <- length(trID_query)
  trID_target.nb <- length(trID_target)
  
  cat("tr blast found:","\n" )
  trID_query_hit.nb <- length(unique(blast$query_id))
  trID_target_hit.nb <- length(unique(blast$subject_id))
  cat( "tr query: ", trID_query_hit.nb,"\n")
  cat( "tr target: ", trID_target_hit.nb,"\n")
  
  
  length(trID_query[!(trID_query %in% blast_bestHit$query_id)])
  length(trID_target[!(trID_target %in% unique(blast$subject_id))])
  
  ## Add fake prct identity of the transcript not in the blast results
  
  blast_bestHit_red <- blast_bestHit[, c("query_id", "subject_id","prct_identity")]
  
  toAdd <- data.frame(query_id = trID_query[!(trID_query %in% blast_bestHit$query_id)], stringsAsFactors = F)
  toAdd$subject_id <- NA
  toAdd$prct_identity <- 0
  
  
  blast_bestHit_red <- rbind(blast_bestHit_red, toAdd)
  
  tmp <- NULL
  for (thr in c(seq(0,100,1), seq(99,100,0.1), seq(99.9,100,0.01))){
    counts <- sum(blast_bestHit_red$prct_identity >= thr)
    tmp <- rbind(tmp, c(thr, counts))
  }
  tmp <- as.data.frame(tmp)
  tmp <- unique(tmp)
  colnames(tmp) <- c("thr", "counts")
  
  tmp$prct  <- (tmp$counts / nrow(blast_bestHit_red))*100
  tmp$genome <- genome
  tmp$ref <- ref
  ## Sub select: keep only the one with the highest % of identity
  #res_q <- as.numeric((quantile(blast_bestHit$prct_identity, c(0,0.001,0.01,0.02, 0.03, 0.04, 0.05,0.1))))
  #res_q <- as.numeric((quantile(blast_bestHit_red$prct_identity, c(0,0.001,0.01,0.02, 0.03, 0.04, 0.05,0.1))))
  
  #toAdd <- c(genome, ref, trID_query.nb, trID_target.nb, trID_query_hit.nb, trID_target_hit.nb, res_q)
  
  all <- rbind(all, tmp)
} 

#colnames(all) <- c("genome", "ref", "trID_query.nb", "trID_target.nb", "trID_query_hit.nb", "trID_target_hit.nb",
                 #  paste0("quant_",c(0,0.001,0.01,0.02, 0.03, 0.04, 0.05,0.1)))
write.table(all, "blastSeqIdentity.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

all.hg38 <- all[all$ref == "T2T", ]
ggplot(all.hg38) +
  aes(x = thr, y = prct, colour = genome) +
  geom_line() +
  scale_color_hue(direction = 1) +
  pauTheme+
  xlab("% sequence identity") + 
  ylab("% of transripts (cumulative)")


ggplot(all.hg38) +
  aes(x = thr, y = prct, colour = genome) +
  geom_line() +
  scale_color_hue(direction = 1) +
  pauTheme+
  xlab("% sequence identity") + 
  ylab("% of transripts (cumulative)") +
  coord_cartesian(xlim = c(98,100))

# 
# 
# all <- as.data.frame(all, stringsAsFactors = F)
# 
# all[, c(3:16)] <- apply(all[, c(3:16)], 2, as.numeric)
# 
# all$trID_query_hit.prct <- (all$trID_query_hit.nb / all$trID_query.nb) * 100
# all$trID_target_hit.prct <- (all$trID_target_hit.nb / all$trID_target.nb) * 100
# 
# all.hg38 <- all[all$ref == "hg38",]
# median(all.hg38$trID_query_hit.prct)
# median(all.hg38$trID_target_hit.prct)
# median(all.hg38$quant_0.02)
# 
# numbers <- all.hg38[, -c(3:6)]
# 
# numbers <- numbers %>%
#   pivot_longer(!c(1,2), names_to = "type", values_to = "count")
# 
# ## Reordering numbers$type
# numbers$type <- factor(numbers$type,
#                        levels = c(
#                          "trID_query_hit.prct", "trID_target_hit.prct", "quant_0", "quant_0.001",
#                          "quant_0.01", "quant_0.02", "quant_0.03", "quant_0.04", "quant_0.05",
#                          "quant_0.1"
#                        )
# )
# 
# ggplot(numbers) +
#   aes(x = type, y = count, fill = type) +
#   geom_boxplot() +
#   scale_fill_manual(
#     values = c(trID_query_hit.prct = "#FF0000",
#                trID_target_hit.prct = "#D78E0C",
#                quant_0 = "#23571E",
#                quant_0.001 = "#32812A",
#                quant_0.01 = "#38A32E",
#                quant_0.02 = "#42CB35",
#                quant_0.03 = "#42D633",
#                quant_0.04 = "#38E527",
#                quant_0.05 = "#32FF1F",
#                quant_0.1 = "#2AFF15")
#   ) +
#   geom_jitter(aes(color = genome)) +
#   pauTheme     + 
#   theme(legend.position  = "none")
# 
# 


# Query ID: The identifier of the query sequence.
# Subject ID: The identifier of the subject sequence (the sequence in the database).
# % Identity: The percentage of identical matches between the query and subject sequences.
# Alignment Length: The length of the alignment between the query and subject.
# Number of MisMatches: The number of mismatched bases in the alignment.
# Number of Gap Opens: The number of gaps introduced in the alignment.
# Query Start: The starting position of the alignment in the query sequence.
# Query End: The ending position of the alignment in the query sequence.
# Subject Start: The starting position of the alignment in the subject sequence.
# Subject End: The ending position of the alignment in the subject sequence.
# E-value: The expected value, which indicates the number of matches you would expect to see by chance.
# # Bit Score: A score that reflects the quality of the alignment; higher scores indicate better alignments.
# 
# 
# ## Create a key of couple
# blast$query_subject_ids <- paste0(blast$query_id, ":", blast$subject_id)
# blast_bestHit$query_subject_ids <- paste0(blast_bestHit$query_id, ":", blast_bestHit$subject_id)
# 
# # Find the best alignment for each query
# #best_alignments <- blast[order(blast$E_value), ]
# #best_alignments <- best_alignments[!duplicated(best_alignments$query_id), ]
# #length(unique(blast$query_id)); length(unique(blast$subject_id))
# #length(unique(blast_bestHit$query_id)); length(unique(blast_bestHit$subject_id))
# 
# unique(best_alignments$subject_id) %in% unique(blast$subject_id)
# 
# 
# unique(blast$subject_id)[!(unique(blast$subject_id) %in% unique(best_alignments$subject_id))]

