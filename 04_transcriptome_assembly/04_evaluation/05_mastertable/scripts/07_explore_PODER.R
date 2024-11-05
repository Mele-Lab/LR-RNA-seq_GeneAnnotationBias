## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/05_mastertable"
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
# Load data
data <- fread("data/29102024_PODER_mastertable.tsv")
metadata <- fread("../../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,4,2,6,7,3,5)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
data[, associated_gene_biotype := factor(associated_gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene"))]

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
median_fun <- function(x,y) {
  return(data.frame(y = y, label = paste0("Median= ", round(10^median(x), 0))))
}

# ADAPTATIONS
data[, predicted_NMD:=fifelse(predicted_NMD==TRUE, "Predicted NMD", "Not NMD")]
data[, proteinv47_protein_is_nmd:=fifelse(proteinv47_protein_is_nmd==TRUE, "Predicted NMD", "Not NMD")]

# LOOK FOR NMD tags
annot <- fread("../../../../../../Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf")
annot <- annot[V3=="transcript", .(V9)]
annot[, transcriptid.v:=tstrsplit(V9, "\"")[[4]]]
annot[, gene_biotype:=tstrsplit(V9, "\"")[[6]]]
annot <- annot[gene_biotype=="protein_coding"]
annot[, transcript_biotype:=tstrsplit(V9, "\"")[[10]]]
nmd <- annot[, .(transcriptid.v, transcript_biotype)]

data <- nmd[data, on=c("transcriptid.v"="associated_transcript")]
colnames(data)[grep("transcriptid.v", colnames(data))] <- "associated_transcript"
data[, structural_category_nmd:=factor(fifelse(transcript_biotype=="nonsense_mediated_decay", 
                                        "FSM\nNMD\nannotated", 
                                        fifelse(structural_category=="FSM", "FSM\nnot NMD\nannotated",as.character(structural_category)), as.character(structural_category)),
                                       levels=c("FSM\nNMD\nannotated", "FSM\nnot NMD\nannotated", "NIC", "NNC", "Fusion"))]

# # plot result from PODER
ggplot(data, aes(x=structural_category, fill=structural_category))+
  geom_bar(aes(alpha=annotated))+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Transcripts in PODER", alpha="Transcript Source")+
  mytheme+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-0.5)+
  scale_alpha_manual(values=c(0.6, 1),
                     labels=c("annotated"="GENCODEv47\n(ISM replacement)", "discovered"="Transcript Discovery Tools"))+
  annotate("text", x="FSM", y=112000, label=paste0("(",nrow(data[annotated=="annotated"]), ")"), size=3)+
  theme(legend.position= c(0.8, 0.9))
       # legend.background = element_rect(fill = "white", color = "black"))
ggsave("../../../10_figures/fig_02/barplot_PODER_PerSqantiCategory_COUNT.pdf", dpi=700, width = 18, height = 14,  units = "cm")


ggplot(unique(data[,.(isoform, associated_gene_biotype, structural_category)]),
       aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#0a9396","#ee9b00","darkgrey"),
                    labels=c("protein_coding"="Protein Coding", "novel/ambiguous gene"="Novel/Ambiguous Gene"))+
  labs(x="", y="Proportion PODER Transcripts", fill="Gene Biotype")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))+
  theme(legend.position = "top")
ggsave("../../../10_figures/fig_02/barplot_PODER_PerSqantiCategory_FreqBiotype_trxlvl.pdf", dpi=700, width = 16, height = 14,  units = "cm")

ggplot(unique(data[,.(associated_gene, associated_gene_biotype, structural_category)]),
       aes(x=structural_category, fill=associated_gene_biotype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#0a9396","#ee9b00","darkgrey"),
                    labels=c("protein_coding"="Protein Coding", "novel/ambiguous gene"="Novel/Ambiguous Gene"))+
  labs(x="", y="Proportion PODER Genes", fill="Gene Biotype")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_fill(vjust = 0.5))+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/barplot_PODER_PerSqantiCategory_FreqBiotype_genelvl.pdf", dpi=700, width = 16, height = 14,  units = "cm")

# exons/ transcript
ggplot(data[, .(num_exons = factor(ifelse(exons >= 25, ">=25", as.character(exons)), levels = c(as.character(1:24), ">=25")),
                structural_category,
                isoform)],
       aes(x = num_exons, fill = structural_category)) +
  geom_bar() +
  mytheme +
  facet_wrap(~structural_category, scales = "free_y", nrow=2) +
  scale_fill_manual(values = colsqanti)+
  scale_x_discrete(breaks = c(as.character(seq(0, 21, by = 5)), ">=25"))+
  guides(fill="none")+
  labs(x="Exons/Transcript", y="# Transcripts in PODER")+
  theme(plot.margin = margin(20, 20, 30, 20))
ggsave("../../../10_figures/suppfig_06/barplot_PODER_PerSqantiCategory_exonsPerTranscript.pdf", dpi=700, width = 20, height = 14,  units = "cm")

# PLOT LENGTH
ggplot(data, aes(y = length, x = structural_category, fill = structural_category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(outliers = FALSE, width = 0.1) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = colsqanti) +
  mytheme +
  labs(y = "Transcript length", x = "") +
  annotation_logticks(sides = "l") +
  guides(fill = "none") +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = 1.2)) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 1.4)) # Adjust y for median labels
ggsave("../../../10_figures/suppfig_07/violin_TranscriptLenght_PerSqantiCategory.pdf", dpi=700, width = 21, height = 14,  units = "cm")




# data[, ORF_seq_match:=fifelse(ORF_seq==proteinv47_protein_sequence, "TRUE", "FALSE")]
# data[, ORF_prediction_match:=fifelse(predicted_ORF==proteinv47_predicted_ORF, "TRUE", "FALSE")]
# data[, NMD_prediction_match:=fifelse(predicted_NMD==proteinv47_protein_is_nmd, "TRUE", "FALSE")]
data[, associated_gene_biotype_sub:=factor(associated_gene_biotype_sub, levels=c("Protein Coding", "lncRNA", "Fusion Gene", "Novel Gene"))]

# PLOT ORF PREDICTION
ggplot(data, aes(x=associated_gene_biotype_sub, fill=proteinv47_predicted_ORF))+
  geom_bar(position="fill", alpha=.95)+
  mytheme+
  labs(x="", y="Proportion PODER Transcripts", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))+
  scale_fill_manual(values=c("#c44536", "#457b9d"))+
  theme(legend.position = "top")  # Change legend position to top
ggsave("../../../10_figures/suppfig/barplot_ORFPrediction_PerBiotype.pdf", dpi=700, width = 12, height = 14,  units = "cm")


# PLOT BLASTP BITSCORE
ggplot(data, aes(y=proteinv47_blastp_bitscore, x=associated_gene_biotype_sub, fill=associated_gene_biotype))+
  geom_violin(alpha=0.7, scale="width")+
  geom_boxplot(outliers=F, width=0.07, show.legend = FALSE)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=4.3))+
  scale_fill_manual(values=c("#ee9b00","#0a9396", "darkgrey"))+
  labs(x="", y="blastp Bitscore", fill="Gene Biotype")+
  scale_y_continuous(trans="log10")+
  annotation_logticks(sides="l")+
  theme(legend.position = "top")+
  guides(fill = guide_legend(override.aes = list(color = NA)))
ggsave("../../../10_figures/suppfig/violin_blastpBitscore_PerBiotype.pdf", dpi=700, width = 16, height = 12,  units = "cm")

# PLOT BLASTP IDENTITY
ggplot(data, aes(y=proteinv47_blastp_identity, x=associated_gene_biotype_sub, fill=associated_gene_biotype))+
  geom_violin(alpha=0.7, scale="width")+
  geom_boxplot(outliers=F, width=0.07, show.legend = FALSE)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=25))+
  scale_fill_manual(values=c("#ee9b00","#0a9396", "darkgrey"))+
  labs(x="", y="% blastp Identity", fill="Gene Biotype")+
  theme(legend.position = "top")+
  guides(fill = guide_legend(override.aes = list(color = NA)))
ggsave("../../../10_figures/suppfig/violin_blastpIdentity_PerBiotype.pdf", dpi=700, width = 16, height = 12,  units = "cm")



# PLOT PREDICTED NMD
ggplot(data[!is.na(proteinv47_protein_is_nmd) & associated_gene_biotype_sub%in%c("Protein Coding", "Fusion Gene") & proteinv47_predicted_ORF=="Predicted ORF"],
       aes(x=structural_category_nmd, fill=proteinv47_protein_is_nmd))+
  geom_bar(position="fill", alpha=.8)+
  mytheme+
  labs(x="", y="Proportion PODER Transcripts", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))+
  scale_fill_manual(values=rev(c("#c44536", "#457b9d")),
                    labels=c("FALSE"="Not NMD Predicted", "TRUE"="NMD Predicted"))+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/barplot_NMDPrediction_PerBiotype.pdf", dpi=700, width = 12, height = 14,  units = "cm")

# POPULATION SPECIFIC TRANSCRIPTS
popsp <-data[sample_sharing>=3 & population_sharing==1, ]
cols_to_keep <- colnames(popsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(popsp))]
subpopsp <- popsp[, c(cols_to_keep, "sample_sharing", "structural_category", "isoform"), with = FALSE]
subpopsplong <- melt(subpopsp, measure.vars= colnames(subpopsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(subpopsp))], variable.name = "sample", value.name = "detected")
subpopsplong <- subpopsplong[!is.na(detected)]
subpopsplong <- subpopsplong[detected==1]

subpopsplongmeta <- metadata[, .(sample, map_reads_assemblymap,population )][, total_throughput:=sum(map_reads_assemblymap), by="population"][subpopsplong, on="sample"]
subpopsplongmeta[, trx_per_cat_per_pop := uniqueN(isoform), by=c("population", "structural_category")]
subpopsplongmeta[, trx_per_cat_per_pop_norm :=trx_per_cat_per_pop/total_throughput*10^6]
subpopsplongmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]

# LINE PLOT FOR 2 SAMPLES SHARING FACETED
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  geom_segment(
    aes(xend = total_throughput/10^6, # Adjust for arrow length
        yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
    arrow = arrow(length = unit(0.2, "cm")), # Arrow size
    color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(
    aes(x = total_throughput/10^6,  # Position text away from point
        y = ifelse(eur == "European", trx_per_cat_per_pop + 30, trx_per_cat_per_pop - 30),
        label = population,
        color=population),
    fontface="bold",
    family="Helvetica",
    size = 4.5,
    alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/line_PODER_popSpecific_transcripts.2samplesSharing.faceted.pdf", dpi=700, width = 24, height = 20,  units = "cm")

# LINE PLOT FOR 2 SAMPLES SHARING
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(col="", alpha="", size="")+
  geom_segment(data = unique(subpopsplongmeta[structural_category == "FSM" ,.(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
    aes(xend = total_throughput/10^6, # Adjust for arrow length
        yend = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
    arrow = arrow(length = unit(0.2, "cm")), # Arrow size
    color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data = unique(subpopsplongmeta[structural_category == "FSM", .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
    aes(x = total_throughput/10^6,  # Position text away from point
        y = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 30, trx_per_cat_per_pop - 30),
        label = population,
        color=population),
    fontface="bold",
    family="Helvetica",
    size = 4.5,
    alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/line_PODER_popSpecific_transcripts.2samplesSharing.pdf", dpi=700, width = 16, height = 14,  units = "cm")



# LINE PLOT FOR 3 SAMPLES SHARING 
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(col="", alpha="", size="")+
  geom_segment(data = unique(subpopsplongmeta[structural_category == "FSM" ,.(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 1, trx_per_cat_per_pop - 1)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data = unique(subpopsplongmeta[structural_category == "FSM", .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 1.2, trx_per_cat_per_pop - 1.2),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 4.5,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/line_PODER_popSpecific_transcripts.3samplesSharing.pdf", dpi=700, width = 16, height = 14,  units = "cm")
