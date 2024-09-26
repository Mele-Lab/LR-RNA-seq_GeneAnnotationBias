## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell"
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
library(ggpmisc)

data <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo.tsv")
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
popcol <- metadata$color_pop
names(popcol) <- metadata$population
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C", "#B5B5B5", "#4F4F4F", "#6E5353", "#000000")
names(colsqanti) <- unique(data$structural_category)[c(1,2,3,5,8,6,7,4)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}


# compute total mapped reads per population
metadata[, population_throughput := sum(map_reads_assemblymap), by=population]


# SAMPLES SHARING------------------------------------------------------------------------------------
sample_sharing <- as.data.frame(data[,.SD, .SDcols = c(grep("isoform", colnames(data)),grep("^[A-Z]{3}[0-9]{1}$",colnames(data)))])
sample_sharing <- column_to_rownames(sample_sharing, var="isoform")
UpSetR::upset(sample_sharing,
              sets = colnames(sample_sharing),
              order.by = "freq")

data$sample_sharing_factor <- ifelse(rowSums(sample_sharing)==1, "1", 
       ifelse(rowSums(sample_sharing)==2, "2", ">=3"))
data[, tools:=factor(ifelse(tool_sharing>1, ">=2", "1"), levels=c("1", ">=2"))]

data[,sample_sharing_factor:=factor(sample_sharing_factor, levels=c("1", "2", ">=3"))]
ggplot(data, aes(x=tools, fill=sample_sharing_factor))+
  geom_bar()+
  mytheme+
  facet_wrap(~structural_category)+
  scale_fill_manual(values=rev(c("#393E41", "#3F88C5", "#44BBA4")))+
  labs(x="Tool sharing", fill="Sample sharing", y="# Transcripts")

mytools <- data[, .(espresso, flair, isoquant, lyric, isoform)]
ggplot(unique(data[, .(isoform, structural_category)][melt(mytools, variable.name = "tool", value.name = "detected", id.vars="isoform"), on="isoform"][detected!=0]),
       aes(x=tool, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  labs(x="", y="# Transcripts")


ggplot(data, aes(x=tools, fill=sample_sharing_factor))+
  geom_bar()+
  mytheme+
  facet_wrap(~structural_category)+
  scale_fill_manual(values=rev(c("#393E41", "#3F88C5", "#44BBA4")))+
  labs(x="Tool sharing", fill="Sample sharing", y="# Transcripts")

# this plots make me wonder about the correlation between mapped reads and total discovered isoforms
sample_discovery <-as.data.frame(colSums(sample_sharing))
colnames(sample_discovery) <- "detected_isoforms"
sample_discovery <- rownames_to_column(sample_discovery, var="sample")
sample_discovery <- metadata[, .(sample, map_reads_assemblymap, population)][sample_discovery, on="sample"]

ggplot(sample_discovery, aes(x=map_reads_assemblymap/10^6, y=detected_isoforms/10^3, color=population))+
  stat_poly_line(color="darkgrey") +
  geom_point()+
  scale_color_manual(values=popcol)+
  xlab("Mapped Reads (M)")+
  ylab("# Detected Transcripts (x1000)")+
  mytheme+
  stat_poly_eq(use_label(c("eq", "R2")), color="black")
# what about transcript sharing across samples
transcript_share <- as.data.frame(table(rowSums(sample_sharing)))
colnames(transcript_share) <- c("samples", "transcripts")
transcript_share$samples <-as.numeric(transcript_share$samples)
ggplot(transcript_share, aes(x=samples, y=transcripts))+geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme

# Define custom breaks
breaks <- c(0, 1,3, 43)
# Add the 'bins' column
transcript_share$bins <- cut(transcript_share$samples, breaks = breaks, include.lowest = F, right = TRUE, labels=c("1", "2-3", ">3"))
setDT(transcript_share)
transcript_share[, transcripts_bin := sum(transcripts), by="bins"]
ggplot(unique(transcript_share[, .(bins, transcripts_bin)]), aes(x=bins, y=transcripts_bin))+geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  geom_text(aes(label=transcripts_bin), vjust=-0.5)+
  ylim(0, max(transcript_share$transcripts_bin) * 1.1)  # Increase the y-axis limit by 10%

# now the same by sqanti category
transcript_sharing_sqanti <- as.data.frame(rowSums(sample_sharing))
colnames(transcript_sharing_sqanti) <- c("samples")
transcript_sharing_sqanti$isoform <- rownames(transcript_sharing_sqanti)

transcript_sharing_sqanti <- data[, .(isoform, structural_category)][transcript_sharing_sqanti, on="isoform"]
transcript_sharing_sqanti[, transcripts := .N, by="structural_category"]

transcript_sharing_sqanti$structural_category <- factor(transcript_sharing_sqanti$structural_category, levels=names(colsqanti))

# plot
ggplot(unique(transcript_sharing_sqanti[, .(samples,isoform,structural_category)]), aes(x=samples, fill=structural_category))+
  geom_bar(position="stack")+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  scale_fill_manual(values=colsqanti)
ggplot(unique(transcript_sharing_sqanti[, .(samples,isoform,structural_category)]), aes(x=samples, fill=structural_category))+
  geom_bar(position="fill")+
  labs(x="# Samples where a transcript is found", y="Proportion of  Transcripts")+
  mytheme+
  scale_fill_manual(values=colsqanti)

# Define custom breaks
breaks <- c(0, 1,3, 43)
# Add the 'bins' column
transcript_sharing_sqanti$bins <- cut(transcript_sharing_sqanti$samples, breaks = breaks, include.lowest = F, right = TRUE, labels=c("1", "2-3", ">3"))
setDT(transcript_sharing_sqanti)
transcript_sharing_sqanti[, transcripts_bin := .N, by=c("bins", "structural_category")]
toplot<- unique(transcript_sharing_sqanti[, .(bins, transcripts_bin, structural_category)])
ggplot(toplot, aes(x=bins, y=transcripts_bin, fill=structural_category))+
  geom_col()+
  labs(x="# Samples where a transcript is found", y="# Transcripts")+
  mytheme+
  geom_text(aes(label=transcripts_bin),position = position_stack(vjust = 0.5))+
  ylim(0, 190000)+
  scale_fill_manual(values=colsqanti)

ggplot(data, aes(x=sample_sharing, fill=factor(tool_sharing)))+
  geom_bar()+
  mytheme+
  labs(x="Sample sharing", fill="Tool sharing", y="# transcripts")+
  scale_fill_manual(values=c("#B44F4F", "#D2C63D", "#579A1D", "#52339A"))


###### POPULATION SPECIFICITY----------------------------------------------------------------------------------------
# compute total mapped reads per population
metadata[, population_throughput := sum(map_reads_assemblymap), by=population]
# Are there population specific transcripts?
popsp <-data[population_sharing==1 & sample_sharing>1, .(isoform,structural_category, AJI, CEU, ITU, HAC, PEL, MPC, YRI, LWK)]
popsp_long <- melt(popsp, variable.name = "population", value.name = "samples")
popsp_long <- popsp_long[samples!=0]
popsp_long <- unique(metadata[,.(population, population_throughput)])[popsp_long, on="population"]
popsp_long[, pop_specific_transcripts:=.N, by=.(population, structural_category)]
popsp_long$structural_category <- factor(popsp_long$structural_category , levels=names(colsqanti))
ggplot(popsp_long, aes(x=reorder(population, population_throughput), fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(x="", y="# Population specific transcripts (>1 sample)")+
  facet_wrap(~structural_category)+
  geom_text(aes(label = after_stat(count), x=population), stat="count", vjust=-0.5)+
  ylim(c(0,420))
ggplot(popsp_long, aes(x=reorder(population, population_throughput), fill=structural_category))+
  geom_bar()+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(x="", y="# Population specific transcripts (>1 sample)")+
  geom_text(aes(label = after_stat(count), x=population), stat="count", position = position_stack(vjust = 0.5))

ggplot(popsp_long, aes(x=reorder(population, population_throughput), fill=structural_category))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(x="", y="# Population specific transcripts (>1 sample)")+
  geom_text(aes(label = after_stat(count), x=population), stat="count", position = position_fill(vjust = 0.5))

popsp_long[, isoform:=NULL]
ggplot(unique(popsp_long), aes( fill=structural_category, y=pop_specific_transcripts, x=population_throughput))+
  geom_area()+
  scale_fill_manual(values=colsqanti)+
  mytheme
ggplot(popsp_long, aes(x=reorder(population, population_throughput), fill=structural_category))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(x="", y="Proportion population specific\ntranscripts (>2 sample)")
rm(popsp, popsp_long, popspupset, pop_specific_transcripts)
####--------------------------------------------------------------------------------------------------------------

### TRANSCRIPT DISCOVERY BY GENE
# is the number of annontated and discovered isoforms related?
trxpergene <- unique(data[!is.na(associated_transcripts_per_gene), .(discovered_transcripts_per_gene, associated_transcripts_per_gene, associated_gene_biotype)])

# Create bins for discovered and associated transcripts using data.table
trxpergene[, All_Discovered_Transcripts := cut(discovered_transcripts_per_gene,
                                               breaks = c(0, 1, 2, 6, seq(10,100, by=10) ,Inf),
                                               labels = c("0", "1", "2-5", "6-9", "10-19", "20-29", "30-39","40-49", "50-59", "60-69", "70-79", "80-89","90-99" ,"100+"),
                                               include.lowest = TRUE)]

trxpergene[, Annotated_Transcripts := cut(associated_transcripts_per_gene,
                                          breaks = c(0, 1, 2, 6, seq(10,100, by=10) ,Inf),
                                          labels = c("0", "1", "2-5", "6-9", "10-19", "20-29", "30-39","40-49", "50-59", "60-69", "70-79", "80-89","90-99" ,"100+"),
                                          include.lowest = TRUE)]

# Count the number of genes in each combination of novel and annotated transcript bins
binned_data_pc <- trxpergene[associated_gene_biotype=="protein_coding"][, .N, by = .(Annotated_Transcripts, All_Discovered_Transcripts)]
binned_data_lnc <- trxpergene[associated_gene_biotype=="lncRNA"][, .N, by = .(Annotated_Transcripts, All_Discovered_Transcripts)]



# Create the heatmap plot
ggplot(binned_data_pc, aes(x = Annotated_Transcripts, y = All_Discovered_Transcripts, fill = N)) +
  geom_tile(color = "white") +                             # Heatmap tiles
  viridis::scale_fill_viridis( breaks = c(0, 25,50, 75,100), na.value = "white") +  # Log scale for colors
  labs(x = "Annotated transcripts/gene", 
       y = "Discovered transcripts/gene") +        # Axis labels
  theme(legend.position = "right")+
  mytheme+
  labs(fill="Protein-coding\ngenes")
ggplot(binned_data_lnc, aes(x = Annotated_Transcripts, y = All_Discovered_Transcripts, fill = N)) +
  geom_tile(color = "white") +                             # Heatmap tiles
  viridis::scale_fill_viridis( breaks = c(0, 25,50, 75,100), na.value = "white") +  # Log scale for colors
  labs(x = "Annotated transcripts/gene", 
       y = "Discovered transcripts/gene") +        # Axis labels
  theme(legend.position = "right")+
  mytheme+
  labs(fill="lncRNA\ngenes")

ggplot(unique(data[associated_gene_biotype%in%c("protein_coding", "lncRNA"), .(associated_transcripts_per_gene,associated_gene_biotype, associated_geneid.v)]),
       aes(x=associated_transcripts_per_gene+1, fill=associated_gene_biotype))+
  geom_density(alpha=0.4)+mytheme+
  labs(x="Log10(annotated transcripts/gene+1)", fill="")+
  scale_x_continuous(
    trans = "log10",
  )
data[, associated_gene_category:=ifelse(structural_category%in%c("full-splice_match", "incomplete-splice_match","novel_in_catalog", "novel_not_in_catalog") & associated_gene_biotype%in%c("protein_coding", "lncRNA"), associated_gene_biotype,
                                        ifelse(structural_category%in%c("full-splice_match", "incomplete-splice_match","novel_in_catalog", "novel_not_in_catalog") & !(associated_gene_biotype%in%c("protein_coding", "lncRNA")), "other_annotated_gene",
                                               ifelse(structural_category=="intergenic", "novel intergenic gene", "fusion/antisense/genic")))]
data[,associated_gene_category:=factor(associated_gene_category, levels=rev(c("protein_coding", "lncRNA", "other_annotated_gene", "fusion/antisense/genic", "novel intergenic gene")))]


ggplot(unique(data[, .(discovered_transcripts_per_gene, associated_gene_category, associated_geneid.v)]),
       aes(x=discovered_transcripts_per_gene+1, y=associated_gene_category, fill=associated_gene_category))+
  ggridges::geom_density_ridges(alpha=0.4)+mytheme+
  guides(fill="none")+
  labs(x="Log10(discovered transcripts/gene+1)", fill="", y="")+
  scale_x_continuous(
    trans = "log10",
  )
ggplot(unique(data[, .(discovered_transcripts_per_gene, associated_gene_biotype, associated_geneid.v,associated_gene_category)]),
       aes(x=discovered_transcripts_per_gene+1, y=associated_gene_biotype, fill=associated_gene_category))+
  ggridges::geom_density_ridges(alpha=0.4)+mytheme+
  labs(x="Log10(discovered transcripts/gene+1)", fill="", y="")+
  scale_x_continuous(
    trans = "log10",
  )
ggplot(data[structural_category=="intergenic", .(discovered_transcripts_per_gene)],
       aes(x=discovered_transcripts_per_gene+1))+
  geom_density(alpha=0.4)+mytheme

### TOOL SHARING ------------------------------------------------------------------------------------------
ggplot(data, aes(x=tool_sharing, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  geom_text(stat="count",aes( label=after_stat(count)),position= position_stack(v=0.5))+
  labs(x="Tool sharing", y="# Transcripts")

tool_sharing <- as.data.frame(data[,.(isoform, lyric, isoquant, flair, espresso)])
tool_sharing <- column_to_rownames(tool_sharing, var="isoform")
UpSetR::upset(tool_sharing,
              sets = colnames(tool_sharing),
              order.by = "freq")
rm(tool_sharing)
tool_sharing2 <- data[, .(isoform, tool_sharing, structural_category)]
tool_sharing2[, tools:=factor(ifelse(tool_sharing>1, ">=2", "1"), levels=c("1", ">=2"))]
tool_sharing2[, structural_category:=factor(structural_category, levels=names(colsqanti))]
ggplot(tool_sharing2, aes(x=tools, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  geom_text(stat="count",aes( label=after_stat(count)),position= position_stack(v=0.5))+
  labs(x="Tool sharing", y="# Transcripts")



# ABOUT SJ----------------------------------------------------------------------------
# is recount support correlated with sample sharing?
library(ggpmisc)
lncrnacutoff <- 25
pccutoff <- 25
cutoff <- 25
ggplot(unique(data[associated_gene_biotype%in%c("protein_coding", "lncRNA"), .(associated_gene_biotype, sj_less_recountsupported_counts, isoform)]), 
       aes(x=associated_gene_biotype, y=sj_less_recountsupported_counts+1, fill=associated_gene_biotype))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.15)+
  mytheme+
  labs(x="", y="Recounts Counts of\nleast supported SJ/transcript")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("darkred", "#356CA1"))+
  geom_hline(yintercept=lncrnacutoff, col="darkred", linetype="dashed")+
  geom_hline(yintercept=pccutoff, col="#356CA1", linetype="dashed")+
  # Proportion above threshold
  annotate("text", col = "#356CA1", x = "protein_coding", y = 1e9, label = paste0("pass ",round(nrow(data[associated_gene_biotype %in% c("protein_coding") & sj_less_recountsupported_counts >= pccutoff]) / nrow(data[associated_gene_biotype %in% c("protein_coding")]), digits = 2) * 100, " %"), vjust = 0) +
  annotate("text", col = "darkred", x = "lncRNA", y = 1e9, label = paste0("pass ",round(nrow(data[associated_gene_biotype %in% c("lncRNA") & sj_less_recountsupported_counts >= lncrnacutoff]) / nrow(data[associated_gene_biotype %in% c("lncRNA")]), digits = 2) * 100, " %"), vjust = 0) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=10))+
  # Threshold labels
  annotate("text", col = "#356CA1", x = 2.15, y = pccutoff, label = pccutoff, vjust = -0.5, hjust = 0.5) +
  annotate("text", col = "darkred", x = 1.15, y = lncrnacutoff, label = lncrnacutoff, vjust = -0.5, hjust = 0.5)+
  annotation_logticks(sides="l")

data[, structural_category:=factor(structural_category, levels=names(colsqanti))]
data[, pass_percentage := round(
  sum(sj_less_recountsupported_counts >= cutoff) / .N * 100, 2), 
  by = structural_category]

ggplot(unique(data[, .(structural_category, sj_less_recountsupported_counts, isoform)]), 
       aes(x=structural_category, y=sj_less_recountsupported_counts+1, fill=structural_category))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.05)+
  mytheme+
  labs(x="", y="Recounts Counts of\nleast supported SJ/transcript")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=colsqanti)+
  annotation_logticks(sides="l")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=10))+
  geom_hline(yintercept=cutoff, col="black", linetype="dashed")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  geom_text(data = unique(data[, .(structural_category, pass_percentage)]), 
            aes(x = structural_category, y = 1e9, label = paste0("pass ", pass_percentage, " %"), 
                color = structural_category), 
            vjust = 0)+
  scale_color_manual(values = colsqanti)+
  annotate("text", label=">=25", y=85, x=0.75)

ggplot(data, aes(x=sample_sharing, y=log10(sj_less_recountsupported_counts+1), color=structural_category))+
  stat_poly_line(color="darkred") +
  geom_point(alpha=0.2)+
  mytheme+
  ylab("Recount3 counts of\nleast supported sj")+
  scale_color_manual(values=colsqanti)+
  facet_wrap(~structural_category)+
  stat_poly_eq(use_label(c("eq", "R2")), color="black")
library(ggridges)

# how is the distribution of recount support by structural category?
data$structural_category <- factor(data$structural_category, levels=rev(names(colsqanti)))
ggplot(data, aes( x=log10(sj_less_recountsupported_counts+1), y=structural_category, fill=structural_category))+
  ggridges::geom_density_ridges() +  
  mytheme+
  xlab("Recount3 log10(counts) of\nleast supported SJ")+
  ylab("")+
  scale_fill_manual(values=colsqanti)

ggplot(data, aes( x=log10(sj_less_recountsupported_counts+1), y=structural_category, fill=structural_category))+
  ggridges::geom_density_ridges() +  
  mytheme+
  xlab("Recount3 log10(counts) of\nleast supported SJ")+
  ylab("")+
  scale_fill_manual(values=colsqanti)

breaks <- c(0, 1,3, 43)


data[, structural_category:=factor(structural_category, levels=names(colsqanti))]
# filter by expression levels

ggplot(data, aes(x=log10(flair_total_counts+1), y=log10(sj_less_recountsupported_counts+1), col=structural_category))+
  geom_point(alpha=0.3)+
  stat_poly_line(color="darkred") +
  theme_bw()+
  labs(y="Recount counts of less supported SJ\nlog10(y+1)", x="Total counts\nlog10(x+1)")+
  facet_wrap(~structural_category)+
  scale_color_manual(values=colsqanti)+
  stat_poly_eq(use_label(c("eq", "R2")), color="black")

data[, sj_less_recount_counts_info:=factor(sj_less_recount_counts_info, levels=c("knownSJ", "novelSJ_knownSS_2", "novelSJ_knownSS_1", "novelSJ_knownSS_0"))]
ggplot(data, aes(x=sj_less_recount_counts_info, fill=all_canonical))+
  geom_bar()+
  mytheme+
  labs(x="SJ with the least counts\nin recount3 per transcript", y="# transcripts", fill="")+
  scale_fill_manual(values=c("#3C6822", "darkred"))+
  facet_wrap(~structural_category)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggplot(data, aes(x=sj_less_recount_counts_info, fill=all_canonical, y=log10(flair_total_counts+1)))+
  geom_boxplot(outliers = F)+
  mytheme+
  labs(x="SJ with the least counts\nin recount3 per transcript", y="Total counts log10(y+1)", fill="")+
  scale_fill_manual(values=c("#3C6822", "darkred"))+
  facet_wrap(~structural_category)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# GENCODE junctions to find recount3 threshold
# load and parse data
gencodesj <- fread("04_transcriptome_assembly/04_evaluation/02_sqanti/data/gencodev47/gencodev47_junctions.txt")
annot <- fread("../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
annot <- annot[, .(transcriptid.v, geneid.v, gene_biotype)]
gencodesj[, junction:=paste(chrom, genomic_start_coord, genomic_end_coord, strand, sep=":")]
gencodesj <- annot[gencodesj[, .(isoform, junction, junction_category, start_site_category, end_site_category, canonical)], on=c("transcriptid.v"="isoform")]
recount <- fread("04_transcriptome_assembly/04_evaluation/03_recount3/data/recount3_srv3h_morethan10counts.tsv")
recount <- recount[, .(V1, V2, V3, V4, V8, V9)]
colnames(recount) <- c("contig", "start", "end", "strand", "samples", "counts")

# # create correspondence between scaffold names and substitute them in recount data
# mycorrespondence <- unique(recount$V1)
# names(mycorrespondence) <- unique(recount$V1)
# mycorrespondence <- gsub("v","\\.",gsub(".*_","",gsub("_random", "", mycorrespondence)))
# fwrite(data.frame(mycorrespondence, names(mycorrespondence)), "04_transcriptome_assembly/04_evaluation/01_adapt_lyric_output/data/ref_contigs_correspondence.tsv",quote=F, row.names=F, col.names = F, sep="\t")
# recount[, contig := mycorrespondence[contig]]
# recount[, junction := paste(contig, start, end, strand, sep=":")]
# 
# # add recount data to sj
# gencodesj_recount <- recount[gencodesj, on="junction"]
# gencodesj_recount$counts <-replace_na(gencodesj_recount$counts, 0)
# fwrite(gencodesj_recount, "04_transcriptome_assembly/04_evaluation/03_recount3/data/gencodev47_junctions_with_recount_counts.tsv", row.names = F, quote=F, sep="\t")
gencodesj_recount <- fread("04_transcriptome_assembly/04_evaluation/03_recount3/data/gencodev47_junctions_with_recount_counts.tsv")

lncrnacutoff <- 25
pccutoff <- 1000
ggplot(gencodesj_recount[gene_biotype%in%c("protein_coding", "lncRNA")], aes(x=gene_biotype, y=counts+1, fill=gene_biotype))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.15)+
  mytheme+
  labs(x="", y="Recounts Counts")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("darkred", "#356CA1"))+
  geom_hline(yintercept=lncrnacutoff, col="darkred", linetype="dashed")+
  geom_hline(yintercept=pccutoff, col="#356CA1", linetype="dashed")+
  annotate("text", col="#356CA1", x="protein_coding", y=pccutoff, label=round(nrow(gencodesj_recount[gene_biotype%in%c("protein_coding") & counts>pccutoff])/nrow(gencodesj_recount[gene_biotype%in%c("protein_coding")]), digits=2), vjust=0, hjust=0)+
  annotate("text", col="darkred", x="lncRNA", y=lncrnacutoff, label=round(nrow(gencodesj_recount[gene_biotype%in%c("lncRNA") & counts>lncrnacutoff])/nrow(gencodesj_recount[gene_biotype%in%c("lncRNA")]), digits=2), vjust=0, hjust=0)+
  annotate("text", col="#356CA1", x=0, y=pccutoff, label=pccutoff, vjust=0, hjust=0)+
  annotate("text", col="darkred", x=0, y=lncrnacutoff, label=lncrnacutoff, vjust=0, hjust=0)+
  annotation_logticks(sides="l")

gencodesub <- gencodesj_recount[gene_biotype%in%c("protein_coding", "lncRNA")]
gencodesub[, minsupport := min(counts), by = transcriptid.v]

lncrnacutoff <- 25
pccutoff <- 25
ggplot(unique(gencodesub[, .(gene_biotype, minsupport, transcriptid.v)]), aes(x=gene_biotype, y=minsupport+1, fill=gene_biotype))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.15)+
  mytheme+
  labs(x="", y="Recounts Counts of\nleast supported SJ/transcript")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("darkred", "#356CA1"))+
  geom_hline(yintercept=lncrnacutoff, col="darkred", linetype="dashed")+
  geom_hline(yintercept=pccutoff, col="#356CA1", linetype="dashed")+
  annotate("text", col="#356CA1", x="protein_coding", y=pccutoff, label=round(nrow(gencodesub[gene_biotype%in%c("protein_coding") & minsupport>=pccutoff])/nrow(gencodesub[gene_biotype%in%c("protein_coding")]), digits=2), vjust=0, hjust=0)+
  annotate("text", col="darkred", x="lncRNA", y=lncrnacutoff, label=round(nrow(gencodesub[gene_biotype%in%c("lncRNA") & minsupport>=lncrnacutoff])/nrow(gencodesub[gene_biotype%in%c("lncRNA")]), digits=2), vjust=0, hjust=0)+
  annotate("text", col="#356CA1", x=0, y=pccutoff, label=pccutoff, vjust=0, hjust=0)+
  annotate("text", col="darkred", x=0, y=lncrnacutoff, label=lncrnacutoff, vjust=0, hjust=0)+
  annotation_logticks(sides="l")

#### QUANTIFICATION LEVEL EXPLORATION-----------------------------------------------------------------
ggplot(annot, aes(x=sample_sharing, y=flair_expressed_samples))+
  geom_bin2d(bins=42)+
  viridis::scale_fill_viridis(trans = "log")+
  mytheme+
  labs(x="# Samples assembling a transcript", y="# Samples with transcrips expression>=1")

annot[, detection_thresholds:=ifelse(is.na(flair_expressed_samples), "never expressed",
                                     ifelse(sample_sharing==flair_expressed_samples, "equal",
                                            ifelse(sample_sharing>flair_expressed_samples, "more detection\nthan expression",
                                                   ifelse(sample_sharing<flair_expressed_samples, "more expression\nthan detection", NA))))]

ggplot(annot, aes(x=detection_thresholds, fill=detection_thresholds))+
  geom_bar(stat = "count")+
  mytheme+
  labs(x="", y="# transcripts")+
  geom_text(aes(label = after_stat(count)), stat="count", vjust=-0.5)+
  ylim(c(0,200000))+
  scale_fill_manual(values=c("#3C6822", "#7F4089", "#224870", "#BD2D1F"))


ggplot(annot, aes(x=factor(tool_sharing), y=log10(flair_total_counts+1), fill=factor(tool_sharing)))+
  geom_boxplot(width=0.15, outliers = F)+
  geom_violin(scale="width", trim=T, alpha=0.5)+
  labs(x="Tool sharing", y="Total counts\n(only in expressed samples)")+
  mytheme+
  scale_fill_manual(values=c("#B44F4F", "#D2C63D", "#579A1D", "#52339A"))+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0))
ggplot(annot, aes(x=factor(tool_sharing), y=log10(flair_mean_counts+1), fill=factor(tool_sharing)))+
  geom_boxplot(width=0.15, outliers = F)+
  geom_violin(scale="width", trim=T, alpha=0.5)+
  labs(x="Tool sharing", y="Mean counts\n(only in expressed samples)")+
  mytheme+
  scale_fill_manual(values=c("#B44F4F", "#D2C63D", "#579A1D", "#52339A"))+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0))


ggplot(data, aes(x=structural_category, y=flair_expressed_samples, fill=structural_category))+
  geom_violin()+
  scale_fill_manual(values=colsqanti)


#### TEST A SET OF FILTERS TO GET A SENSE OF HOW MUCH WOULD WE BE REMOVING
ggplot(data[c(sample_sharing>1 | tool_sharing>1 & flair_total_counts>=10), ], 
       aes(x=structural_category, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=colsqanti)

# # How many transcripts are there associated to each gene by sqanti category
# geneperstrcount <-unique(data[, .(trxPerGeneAndCategory = .N), by=c("associated_geneid.v", "structural_category")])
# transcriptperstrcount <-unique(data[, .(trxPerAssTrxAndCategory = .N), by=c("associated_transcriptid.v", "structural_category")])
# 
# 
# geneperstrcount_wide <- dcast(unique(geneperstrcount[, .(associated_geneid.v, structural_category, trxPerGeneAndCategory)]), associated_geneid.v~structural_category, fill=0)
# transcriptperstrcount_wide <- dcast(unique(transcriptperstrcount[, .(associated_transcriptid.v, structural_category, trxPerAssTrxAndCategory)]), associated_transcriptid.v~structural_category, fill=0)
# colnames(geneperstrcount_wide)[2:ncol(geneperstrcount_wide)] <- paste0("trx_pergene_count_",colnames(geneperstrcount_wide)[2:ncol(geneperstrcount_wide)] )
# colnames(transcriptperstrcount_wide)[2:ncol(transcriptperstrcount_wide)] <- paste0("trx_per_asstrx_count_",colnames(transcriptperstrcount_wide)[2:ncol(transcriptperstrcount_wide)] )
# 
# data <- geneperstrcount_wide[data, on = .(associated_geneid.v)]
# data <- transcriptperstrcount_wide[data, on = .(associated_transcriptid.v)]
# 
# data[, existsFSMinGene:=ifelse(`trx_pergene_count_full-splice_match`>0, TRUE, FALSE)]
# data[, existsFSMinTranscript:=ifelse(`trx_per_asstrx_count_full-splice_match`>0, TRUE, FALSE)]

# Check if ISM missing FSM have any other transcript in the gene (IOW do we have to include the gene entry?)


# do we have the corresponding FSM of ISM associated transcripts?
ggplot(data[structural_category=="incomplete-splice_match"], aes(x=ref_exons, fill=existsFSMinTranscript))+
  geom_bar()+
  mytheme+
  xlim(c(0,50))+
  labs(x="# Exons of associated transcript", y="# Transcript Models", fill="Matching FSM\ndetected")
# how many genes have FSM?
ggplot(data, aes(x=associated_transcripts_per_gene, y=`trx_pergene_count_full-splice_match`))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +  mytheme+
  scale_y_continuous(transform = "log10")+
  scale_x_continuous(transform = "log10")+
  scale_fill_gradient2(low="white", high="darkred")+
  labs(x="# Annotated transcripts/gene", y="# FSM/gene")

# find filters for quantification
# data$flair_max_counts <- replace_na(data$flair_max_counts, 0)
cutoff <- 3
data[, pass_quanti_percentage := round(sum(flair_max_counts>=cutoff)*100/.N,  digits=2), by=structural_category]
ggplot(data, aes(x=structural_category, y=flair_max_counts+1, fill=structural_category))+
  geom_violin( alpha=0.7,adjust = 4)+
  geom_boxplot(outliers = F, width=0.1)+
  labs(x="", y="FLAIR max counts")+
  guides(fill="none", color="none")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=colsqanti)+
  annotation_logticks(sides="l")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=6))+
  geom_hline(yintercept=cutoff, col="black", linetype="dashed")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  geom_text(data = unique(data[, .(structural_category, pass_quanti_percentage)]), 
            aes(x = structural_category, y = 0.5e6, label = paste0("pass ", pass_quanti_percentage, " %"), 
                color = structural_category), 
            vjust = 0)+
  scale_color_manual(values = colsqanti)+
  annotate("text", label=as.character(cutoff), y=cutoff, x=0.75)+
  mytheme

  
  
## EXPLORE ISM
ggplot(data[structural_category=="incomplete-splice_match" & `trx_per_asstrx_count_full-splice_match`==FALSE], 
       aes(x=flair_mean_counts))+
  mytheme+
  geom_density(adjust=2)+
  xlab("")+
xlim(c(0,10))

ggplot(data[structural_category=="incomplete-splice_match" & `trx_per_asstrx_count_full-splice_match`==FALSE], 
       aes(y=flair_mean_counts, x=sample_sharing))+
  mytheme+
  geom_jitter(alpha=0.5)+
  xlab("")+
  ylim(c(0,10))
ggplot(data[structural_category=="incomplete-splice_match" & `trx_per_asstrx_count_full-splice_match`==FALSE], 
       aes(x=`trx_per_asstrx_count_incomplete-splice_match`))+
  mytheme+
  geom_bar()+
  labs(x="# ISM/ISM associated transcript missing FSM", y="# Transcript Models")
ggplot(data[structural_category=="incomplete-splice_match" & `trx_per_asstrx_count_full-splice_match`==FALSE], 
       aes(y=`trx_per_asstrx_count_incomplete-splice_match`, x=ref_exons))+
  mytheme+
  geom_point(alpha = 0.3, color="grey")+
  geom_density2d(color="black")+
  labs(y="# ISM/ISM associated transcript missing FSM", x="# Exons in Full transcript")
##### EXPLORE LENGTHS
cutoff <- 300
data[, pass_length_percentage := round(sum(length>=cutoff)*100/.N,  digits=2), by=structural_category]
ggplot(data[filter=="pass"], aes(x=structural_category, y=length, fill=structural_category))+
  geom_violin(scale="width", alpha=0.7)+
  geom_boxplot(outliers = F, width=0.1)+
  guides(fill="none", color="none")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=colsqanti)+
  annotation_logticks(sides="l")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=6))+
  geom_hline(yintercept=cutoff, col="black", linetype="dashed")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  geom_text(data = unique(data[, .(structural_category, pass_length_percentage)]), 
            aes(x = structural_category, y = 0.5e6, label = paste0("pass ", pass_length_percentage, " %"), 
                color = structural_category), 
            vjust = 0)+
  scale_color_manual(values = colsqanti)+
  annotate("text", label=as.character(cutoff), y=cutoff, x=0.75)+
  mytheme
