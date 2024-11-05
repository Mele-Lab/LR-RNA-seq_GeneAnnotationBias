## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/08_allele_specifics"
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
library(ggpubr)
metadataraw <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadataraw <- metadataraw[mixed_samples==FALSE]
metadata <- metadataraw[merged_run_mode==TRUE]
metadata[, samplecode:=paste(lab_number_sample, lab_sampleid, cell_line_id, sep="_")]

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)

n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
# load data
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
arrayenh <- fread("array_enhanced_gencode", header = F)
array <- rbind.data.frame(arraygen, arraypan)
array <- rbind.data.frame(array, arrayenh)

asts <- list()

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  temp <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_annotated_taggedresults.tsv"))
  temp[, `:=`(annot=TYPE, samplecode=SAMPLE)]
  asts <- append(asts, list(temp))}

astsraw <- rbindlist(asts, use.names=TRUE)
astsraw <- metadata[, .(samplecode, sample, population, map_reads_assemblymap)][astsraw, on="samplecode"]
astsraw[, annot := ifelse(annot=="gencode", "GENCODEv47", 
                          ifelse(annot=="enhanced_gencode", "Enhanced\nGENCODEv47", "PODER"))]
astsraw[, annot:=factor(annot, levels=c("GENCODEv47", "PODER", "Enhanced\nGENCODEv47"))]
# fwrite(astsraw, "data/ASTS_results_threeannots.tsv", quote=F, row.names = F, sep="\t")
# astsraw <- fread("data/ASTS_results_threeannots.tsv")
asts <- astsraw[gene_testable==TRUE]

# computed number of tested genes
asts[, tested_genes:=uniqueN(geneid.v), by=c("sample", "annot")]
unique_sig_count <- asts[FDR < 0.05, .(sig_genes=uniqueN(geneid.v)), by=c("sample", "annot")]
asts <- unique_sig_count[asts, on=c("sample", "annot")]

asts[, afr:=fifelse(population%in%c("YRI", "LWK", "MPC"), "African", "OOA")]
ggplot(unique(asts[, .(sig_genes, tested_genes, population, map_reads_assemblymap, annot, sample)]), 
       aes(x=tested_genes, y=sig_genes))+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.9, size = 5) +  # Equation and R-squared
  geom_point(aes(col=population,size=map_reads_assemblymap/10^6), alpha=0.7)+
  mytheme+
  labs(y="# ASTU Significant Genes", x="# ASTU Tested Genes", size="Reads (M)", col="Population")+
  scale_color_manual(values=popcols)+
  facet_wrap(~annot)
# number of tested genes per annot?

ggplot(unique(asts[, .(sample, annot, sig_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=sig_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers = F, width=0.1)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Significant ASTU Genes", size="Mapped Reads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("PODER", "GENCODEv47"), 
                                        c("Enhanced\nGENCODEv47", "PODER"),
                                        c("Enhanced\nGENCODEv47", "GENCODEv47")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE))+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")
ggplot(unique(asts[, .(sample, annot, tested_genes, afr, population, map_reads_assemblymap)]), 
       aes(x=annot, y=tested_genes, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers = F, width=0.1)+
  ggbeeswarm::geom_quasirandom(aes(color=population, size=map_reads_assemblymap/10^6),alpha=0.8)+
  scale_color_manual(values=c(popcols))+
  mytheme+
  labs(x="", y="# Tested Genes", size="Mapped Reads (M)", col="Population", fill="")+
  stat_compare_means(comparisons = list(c("PODER", "GENCODEv47"), 
                                        c("Enhanced\nGENCODEv47", "PODER"),
                                        c("Enhanced\nGENCODEv47", "GENCODEv47")),method = "t.test",
                     method.args = list(alternative = "two.sided", paired=TRUE))+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59", "#62531C"))+
  guides(fill="none")



# Explore how many more genes are tested in each annot
astswide <- dcast(unique(asts[, .(sample, annot, tested_genes)]), sample ~ annot, value.var = "tested_genes") 
astswide[, diff_tested_genes:=PODER-GENCODEv47]
astswide[, ratio_tested_genes:=PODER/GENCODEv47]

astswidemeta <- metadata[, .(sample, population, map_reads_assemblymap)][astswide, on="sample"]
astswidemeta[, eur:=ifelse(population %in% c("LWK", "YRI"), "AFR",
                           ifelse(population %in% c("PEL", "HAC", "ITU"), "nonEUR-OOA", "EUR"))]
ggplot(astswidemeta, aes(x=sample, y=diff_tested_genes, fill=population))+
  geom_col()+
  mytheme+
  scale_fill_manual(values=popcols)
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=diff_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Increment of genes tested by PODER against GENCODE")+
  ylim(c(0,115))
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=ratio_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Ratio of genes tested by PODER against GENCODE")+
  geom_hline(yintercept=1, linetype="dashed", color="black")

# are EUR beign the least benefitted?
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=diff_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Increment of genes tested by PODER against GENCODE")+
  ylim(c(0,115))+
  facet_wrap(~eur)
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=ratio_tested_genes))+
  mytheme+
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Ratio of genes tested by PODER against GENCODE")+
  geom_hline(yintercept=1, linetype="dashed", color="black")+
  facet_wrap(~population)


# check  normality and homocedasticity
car::leveneTest(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
bartlett.test(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
shapiro.test(astswidemeta$ratio_tested_genes)
qqnorm(astswidemeta$ratio_tested_genes)
qqline(astswidemeta$ratio_tested_genes, col = "red")
my_comparisons <- list( c("CEU", "HAC"), c("CEU", "ITU"), c("CEU", "LWK"),c("CEU", "PEL"), c("CEU", "YRI") )
ggplot(astswidemeta, aes(x=population, y=ratio_sig_genes))+
  mytheme+
  geom_boxplot(outliers=F)+
  geom_jitter(alpha=0.8, aes(col=population, size=map_reads_assemblymap/10^6))+
  scale_color_manual(values=popcols)+
  labs(x="", y="Ratio of genes tested by PODER against GENCODE", size="Mapped Reads (M)", col="Population")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                     method.args = list(alternative = "less")) # Add pairwise comparisons p-value

ggplot(astswidemeta, aes(x=population, y=map_reads_assemblymap/10^6))+
  mytheme+
  geom_boxplot(outliers=F)+
  geom_jitter(alpha=0.8, aes(col=population, size=ratio_sig_genes))+
  scale_color_manual(values=popcols)+
  labs(x="", y="Mapped reads (M)", size="Ratio tested genes\nPODER/GENCODE", col="Population")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",
                     method.args = list(alternative = "less")) 

# check biotypes of tested genes
astsgenes <- unique(asts[, .(annot, sample, population, map_reads_assemblymap, geneid.v, gene_biotype, gene_testable, sig_genes)])

ggplot(unique(astsgenes[, .(annot, ggplot(unique(astsgenes[, .(annot, ggplot(unique(astsgenes[, .(annot, geneid.v, gene_biotype)]), aes(x=annot, fill=gene_biotype))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "purple","#F79D5C","darkgrey","#297373" ))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  scale_x_discrete(labels = c("gencode" = "GENCODEv47", "poder"="PODER"))+
  labs(x="", y="# ASTS Tested Genes", fill="Gene Biotype")

# which genes are tested with PODER and GENCODE
# astsgenesannot[, tested:="Tested"]

astsgenesannotwide <- dcast(unique(astsgenes[annot%in%c("GENCODEv47", "PODER"), .(annot, geneid.v, gene_biotype, gene_testable)]), ...~annot, fill="Not Tested")
astsgenesannotwide[, wheretested := case_when(
  GENCODEv47 == TRUE & PODER == TRUE ~ "Both",
  is.na(GENCODEv47) & PODER == TRUE ~ "Only PODER",
  GENCODEv47 == TRUE & is.na(PODER) ~ "Only GENCODEv47"
)]
astsgenesannotwide[, `:=`(wheretested=factor(wheretested, levels=c("Only GENCODEv47", "Both", "Only PODER")), 
                          gene_biotype=factor(gene_biotype, levels=c("protein_coding", "lncRNA", "novel/ambiguous gene", "IG_C_gene")))]


astsgenesannotwideplussamples <- astsgenes[, .(sample, population, map_reads_assemblymap, geneid.v, annot)][astsgenesannotwide, on="geneid.v"]
ggplot(unique(astsgenesannotwideplussamples[, .(wheretested, geneid.v, gene_biotype, population)]), aes(x=gene_biotype, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  scale_x_discrete(labels = c("gencode" = "GENCODEv47", "poder"="PODER"))+
  labs(x="", y="# ASTS Tested Genes", fill="Tested in")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  facet_wrap(~population)
ggplot(unique(astsgenesannotwideplussamples[, .(wheretested, geneid.v, gene_biotype)]), aes(x=gene_biotype, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  scale_x_discrete(labels = c("gencode" = "GENCODEv47", "poder"="PODER"))+
  labs(x="", y="# ASTS Tested Genes", fill="Tested in")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))



# # are we testing (and finding sig) more novel genes than expected by whats in PODER?
# chiinput <-rbind.data.frame(summary(factor(unique(poder[, .(geneid.v, gene_biotype)])$gene_biotype)),
#                             summary(factor(unique(astsgenesannot[annot=="poder", .(geneid.v, gene_biotype)])$gene_biotype)))
# colnames(chiinput) <- c("lncRNA", "novel/ambiguous", "protein_coding")
# rownames(chiinput) <- c("background", "tested")
# chiinput$novellnc <- apply(chiinput[, c("lncRNA", "novel/ambiguous")], 1, sum)
# 
# fisher.test(chiinput[, c("novellnc", "protein_coding")])
# mosaicplot(chiinput[, c("novellnc", "protein_coding")])
# mosaicplot(chiinput[, c("protein_coding", "lncRNA","novel/ambiguous")])

# intersection of tested genes across samples 
astsgenesannotwideplussamples[, samplesharing := uniqueN(sample), by=c("geneid.v","wheretested")]
astsgenesannotwideplussamples[, samplesharingbin :=cut(samplesharing,
                                                       breaks = c(0,  1, 3, 10, Inf),
                                                       labels = c("1", "2-3", "4-10", ">10"),
                                                       right = TRUE), by="wheretested"]
ggplot(unique(astsgenesannotwideplussamples[, .(samplesharingbin, wheretested, geneid.v)]), aes(x=samplesharingbin, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  labs(x="# Samples", y="# Genes Tested", fill="Tested in")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))

ggplot(unique(astsgenesannotwideplussamples[, .(samplesharing, wheretested, geneid.v)]), aes(x=samplesharing, fill=wheretested))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B","#D2D4C8", "#A2AD59"))+
  labs(x="# Samples", y="# Genes Tested", fill="Tested in")
# intersection of tested genes across populations

# What do the annotation specific genes have in particular so they are tested in higher proportion?
astsraw[, gene_morethanonetrx := gene_numtrx>1]

ggplot(unique(astsraw[, .(annot, trx_tencounts, transcript_variant, sample)]), aes(x=annot, fill=trx_tencounts))+
  geom_bar()+
  mytheme+
  labs(x="", y="# Transcript-Variant pairs", fill=">=10 counts")+
  scale_fill_manual(values=rev(c("#3C7A89", "#A84944")))+
  geom_text(aes(label=scales::percent(after_stat(count) / tapply(after_stat(count) , ..x.., sum)[..x..], accuracy = 1)),
            stat="count", position=position_stack(vjust=0.5))

ggplot(unique(astsraw[!is.na(gene_morethanonetrx), .(annot, gene_morethanonetrx, gene_variant,sample)]), aes(x=annot, fill=gene_morethanonetrx))+
  geom_bar(position="stack")+
  mytheme+
  labs(x="", y="# Gene-Variant pairs", fill=">1 Transcript/Gene")+
  scale_fill_manual(values=rev(c("#3C7A89", "#A84944")))+
  geom_text(aes(label=scales::percent(after_stat(count) / tapply(after_stat(count) , ..x.., sum)[..x..], accuracy = 1)),
            stat="count", position=position_stack(vjust=0.5))
ggplot(unique(astsraw[!is.na(gene_twentycount), .(annot, gene_twentycount, gene_variant, sample)]), aes(x=annot, fill=gene_twentycount))+
  geom_bar()+
  mytheme+
  labs(x="", y="# Gene-Variant pairs", fill=">=20 Counts/Gene")+
  scale_fill_manual(values=rev(c("#3C7A89", "#A84944")))+
  geom_text(aes(label=scales::percent(after_stat(count) / tapply(after_stat(count) , ..x.., sum)[..x..], accuracy = 1)),
            stat="count", position=position_stack(vjust=0.5))


ggplot(unique(astsraw[!is.na(gene_heterozygous), .(annot, gene_heterozygous, gene_variant, sample)]), aes(x=annot, fill=gene_heterozygous))+
  geom_bar()+
  mytheme+
  labs(x="", y="Proportion Gene-Variant pairs", fill="Both Alleles\nExpression")+
  scale_fill_manual(values=rev(c("#3C7A89", "#A84944")))+
  geom_text(aes(label=scales::percent(after_stat(count) / tapply(after_stat(count) , ..x.., sum)[..x..], accuracy = 1)),
            stat="count", position=position_stack(vjust=0.5))
#### IS PODER TESTING MORE GENES; TRANSCRIPTS; AND TRANSCRIPTS PER GENE??
ggplot(unique(unique(astsraw[gene_testable==TRUE, 
                             .(annot, gene_testable, geneid.v, sample)]
                     [, unique_tested_genes_per_sample:=uniqueN(geneid.v), by=c("sample", "annot")]
                     [, .(unique_tested_genes_per_sample, annot, sample)])), aes(x=annot, y=unique_tested_genes_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Tested Genes/Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(alternative = "less"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=550), vjust=0.5)

ggplot(unique(unique(astsraw[gene_testable==TRUE, 
                             .(annot, gene_testable, transcriptid.v, sample)]
                     [, unique_tested_transcripts_per_sample:=uniqueN(transcriptid.v), by=c("sample", "annot")]
                     [, .(unique_tested_transcripts_per_sample, annot, sample)])), aes(x=annot, y=unique_tested_transcripts_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Tested Transcripts/Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(alternative = "less"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=2000), vjust=0.5)

ggplot(unique(unique(astsraw[gene_testable==TRUE, 
                             .(annot, gene_testable, transcriptid.v, geneid.v, sample)]
                     [, unique_tested_transcripts_per_sample_per_gene:=uniqueN(transcriptid.v), by=c("sample", "geneid.v","annot")]
                     [, .(unique_tested_transcripts_per_sample_per_gene, annot, sample, geneid.v)])), aes(x=annot, y=unique_tested_transcripts_per_sample_per_gene, fill=annot))+
  geom_violin(alpha=0.7, adjust=3)+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  labs(x="", y="# Unique Tested Transcripts/Gene per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(alternative = "less"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=45), vjust=0.5)

### IS PODER DETECTING MORE TRANSCRIPTS; GENES AND TRANSCRPT/GENE WITH ALLELIC COUNTS;
ggplot(unique(unique(astsraw[, 
                             .(annot, gene_testable, geneid.v, sample)]
                     [, unique_tested_genes_per_sample:=uniqueN(geneid.v), by=c("sample", "annot")]
                     [, .(unique_tested_genes_per_sample, annot, sample)])), aes(x=annot, y=unique_tested_genes_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Detected Genes/Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args=list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=2500), vjust=0.5)

ggplot(unique(unique(astsraw[, 
                             .(annot, gene_testable, transcriptid.v, sample)]
                     [, unique_tested_transcripts_per_sample:=uniqueN(transcriptid.v), by=c("sample", "annot")]
                     [, .(unique_tested_transcripts_per_sample, annot, sample)])), aes(x=annot, y=unique_tested_transcripts_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Detected Transcripts/Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args=list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=9000), vjust=0.5)

ggplot(unique(unique(astsraw[, 
                             .(annot, gene_testable, transcriptid.v, geneid.v, sample)]
                     [, unique_tested_transcripts_per_sample_per_gene:=uniqueN(transcriptid.v), by=c("sample", "geneid.v","annot")]
                     [, .(unique_tested_transcripts_per_sample_per_gene, annot, sample, geneid.v)])), aes(x=annot, y=unique_tested_transcripts_per_sample_per_gene, fill=annot))+
  geom_violin(alpha=0.7, adjust=3)+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  labs(x="", y="# Unique Detected Transcripts/Gene per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(alternative = "less"))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=80), vjust=0.5)


#### CHECK FILTERS
ggplot(unique(unique(astsraw[, 
                             .(annot, trx_tencounts, transcriptid.v, sample,population)]
                     [, unique_trx_tencounts_per_sample:=uniqueN(transcriptid.v), by=c("sample","annot", "trx_tencounts")]
                     [, .(unique_trx_tencounts_per_sample, annot, sample, trx_tencounts, population)])), aes(x=annot, y=unique_trx_tencounts_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Transcripts per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)+
  facet_wrap(~trx_tencounts, labeller=labeller(trx_tencounts=(c("FALSE"="<10 counts", "TRUE"=">=10 counts"))))



toplot1 <-unique(unique(astsraw[, 
                                .(annot, gene_numtrx, geneid.v, sample,population)][, gene_morethanonetrx:=gene_numtrx>1])
                 [, unique_gene_with2trx_per_sample:=uniqueN(geneid.v), by=c("sample","annot", "gene_morethanonetrx")]
                 [, .(unique_gene_with2trx_per_sample, annot, sample, gene_morethanonetrx, population)])
toplot1$gene_morethanonetrx[is.na(toplot1$gene_morethanonetrx)] <- "NA"   
toplot1[, gene_morethanonetrx:=factor(gene_morethanonetrx, levels=c("NA", "FALSE", "TRUE"))]
ggplot(toplot1, aes(x=annot, y=unique_gene_with2trx_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Genes per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)+
  facet_wrap(~gene_morethanonetrx, 
             labeller = labeller(gene_morethanonetrx = as_labeller(c("FALSE" = "1 transcript", 
                                                                     "TRUE" = ">=2 transcripts",
                                                                     "NA" = "No transcripts >=10 counts"))))

toplot1 <-unique(unique(astsraw[,.(annot, gene_numtrx, gene_twentycount, geneid.v, sample,population)]
                        [, gene_morethanonetrx:=gene_numtrx>1])
                 [gene_morethanonetrx==TRUE & gene_twentycount==TRUE, gene_enoughcounts_and_trx:=TRUE]
                 [, .(gene_enoughcounts_and_trx, annot, geneid.v, sample)])[,gene_enoughcounts_and_trx_count:=uniqueN(geneid.v), by=c("sample", "annot", "gene_enoughcounts_and_trx") ]

toplot1$gene_enoughcounts_and_trx[is.na(toplot1$gene_enoughcounts_and_trx)] <- "NA"   
toplot1[, gene_enoughcounts_and_trx:=factor(gene_enoughcounts_and_trx, levels=c("NA", "TRUE"))]
ggplot(unique(toplot1[, .(annot, gene_enoughcounts_and_trx_count, sample, gene_enoughcounts_and_trx)]), aes(x=annot, y=gene_enoughcounts_and_trx_count, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Genes per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)+
  facet_wrap(~gene_enoughcounts_and_trx, 
             labeller = labeller(gene_enoughcounts_and_trx = as_labeller(c("TRUE" = ">=20 gene counts\n&>1 transcripts >=10 counts",
                                                                     "NA" = "Not enough transcripts or counts"))))



toplot1 <-unique(unique(astsraw[, 
                                .(annot, gene_heterozygous, geneid.v, sample,population)])
                 [, unique_gene_heterozygous_per_sample:=uniqueN(geneid.v), by=c("sample","annot", "gene_heterozygous")]
                 [, .(unique_gene_heterozygous_per_sample, annot, sample, gene_heterozygous, population)])
toplot1$gene_heterozygous[is.na(toplot1$gene_heterozygous)] <- "NA"   
toplot1[, gene_heterozygous:=factor(gene_heterozygous, levels=c("NA", "FALSE", "TRUE"))]
ggplot(toplot1, aes(x=annot, y=unique_gene_heterozygous_per_sample, fill=annot))+
  geom_violin(alpha=0.7)+
  geom_boxplot(outliers=F, width=0.1)+
  geom_jitter()+
  mytheme+
  labs(x="", y="# Unique Genes per Sample")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("GENCODEv47", "PODER")),method = "t.test",
                     method.args = list(paired=T))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=0), vjust=0.5)+
  facet_wrap(~gene_heterozygous, 
             labeller = labeller(gene_heterozygous = as_labeller(c("FALSE" = "Unreliable heterozygous variants", 
                                                                  "TRUE" = "At least 1 heterozygous variant",
                                                                  "NA" = "No transcripts >=10 counts"))))






# number of detected transcripts
astsraw[gene_testable==TRUE, trx_tested_per_gene:=uniqueN(transcriptid.v), by=c("gene_variant", "sample", "annot")]
astsraw[, trx_unfilt_per_gene:=uniqueN(transcriptid.v), by=c("gene_variant", "sample", "annot")]

astsrawsublong <-melt(unique(astsraw[, .(gene_variant, sample, annot, gene_testable, trx_tested_per_gene, trx_unfilt_per_gene)]), 
                      value.name = "trxperGene", 
                      variable.name = "set",
                      measure.vars=c("trx_tested_per_gene", "trx_unfilt_per_gene"))

ggplot(astsrawsublong[!is.na(trxperGene) & sample=="YRI1"], aes(x=trxperGene,  fill=annot))+
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"))+
  mytheme+
  facet_wrap(~set, scales="free")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="# Expressed Transcripts/Gene", y="# Genes")

ggplot(astsrawsublong[!is.na(trxperGene) & sample=="YRI1"], aes(x=trxperGene,  fill=annot))+
  geom_area(stat="count")+
  mytheme+
  facet_wrap(~set, scales="free")+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="# Expressed Transcripts/Gene", y="# Genes")

astsraw[, detectedtrxmorethanone := uniqueN(transcriptid.v)>1, by=c("sample", "gene_variant", "annot")]
newlab <- c("Single Transcript Gene", "Multiple Transcripts Gene")
names(newlab) <- c("FALSE", "TRUE")

ggplot(astsraw[detectedtrxmorethanone==TRUE, minTrxCountperGene:=min(Count), by=c("sample", "gene_variant", "annot")], 
       aes(x=minTrxCountperGene,  fill=annot))+
  geom_density(alpha=0.5)+
  mytheme+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="# Min Transcript Count/Gene")+
  scale_x_continuous(trans="log10")


ggplot(astsraw, aes(x=trx_tencounts, fill=annot))+
  geom_bar()+
  mytheme+
  facet_wrap(~detectedtrxmorethanone, labeller=labeller(detectedtrxmorethanone=newlab))+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x=">=10 Transcript Counts", y="# Transcript-Variants")



ggplot(unique(astsraw[, .(transcript_variant, annot, sample)]), aes(x=annot, fill=annot))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="", y="# Unique Transcripts with Allelic Counts")+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-0.5)
ggplot(unique(astsraw[, .(transcript_variant, annot, sample, trx_tencounts)]), aes(x=annot, fill=trx_tencounts))+
  geom_bar()+
  mytheme+
  labs(x="", y="# Unique Transcripts with Allelic Counts")+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-0.5)
ggplot(unique(astsraw[, .(geneid.v, annot, sample)]), aes(x=annot, fill=annot))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="", y="# Unique Genes with Allelic Counts")+
  guides(fill="none")+
  geom_text(aes(label=after_stat(count)), stat="count", vjust=-0.5)


ggplot(astsraw[, trxperGenenew := uniqueN(transcriptid.v), by=c("annot", "geneid.v")], aes(x=annot, y=trxperGenenew, fill=annot))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.05)+
  mytheme+
  scale_fill_manual(values=c( "#7F675B", "#A2AD59"))+
  labs(x="", y="# Transcripts/Gene (in each sample)")+
  guides(fill="none")+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=150), vjust=0.5)





### explore the transcripts that have more than 10 counts and later those than have more than 10 counts and have been tested
astsraw[, trx_tencounts_tested:= paste(trx_tencounts, gene_testable, sep="_")]

# load transcripts data
master <- fread("../04_transcriptome_assembly/04_evaluation/05_mastertable/data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv")
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(master$structural_category)[c(1,2,3,5,4,8,6,7)]
submaster <-unique(master[filter=="pass", .(isoform, structural_category)])
submaster[, `:=`(transcriptid.v=isoform, isoform=NULL, annot="PODER")]

astsmod <- submaster[astsraw, on=c("transcriptid.v", "annot")]
astsmod[annot=="GENCODEv47", structural_category:="FSM"]
astsmod[is.na(structural_category), structural_category:="FSM"]


astsmod[, unique_trx_tencounts_tested:=uniqueN(transcriptid.v), by=c("sample", "annot", "structural_category", "trx_tencounts_tested")]
ggplot(unique(astsmod[structural_category %in% c("FSM", "NIC", "NNC"), 
                      .(annot, sample, unique_trx_tencounts_tested, structural_category, trx_tencounts_tested)]), 
       aes(x = annot, y = unique_trx_tencounts_tested, fill = structural_category)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.5, scale="width") +  # Transparent violin plot for better overlap
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Narrow boxplot to overlap
  mytheme +
  scale_fill_manual(values = colsqanti) +
  facet_wrap(~trx_tencounts_tested, labeller=labeller(trx_tencounts_tested=c("FALSE_NA"="Transcripts\n<10 counts",
                                                                             "TRUE_NA"="Transcripts\n>=10 counts\nNot Tested Genes",
                                                                             "TRUE_TRUE"="Tested Transcripts")))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=3750), vjust=0.5)+
  labs(x="", fill="Structural Category", y="# Transcripts")
ggplot(unique(astsmod[structural_category %in% c("FSM", "NIC", "NNC") & trx_tencounts_tested=="TRUE_TRUE", 
                      .(annot, sample, unique_trx_tencounts_tested, structural_category, trx_tencounts_tested, population)]), 
       aes(x = population, y = unique_trx_tencounts_tested, fill = structural_category)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.5, scale="width") +  # Transparent violin plot for better overlap
  geom_boxplot(aes(col=structural_category),width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Narrow boxplot to overlap
  mytheme +
  scale_fill_manual(values = colsqanti) +
  scale_color_manual(values = colsqanti) +
  facet_wrap(~trx_tencounts_tested, labeller=labeller(trx_tencounts_tested=c("FALSE_NA"="Transcripts\n<10 counts",
                                                                             "TRUE_NA"="Transcripts\n>=10 counts\nNot Tested Genes",
                                                                             "TRUE_TRUE"="Tested Transcripts")))+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y=3750), vjust=0.5)+
  labs(x="", fill="Structural Category", y="# Transcripts")

ggplot(unique(astsmod[annot=="PODER" & gene_testable==TRUE, .(sample, transcriptid.v, population, structural_category)]),
       aes(x=population, fill=structural_category))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colsqanti)+
  mytheme+
  labs(x="", y="Proportion of tested Transcripts", fill="Structural\nCategory")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))

astsmod[gene_testable==TRUE, gene_with_noveltrx := 
     any(structural_category %in% c("NIC", "NNC")), 
   by = .("sample", "annot", "geneid.v")]
astsmod[, morethan2trx :=ifelse(gene_numtrx>2, ">2", "=2")]
astsmod[!is.na(gene_with_noveltrx) & annot=="PODER", numgenes_per2trx_per_noveltrx_sample :=uniqueN(geneid.v), by=c("morethan2trx", "annot", "sample", "gene_with_noveltrx")]

# start from a subset dt
asts <- unique(astsmod[annot=="PODER" & gene_testable==TRUE, .(sample, population, structural_category, geneid.v, transcriptid.v)])
asts[, gene_with_novel := 
         any(structural_category %in% c("NNC", "NIC")), 
       by = c("geneid.v", "sample")][, gene_with_novel:=ifelse(gene_with_novel==FALSE, "Genes without\nnovel Transcripts", "Genes with\nnovel Trasncripts")]
asts[, trxPerGene := uniqueN(transcriptid.v), by=c("sample", "geneid.v")][, twotrx := ifelse(trxPerGene>2, ">2 transcript/gene", "=2 transcripts/gene")]
asts <- unique(asts[, .(geneid.v, sample, population, gene_with_novel, twotrx)])
asts[, genes_per_cat := uniqueN(geneid.v), by=c("sample","gene_with_novel", "twotrx")]
asts[, propgenes_per_cat := genes_per_cat/sum(uniqueN(geneid.v)), by=c("sample")]
ggplot(unique(asts[, .(sample, propgenes_per_cat, population, gene_with_novel, twotrx)]), 
       aes(x=population, y=propgenes_per_cat, fill=population))+
  geom_violin(alpha=0.8, scale="width")+
  geom_jitter()+
  scale_fill_manual(values=popcols)+
  labs(x="", y="# Tested Genes", fill="Population")+
  mytheme+
  facet_grid(gene_with_novel~twotrx)








ggplot(unique(astsmod[!is.na(gene_with_noveltrx) & annot=="PODER", .(numgenes_per2trx_per_noveltrx_sample, morethan2trx,sample, population,gene_with_noveltrx)])[, numgenes_per2trx_per_noveltrx_sample_prop :=numgenes_per2trx_per_noveltrx_sample/sum(numgenes_per2trx_per_noveltrx_sample), by=c("sample")], 
       aes(x=morethan2trx,y=numgenes_per2trx_per_noveltrx_sample_prop, color=population))+
  facet_wrap(~gene_with_noveltrx, labeller=labeller(gene_with_noveltrx=c("FALSE"="Genes without novel transcripts", "TRUE"="Genes with novel transcripts")))+
  geom_jitter(position=position_dodge(width = 0.5))+
  geom_boxplot(width=0.1, position=position_dodge(width = 0.5))+
  scale_color_manual(values=popcols)+
  labs(x="# Tested Transcripts/Gene", y="Proportion of Tested Genes")+
  mytheme
ggplot(unique(astsmod[!is.na(gene_with_noveltrx) & annot=="PODER", .(numgenes_per2trx_per_noveltrx_sample, structural_category,morethan2trx,sample, population,gene_with_noveltrx)])[, numgenes_per2trx_per_noveltrx_sample_prop :=numgenes_per2trx_per_noveltrx_sample/sum(numgenes_per2trx_per_noveltrx_sample), by=c("sample")], 
       aes(x=morethan2trx,y=numgenes_per2trx_per_noveltrx_sample_prop, fill=structural_category))+
  facet_grid(gene_with_noveltrx~population, labeller=labeller(gene_with_noveltrx=c("FALSE"="Genes without novel transcripts", "TRUE"="Genes with novel transcripts")))+
  geom_col()+
  scale_fill_manual(values=colsqanti)+
  labs(x="# Tested Transcripts/Gene", y="Proportion of Tested Genes")+
  mytheme
ggplot(unique(astsmod[!is.na(gene_with_noveltrx) & annot=="PODER", .(morethan2trx, gene_with_noveltrx,numgenes_per2trx_per_noveltrx_sample, sample, population)]), aes(x=morethan2trx,y=numgenes_per2trx_per_noveltrx_sample, color=population))+
  facet_wrap(~gene_with_noveltrx, labeller=labeller(gene_with_noveltrx=c("FALSE"="No novel transcripts", "TRUE"="Gene with novel transcripts")))+
  geom_jitter(position=position_dodge(width = 0.5))+
  geom_boxplot(width=0.1, position=position_dodge(width = 0.5))+
  mytheme+
  scale_color_manual(values=popcols)+
  labs(x="# Tested Transcripts/Gene", y="# Tested Genes")

ggplot(unique(astsmod[structural_category %in% c("FSM", "NIC", "NNC"), 
                      .(annot, structural_category,trx_tencounts_tested,transcriptid.v)][, structural_category:=factor(structural_category, levels=c("NNC", "NIC", "FSM"))]), 
       aes(x = annot, fill = structural_category))+
  geom_bar()+
  facet_wrap(~trx_tencounts_tested, labeller=labeller(trx_tencounts_tested=c("FALSE_NA"="Transcripts\n<10 counts",
                                                                             "TRUE_NA"="Transcripts\n>=10 counts\nNot Tested Genes",
                                                                             "TRUE_TRUE"="Tested Transcripts")))+
  mytheme +
  scale_fill_manual(values = colsqanti)+
  labs(x="", y="# Unique Transcripts across samples", fill="Structural\nCategory")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))
asts_tested_noveltrx <-unique(astsmod[gene_testable==TRUE & structural_category%in%c("NIC", "NNC") , transcriptid.v])




astsnoveltested <-melt(master[isoform%in%asts_tested_noveltrx], 
     measure.vars=c("CEU", "AJI", "YRI", "LWK", "MPC", "ITU", "HAC", "PEL"),
     value.name="detected",
     variable.name="population")[detected>=1,]
ggplot(astsnoveltested[!population%in%c("AJI", "MPC")],
       aes(x=population, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values = colsqanti) +
  labs(x="Discovered in", y="# ASTS tested novel transcripts", fill="Structural\nCategory")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))


ggplot(master[isoform%in%asts_tested_noveltrx, ceu_detected:=CEU>0][!is.na(ceu_detected)],
       aes(x=ceu_detected, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values = colsqanti) +
  labs(x="CEU discovered", y="# ASTS tested novel transcripts", fill="Structural\nCategory")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))
ggplot(master[isoform%in%asts_tested_noveltrx, ],
       aes(x=ceu_detected, fill=structural_category))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values = colsqanti) +
  labs(x="CEU discovered", y="# ASTS tested novel transcripts")
astsrawwide <- dcast(astsraw, ... ~ annot, 
                     value.var = c("trx_tested_per_gene", "trx_unfilt_per_gene"))
subastsrawide <- astsrawwide[, .(geneid.v, gene_biotype, population, sample, 
                                 map_reads_assemblymap, gene_testable, 
                                 trx_tested_per_gene_GENCODEv47, trx_tested_per_gene_PODER, 
                                 trx_unfilt_per_gene_GENCODEv47, trx_unfilt_per_gene_PODER)]
subastsrawide[, `:=`(trxTestedPerGene_diff=trx_tested_per_gene_PODER-trx_tested_per_gene_GENCODEv47,
                     trxUnfiltPerGene_diff=trx_unfilt_per_gene_PODER-trx_unfilt_per_gene_GENCODEv47)]
ggplot(subastsrawide, aes(x=population, y=trxTestedPerGene_diff, col=population))+
  geom_jitter()+
  mytheme+
  scale_color_manual(values=popcols)

# number of tested and significant variables
variantsub <-unique(asts[, .(variant, sample, population, FDR, annot)])[, sig:=FDR<0.05][, num_var := uniqueN(variant), 
                                                           by=c("sample", "annot", )]
ggplot()




























###################### CHECK THE SAME FOR ENHANCED GENCODE THAN FOR PODER
# Explore how many more genes are tested in each annot
astswide <- dcast(unique(asts[, .(sample, annot, tested_genes)]), sample ~ annot, value.var = "tested_genes") 
astswide[, diff_tested_genes:=`Enhanced GENCODEv47`-GENCODEv47]
astswide[, ratio_tested_genes:=`Enhanced GENCODEv47`/GENCODEv47]

astswidemeta <- metadata[, .(sample, population, map_reads_assemblymap)][astswide, on="sample"]
astswidemeta[, eur:=ifelse(population %in% c("LWK", "YRI"), "AFR",
                           ifelse(population %in% c("PEL", "HAC", "ITU"), "nonEUR-OOA", "EUR"))]
ggplot(astswidemeta, aes(x=sample, y=diff_tested_genes, fill=population))+
  geom_col()+
  mytheme+
  scale_fill_manual(values=popcols)
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=diff_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Increment of genes tested by\nEnhanced GENCODE against GENCODE", color="Population")+
  ylim(c(0,115))
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=ratio_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Ratio of genes tested by PODER against GENCODE")+
  geom_hline(yintercept=1, linetype="dashed", color="black")

# are EUR beign the least benefitted?
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=diff_tested_genes))+
  mytheme+
  stat_poly_line(color="darkgrey")+
  stat_poly_eq(use_label(c( "p", "n")),label.y=0.9)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95, size = 5) +  # Equation and R-squared
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Increment of genes tested by PODER against GENCODE")+
  ylim(c(0,115))+
  facet_wrap(~eur)
ggplot(astswidemeta, aes(x=map_reads_assemblymap/10^6, y=ratio_tested_genes))+
  mytheme+
  geom_point(alpha=0.8, size=3, aes(col=population))+
  scale_color_manual(values=popcols)+
  labs(x="Mapped Reads (M)", y="Ratio of genes tested by PODER against GENCODE")+
  geom_hline(yintercept=1, linetype="dashed", color="black")+
  facet_wrap(~population)

# check  normality and homocedasticity
car::leveneTest(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
bartlett.test(ratio_tested_genes ~ population, data = astswidemeta[, .(sample, population, ratio_tested_genes)])
shapiro.test(astswidemeta$ratio_tested_genes)
qqnorm(astswidemeta$ratio_tested_genes)
qqline(astswidemeta$ratio_tested_genes, col = "red")
my_comparisons <- list( c("CEU", "HAC"), c("CEU", "ITU"), c("CEU", "LWK"),c("CEU", "PEL"), c("CEU", "YRI") )
ggplot(astswidemeta, aes(x=population, y=ratio_tested_genes))+
  mytheme+
  geom_boxplot(outliers=F)+
  ggbeeswarm::geom_quasirandom(alpha=0.8, aes(col=population, size=map_reads_assemblymap/10^6))+
  scale_color_manual(values=popcols)+
  labs(x="", y="Ratio of genes tested by\nEnhanced GENCODE against GENCODE", size="Mapped Reads (M)", col="Population")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                     method.args = list(alternative = "less")) # Add pairwise comparisons p-value

ggplot(astswidemeta, aes(x=population, y=map_reads_assemblymap/10^6))+
  mytheme+
  geom_boxplot(outliers=F)+
  geom_jitter(alpha=0.8, aes(col=population, size=ratio_tested_genes))+
  scale_color_manual(values=popcols)+
  labs(x="", y="Mapped reads (M)", size="Ratio tested genes\nPODER/GENCODE", col="Population")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",
                     method.args = list(alternative = "less")) 
