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

popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)


nametype <- c("gencode", "poder", "enhanced_gencode")
# load data
ase <-list(fread("data/ASE_results_gencode.tsv"))
ase <- append(ase,
              list(fread("data/ASE_results_pantrx.tsv")))
ase <- append(ase,
              list(fread("data/ASE_results_enhanced.tsv")))
asts <- list(fread("data/ASTS_results_gencode_variantgenefilt.tsv"))  
asts <- append(asts, 
               list(fread("data/ASTS_results_pantrx_variantgenefilt.tsv")))
asts <- append(asts, 
               list(fread("data/ASTS_results_enhanced_variantgenefilt.tsv")))

names(ase) <- nametype
names(asts) <- nametype

ase <- rbindlist(ase, idcol="annot")
asts <- rbindlist(asts, idcol="annot")

# load annotations
gencode <- fread("../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
poder <- fread("../../novelannotations/merged/240926_filtered_with_genes.transcript2gene_with_biotypes.tsv")
gencode[, annot:="gencode"]
poder[, annot:="poder"]
annot <- unique(rbind.data.frame(poder, gencode[, .(geneid.v, transcriptid.v, gene_biotype, annot)]))


# Explore how many more genes are tested in each annot
astswide <- dcast(unique(asts[, .(sample, annot, tested_genes)]), sample ~ annot, value.var = "tested_genes") 
astswide[, diff_tested_genes:=poder-gencode]
astswide[, ratio_tested_genes:=poder/gencode]

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
ggplot(astswidemeta, aes(x=population, y=ratio_tested_genes))+
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
  geom_jitter(alpha=0.8, aes(col=population, size=ratio_tested_genes))+
  scale_color_manual(values=popcols)+
  labs(x="", y="Mapped reads (M)", size="Ratio tested genes\nPODER/GENCODE", col="Population")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",
                     method.args = list(alternative = "less")) 

# check biotypes of tested genes
astsgenes <- unique(asts[, .(annot, sample, population, map_reads_assemblymap, geneid.v)])
astsgenes <- astsgenes[geneid.v!=""]
astsgenesannot <- unique(annot[, .(geneid.v, gene_biotype, annot)])[astsgenes, on=c("geneid.v", "annot")]
astsgenesannot$gene_biotype <- ifelse(grepl("_", astsgenesannot$geneid.v), "novel/ambiguous gene", astsgenesannot$gene_biotype)

ggplot(unique(astsgenesannot[, .(annot, geneid.v, gene_biotype)]), aes(x=annot, fill=gene_biotype))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c( "purple","#F79D5C","darkgrey","#297373" ))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  scale_x_discrete(labels = c("gencode" = "GENCODEv47", "poder"="PODER"))+
  labs(x="", y="# ASTS Tested Genes", fill="Gene Biotype")

# which genes are tested with PODER and GENCODE
astsgenesannot[, tested:="Tested"]

astsgenesannotwide <- dcast(unique(astsgenesannot[, .(annot, geneid.v, gene_biotype, tested)]), ...~annot, fill="Not Tested")
astsgenesannotwide[, wheretested:=ifelse(gencode=="Tested" & poder=="Tested", "Both",
                                         ifelse(gencode=="Tested" & poder=="Not Tested", "Only GENCODEv47",
                                                ifelse(gencode=="Not Tested" & poder=="Tested", "Only PODER", "Error")))]
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

# are we testing (and finding sig) more novel genes than expected by whats in PODER?
chiinput <-rbind.data.frame(summary(factor(unique(poder[, .(geneid.v, gene_biotype)])$gene_biotype)),
                            summary(factor(unique(astsgenesannot[annot=="poder", .(geneid.v, gene_biotype)])$gene_biotype)))
colnames(chiinput) <- c("lncRNA", "novel/ambiguous", "protein_coding")
rownames(chiinput) <- c("background", "tested")
chiinput$novellnc <- apply(chiinput[, c("lncRNA", "novel/ambiguous")], 1, sum)

fisher.test(chiinput[, c("novellnc", "protein_coding")])
mosaicplot(chiinput[, c("novellnc", "protein_coding")])
mosaicplot(chiinput[, c("protein_coding", "lncRNA","novel/ambiguous")])

# intersection of tested genes across samples AIXO ESTA AMLAMENT
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
# load data of calc asts
arraygen <- fread("array_gencode", header = F)
arraypan <- fread("array_pantrx", header = F)
array <- rbind.data.frame(arraygen, arraypan)

rawasts <- list()

for(i in 1:nrow(array)){
  TYPE<-array[i, V3]
  SAMPLE<- array[i, V1]
  print(TYPE)
  print(SAMPLE)
  temp <- fread(paste0("data/", TYPE,"/05_calc_asts/", SAMPLE,"_asts_quant.tsv"))
  temp[, `:=`(annot=TYPE, samplecode=SAMPLE)]
  rawasts <- append(rawasts, list(temp))}

astsraw <- rbindlist(rawasts)

# prepare raw data
astsraw <- annot[astsraw, on=c("transcriptid.v"="transcript", "annot")]
astsraw[, variant := paste(contig, position, refAllele, altAllele, sep="_")]   
astsraw[, `:=`(contig=NULL, position=NULL, refAllele=NULL, altAllele=NULL)]

# Transcripts with at least 10 counts
astsraw[, Count:=refCount+altCount][, trx_count_filter_pass:=Count>=10]
# Gene-variant pairs with less than 20 counts
astsraw[, morethan20GeneCount:=sum(Count)>=20, by=c("variant", "geneid.v", "samplecode", "annot")]

# Keep only genes with more than 1 transcript
astsraw[trx_count_filter_pass==TRUE, passtrxCountperGene := uniqueN(transcriptid.v) , by=c("variant","geneid.v","samplecode", "annot")]

# Remove genes where all the counts are on the same allele because they are most likely wrongly genotyped (or not trustable)
astsraw[trx_count_filter_pass==TRUE & passtrxCountperGene>1, trustable_genevariants:=(sum(refCount)!=0 & sum(Count)!=sum(refCount)), by=c("variant","geneid.v","samplecode", "annot")]


astsraw[passtrxCountperGene>1 & trustable_genevariants==TRUE , gene_filter_pass:=TRUE]
astsraw[is.na(gene_filter_pass), gene_filter_pass:=FALSE]

# Create gene-variant pairs to test all the variants in a gene
astsraw[, gene_variant:=paste(geneid.v, variant, sep=":")]
astsraw[, transcript_variant := paste(transcriptid.v, variant, sep=":")]
astsraw[, cell_line_id:=tstrsplit(samplecode, "_")[[3]]]
astsraw <-metadata[, .(cell_line_id, sample)][astsraw, on="cell_line_id"]
# # merge with previous data results
# astsfull <- unique(astsgenesannotwideplussamples)[astsraw , on=c("sample", "geneid.v")]
# 
# ggplot(unique(astsfull[gene_filter_pass==TRUE, .(sample, annot, geneid.v)]))
# # are we testing a different amount of gene-variant pairs?
# # how many transcripts/gene are we testing?



full <- unique(unique(astsraw[gene_filter_pass==TRUE, .(sample, annot, geneid.v)])$geneid.v)
part <- unique(unique(asts[, .(sample, annot, geneid.v)])$geneid.v)



discordant <- astsraw[geneid.v%in%(part[!part%in%full])]

discordant <- unique(discordant[trx_count_filter_pass==TRUE & gene_filter_pass==TRUE, .(sample, geneid.v, transcriptid.v, refCount, altCount, Count, annot, morethan20GeneCount, passtrxCountperGene, trustable_genevariants,gene_filter_pass )])


part[!part%in%full]
full[!full%in%part]













astsraw[, Count:=refCount+altCount][, trx_count_filter_pass:=Count>=10]
astsraw[, morethan20GeneCount:=sum(Count)>=20, by=c("variant", "geneid.v", "samplecode", "annot")]
astsraw[trx_count_filter_pass==TRUE, passtrxCountperGene := uniqueN(transcriptid.v) , by=c("variant","geneid.v","samplecode", "annot")]
astsraw[passtrxCountperGene>1 , gene_filter_pass:=TRUE]