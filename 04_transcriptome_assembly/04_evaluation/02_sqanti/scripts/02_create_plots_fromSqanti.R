## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/02_sqanti"
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

require(MASS) # to access Animals data sets
require(scales)
# load data
# isoclass <- fread("data/unmasked_genomic/unmasked_genomic_isoquant_guided_classification.txt")
# isoclass <- fread("data/pseudomasked_genomic_espresso/pseudomasked_genomic_espresso_classification.txt")
isoclass <- fread("data/pseudomasked_genomic_isoquant_guided/pseudomasked_genomic_isoquant_guided_classification.txt")
isoclass <- isoclass[!is.na(isoclass$structural_category),]
# Original vector
lookup_vector <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", 
                   "novel_not_in_catalog", "genic", "antisense", "fusion","intergenic")
replacement_vector <- c("FSM", "ISM", "NIC", "NNC", "Genic", "Antisense", "Fusion","Intergenic")
replacement_lookup <- setNames(replacement_vector, lookup_vector)
# Change full-splice_match to FSM etc
isoclass$structural_category <- replacement_lookup[isoclass$structural_category]
isoclass[, structural_category := factor(structural_category, levels = c("FSM", 
                                                                         "ISM", 
                                                                         "NIC", 
                                                                         "NNC", 
                                                                         "Genic", 
                                                                         "Antisense",
                                                                         "Fusion",
                                                                         "Intergenic"))]


# count isoforms
isoformlvl <- unique(isoclass[, .(structural_category, isoform,subcategory)])[, isoformspercat := .N, by=structural_category]
isoformlvl <- unique(isoformlvl[, .(subcategory, isoform, structural_category, isoformspercat)])[, isoformspersubcat := .N, by=subcategory]

# count genes
genelvl <- unique(isoclass[, .(structural_category, associated_gene)])[, genespercat := .N, by=c("structural_category")]

# merge data
bothlvl <- unique(genelvl[, .(structural_category, genespercat)])[isoformlvl, on="structural_category"]


# mycols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# mycols <- c("#009E73", "#0072B2", "#F0E442", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#CC79A7")
# mycols <- c("#90BE6D", "#577590", "#F9C74F", "#F3722C", "#4D908E","#CC79A7", "#F8961E", "#F94144")
mycols <- c("#448AFF", "#0D53A3", "#009688", "#8BC34A", "#FFC107", "#FF9800", "#F44336", "#AD1457")

# Plot number of genes and isoforms
ggplot(unique(bothlvl[,.(isoformspercat, genespercat, structural_category)]),
                  aes(x=structural_category))+
  geom_col(aes(y=isoformspercat, fill=structural_category), alpha=0.5)+
  geom_col(aes(y=genespercat, fill=structural_category))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  annotation_logticks(sides="l")+
  mytheme+
  geom_text(aes(label=isoformspercat, y=isoformspercat), vjust=-0.5)+
  geom_text(aes(label=genespercat, y=genespercat), vjust=1.5)+
  expand_limits(y = max(bothlvl$isoformspercat) * 10)+
  scale_fill_manual(values =mycols)+
  labs(fill="", x="", y="# features")+
  guides(fill="none")



# plot #models / gene by sqanticat
modelspergene <- unique(isoclass[, .(isoform, associated_gene,structural_category)])[, modelspergene:=.N, by=c("associated_gene", "structural_category")]
modelspergene$modelspergenebin <- ifelse(modelspergene$modelspergene>=50, ">=50", 
                                      ifelse(modelspergene$modelspergene>=10, ">=10", modelspergene$modelspergene))

ggplot(modelspergene)+
  geom_histogram(aes(x=factor(modelspergenebin, levels=c(1:9, ">=10", ">=50")), fill=structural_category), stat="count")+
  facet_wrap(~structural_category, scales="free_y")+
  scale_fill_manual(values =mycols)+
  labs(fill="", x="# Models/gene", y="# genes")+
  guides(fill="none")+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  



# plot lengths of the isoforms
ggplot(isoclass, aes(y=length, fill=structural_category, x=structural_category))+
  geom_hline(yintercept = median(isoclass[structural_category=="FSM",length]),
             linetype="dashed", size=1.5, color="darkgrey")+
  geom_violin(alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  scale_fill_manual(values =mycols)+
  labs(fill="", x="", y="Isoform length")+
  guides(fill="none")+
  annotation_logticks(sides="l")

# plot lengths of the FSM isoforms by subcategory
ggplot(isoclass[structural_category=="FSM",], 
       aes(y=length, fill=subcategory, x=subcategory))+
  geom_hline(yintercept = median(isoclass[structural_category=="FSM",length]),
             linetype="dashed", size=1.5, color="darkgrey")+
  geom_violin(alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  scale_fill_manual(values =c("#6A994E", "#386641", "#A7C957", "#BC4749", "#F2E8CF"))+
  labs(fill="", x="", y="Isoform length")+
  guides(fill="none")+
  annotation_logticks(sides="l")

# load annot data
annot <- fread("../../../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
annotgene <- unique(annot[, .(gene_biotype, geneid.v, gene_name)])
annotiso <- unique(annot[,.(transcript_biotype, transcriptid.v)])
colnames(isoclass)[grep("associated_gene",colnames(isoclass))] <- "gene_name"
colnames(isoclass)[grep("associated_transcript",colnames(isoclass))] <- "transcriptid.v"


isoclassannot <- annotgene[isoclass, on= "gene_name"]
isoclassannot <- annotiso[isoclassannot, on= "transcriptid.v"]

isoclassannotgene <-unique(annotgene[isoclass, on= "geneid"][!is.na(gene_biotype)])
ggplot(isoclassannotgene, aes(x=structural_category, fill=gene_biotype))+
  geom_histogram(stat="count")

# plot length of isoform by associated-gene biotype and structural category
ggplot(isoclassannot[structural_category%in%c("FSM", "ISM", "NIC", "NNC") & !is.na(gene_biotype)], 
       aes(x=gene_biotype, y= length, fill=structural_category))+
  geom_hline(yintercept = median(isoclass[structural_category=="FSM",length]),
             linetype="dashed", size=1.5, color="darkgrey")+
  geom_jitter(size=0.5, alpha=0.25)+
  geom_violin(alpha=0.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  geom_boxplot(outliers=F, width=0.1, alpha=0.8)+
  mytheme+
  facet_wrap(~structural_category, scales = "free_x")+
  scale_fill_manual(values =mycols)+
  labs(fill="", x="", y="Isoform length")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  annotation_logticks(sides="l")

# plot length of isoform by protein-coding gene-associated models by transcript
# biotype and structural category

ggplot(isoclassannot[structural_category%in%c("FSM", "ISM") & gene_biotype=="protein_coding"], 
       aes(x=transcript_biotype, y= length, fill=structural_category))+
  geom_hline(yintercept = median(isoclassannot[structural_category=="FSM" & gene_biotype=="protein_coding",length]),
             linetype="dashed", size=1.5, color="darkgrey")+
  geom_jitter(size=0.5, alpha=0.4)+
  geom_violin(alpha=0.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  geom_boxplot(outliers=F, width=0.1)+
  mytheme+
  facet_wrap(~structural_category, scales = "free_x")+
  scale_fill_manual(values =mycols)+
  labs(fill="", x="", y="Isoform length of protein-coding\ngene associated models")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  annotation_logticks(sides="l")

# Plot number of isoforms by coding potential
ggplot(unique(isoclassannot[gene_biotype%in%c("protein_coding", NA)][, .(structural_category,  coding, isoform)]), aes(x=structural_category, alpha=coding, fill=structural_category))+
  geom_histogram(stat="count", position="dodge")+
  mytheme+
  scale_fill_manual(values =mycols)+
  scale_alpha_manual(values=c(1, 0.55))+
  labs(fill="", x="", y="# features")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0.1))+
  annotation_logticks(sides="l")+
  labs(fill="", x="", alpha="Predicted", y="# isoforms of novel genes\nor protein-coding genes associated")
# number of exons

# SJ plots

# load SJ info
#sj <- fread("data/pseudomasked_genomic_espresso/pseudomasked_genomic_espresso_junctions.txt", stringsAsFactors=T)
sj <- fread("data/pseudomasked_genomic_isoquant_guided/pseudomasked_genomic_isoquant_guided_junctions.txt", stringsAsFactors=T)

sj[, numbersj := .N, by=.(junction_category, canonical)][, percentsj := numbersj/nrow(sj)*100]

ggplot(unique(sj[, .(junction_category, numbersj, canonical, percentsj)]), aes(x=junction_category, y=numbersj, fill=junction_category, alpha=canonical))+
  geom_col(position="stack")+
  mytheme+
  scale_fill_manual(values =c("#94B84D", "#9239B3"))+
  scale_alpha_manual(values=c(1, 0.55))+
  labs(y="# SJ", x="SJ", alpha="Splice Sites")+
  geom_text(aes(label=sprintf("%.2f%%", percentsj)),  position = position_stack(vjust = 0.5))+
  guides(fill="none")

sj[, junction_number := as.integer(gsub("junction_", "",junction_number))]
ggplot(sj, aes(x=junction_number))+
  geom_histogram(position="dodge", binwidth = 1)+
  mytheme+
  labs(y="# SJ", )+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0.1))+
  facet_wrap(~canonical+junction_category)

ggplot(sj, aes(x=junction_number,fill=junction_category, alpha=canonical))+
  geom_histogram(position="dodge", binwidth = 1)+
  facet_wrap(~canonical+junction_category)+
  scale_fill_manual(values =c("#94B84D", "#9239B3"))+
  scale_alpha_manual(values=c(1, 0.55))+
  xlim(c(0,50))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand=c(0,0))+
  mytheme
# PSI plots
