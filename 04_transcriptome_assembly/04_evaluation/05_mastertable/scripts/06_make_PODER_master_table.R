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



# load data
mastertable <- fread("data/240909merge_reslongmeta_annot_sqanti_sj_gencodev47_quantification_novellocus_proteinInfo_updatedrecount_disambiguatedGenes_replacedFLAIRna&addedISMinfo_filteredINFO.tsv")
sqanti <- fread("../02_sqanti/data/poder/poder_classification.txt")
protein <- fread("../../../../novelannotations/poder_protein/transcripts_protein_annotation.tsv")
annot <- fread("../../../../../../Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
# change levels of structuralc categories
categories <- c("FSM",  "NNC","Fusion","NIC","Antisense","Intergenic","Genic")
names(categories) <- unique(sqanti$structural_category)
sqanti$structural_category <- categories[sqanti$structural_category]
sqanti[, structural_category :=factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense"))]
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- levels(sqanti$structural_category)
sqanti[, annotated:=ifelse(grepl("ENST", isoform), "annotated", "discovered")]


# keep only poder transcripts
mastertable <- mastertable[filter=="pass"]

# prepare new protein table to be included in the master table
protein[, isoform:=tid]
colnames(protein) <- paste0("proteinv47_", colnames(protein))
colnames(protein)[grepl("isoform", colnames(protein))] <- "isoform"


# now separate ISMs
ism <- mastertable[structural_category=="ISM"]

# ggplot(ism[, .(multiple_subchains=uniqueN(isoform), subcategory, isoform), by=c("associated_transcriptid.v")][, .(subcategory, isoform,multiple=ifelse(multiple_subchains>1, "Multiple", "Unique"))], 
#        aes(x=subcategory, fill=multiple))+
#   geom_bar()+
#   mytheme+
#   scale_x_discrete(labels=c("3prime_fragment"="3' Fragment",
#                             "5prime_fragment"="5' Fragment",
#                             "internal_fragment"="Internal Fragment",
#                             "intron_retention"="Intron Retention"))+
#   labs(x="Structural Subcategory", y="# ISM replaced", fill="# Subchains of Associated\nAnnotated Transcript")+
#   geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
#   theme(legend.position= c(0.8, 0.9))+
#   scale_fill_manual(values=rev(c("#BBBAC6", "#8F7E4F")))
# ggsave("../../../10_figures/suppfig/barplot_ISM_PerSqantSubiCategory_SubchainsPerTranscript.pdf", dpi=700, width = 16, height = 14,  units = "cm")

# summarize ISM subchains
ism_sub <- ism[,.SD, .SDcols = c(which(colnames(ism)%in%c("associated_transcriptid.v","espresso","flair","isoquant","lyric")),grep("isoform", colnames(ism)),grep("^[A-Z]{3}",colnames(ism)))]
numeric_cols <- names(ism_sub)[sapply(ism_sub, is.numeric)]
ism_sum <-ism_sub[, lapply(.SD, sum), by = "associated_transcriptid.v",.SDcols = numeric_cols]
# now replace all >1 by 1
# Apply the replacement only to columns that match the regex pattern "^[A-Z]{3}[0-9]$"
ism_sum[, (grep("^[A-Z]{3}[0-9]$", names(ism_sum), value = TRUE)) := 
          lapply(.SD, function(x) if (is.numeric(x)) fifelse(x > 1, 1, x) else x), 
        .SDcols = grep("^[A-Z]{3}[0-9]$", names(ism_sum), value = TRUE)]
colnames(ism_sum)[grep("associated_transcriptid.v", colnames(ism_sum))] <- "isoform"

# add sharings
# Population sharing
population_cols <- grep("^[A-Z]{3}$", names(ism_sum), value = TRUE)
ism_sum[, population_sharing := rowSums(.SD), .SDcols = population_cols]
# Sample sharing
sample_cols <- grep("^[A-Z]{3}[0-9]{1}$", names(ism_sum), value = TRUE)
ism_sum[, sample_sharing := rowSums(.SD), .SDcols = sample_cols]
# Tool sharing
tool_cols <-c("espresso", "lyric", "flair", "isoquant")
ism_sum[, tool_sharing := rowSums(.SD), .SDcols = tool_cols]
mastertable_nonism <- mastertable[structural_category!="ISM", .SD, .SDcols = c(which(colnames(mastertable)%in%c("isoform", "espresso", "flair", "isoquant", "lyric", "population_sharing","tool_sharing", "sample_sharing")),
                                                                                     grep("^[A-Z]{3}", colnames(mastertable)))]

full_sharing <-rbind.data.frame(ism_sum, mastertable_nonism)
sqanti <-sqanti[, Filter(function(x) !all(is.na(x)), .SD)]

final <- full_sharing[sqanti, on="isoform"]

final <- protein[final, on="isoform"]

final[, predicted_ORF:=fifelse(ORF_length>0, "Predicted ORF", "No ORF", "No ORF")]
final[, proteinv47_predicted_ORF:=fifelse(proteinv47_orf_length_nt>0, "Predicted ORF", "No ORF", "No ORF")]

final <- unique(mastertable[, .(isoform, associated_gene_biotype)])[final, on="isoform"]

utrx <- unique(mastertable[, .(associated_gene_biotype,associated_transcriptid.v)])
biotypesvec <-utrx[, associated_gene_biotype]
names(biotypesvec) <- utrx[, associated_transcriptid.v]
final[, associated_gene_biotype := fifelse(is.na(associated_gene_biotype), 
                                           as.character(biotypesvec[isoform]), 
                                           as.character(associated_gene_biotype))]
final[, associated_gene_biotype:=factor(fifelse(associated_gene_biotype=="protein_coding", "Protein Coding", 
                                                fifelse(associated_gene_biotype=="", "Novel/Ambiguous Gene", associated_gene_biotype)), levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene"))]

# load gene transcript match
gene_trx_correspondence <- fread("../novelannotations/merged/240926_filtered_with_genes.transcript2gene.tsv", header = F)
colnames(gene_trx_correspondence) <- c("isoform", "geneid.v")
final <- gene_trx_correspondence[data, on="isoform"]

fwrite(final,"data/29102024_PODER_mastertable.tsv", row.names = F, quote = F, sep="\t")
