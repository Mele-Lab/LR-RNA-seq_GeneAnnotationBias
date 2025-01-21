## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/06_quantification"
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

data <- fread("../04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- unique(data$structural_category)[c(1,4,2,6,7,3,5)]
data[, structural_category := factor(structural_category, levels=names(colsqanti))]
data[, associated_gene_biotype := factor(associated_gene_biotype, levels=c("Protein Coding", "lncRNA", "Novel/Ambiguous Gene"))]



metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE & mixed_samples==FALSE]
popcol <- metadata$color_pop
names(popcol) <- metadata$population

# ISOQUANT
iq <- fread("01_isoquantify/data/poder/concat_isoquant_counts.tsv")
colnames(iq) <- c("trx", "counts", "sample")
iq <- iq[!grepl("^_", trx)]
iq[, cell_line_id:=gsub(".*_", "", sample)][, sample:=NULL]
iq <- metadata[, .(cell_line_id, sample)][iq, on="cell_line_id"][, cell_line_id:=NULL][, tool:="IsoQuant"]
# FLAIR
fl <- fread("02_flairquantify/data/poder/concat_flair_counts.tsv")
colnames(fl) <- c("trx", "counts", "sample")
fl[, tool:="FLAIR"]
# KALLISTO
ka <- fread("../../novelannotations/quantifications/kallisto_quant/matrix.abundance.tsv")
ka <- melt(ka, id.vars = "transcript_id", value.name = "counts", variable.name = "cell_line_id")
ka[, cell_line_id:=gsub("_1", "", cell_line_id)]
ka <- metadata[, .(cell_line_id, sample)][ka, on="cell_line_id"][, cell_line_id:=NULL]
colnames(ka) <- c("sample", "trx", "counts")
ka[, tool:="Kallisto"]
# concatenate all of them
mix <- rbind.data.frame(iq, fl)
mix <- rbind.data.frame(mix, ka)
mix[, population:=gsub(".$", "", sample)]


# remove counts=0
mix <- mix[counts>0]
mix <- data[, .(isoform, structural_category)][mix, on=c("isoform"="trx")]
addpoder <- data[, .(isoform, structural_category)][, `:=`(sample="mock", counts=1, tool="PODER\nBackground", population="mock")]
mix <- rbind.data.frame(mix, addpoder)

ggplot(unique(mix[, .(tool, isoform, structural_category)]), aes(x=tool, fill=structural_category))+
  geom_bar()+
  mythemen+
  labs(x="", y="# Transcripts with counts>0", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  scale_fill_manual(values=colsqanti)

# do correlations
mix_wide <- dcast(mix[, .(sample, isoform, tool, counts)], sample+isoform~tool)
mix_wide[is.na(mix_wide)] <- 0
mix_wide[mix_wide==0] <- NA

cor(mix_wide$FLAIR, mix_wide$IsoQuant, use="complete", method = "spearman")
cor(mix_wide$Kallisto, mix_wide$IsoQuant,  use="complete", method = "spearman")
cor(mix_wide$Kallisto, mix_wide$FLAIR,  use="complete", method = "spearman")

ggplot(mix_wide, aes(x=FLAIR))+
  geom_point(aes(y=Kallisto), color="black")+
  geom_point(aes(y=IsoQuant), color="darkred")+
  mythemen+
  coord_equal()


summar <- unique(mix[tool!="Background\nPODER" & sample!="mock"][, summ:=sum(counts), by=c("tool", "sample")][, .(summ, tool, sample)])
summar <- metadata[, .(sample, map_reads_generalmap)][summar, on="sample"]

ggplot(summar, 
       aes(x=sample, fill=tool, y=summ))+
  geom_col(position="dodge")+
  geom_point(aes(x=sample, y=map_reads_generalmap))+
  mythemen+
  theme(axis.text.x = element_text(angle=90))+
  labs(x="", y="Total Counts = Reads Used")
