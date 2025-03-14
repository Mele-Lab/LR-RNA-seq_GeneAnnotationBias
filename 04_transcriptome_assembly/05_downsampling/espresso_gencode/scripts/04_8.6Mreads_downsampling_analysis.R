## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling/espresso"
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


n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
metadata <- fread("../../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[merged_run_mode==TRUE]
metadata <- metadata[mixed_samples==FALSE]
popcols <- unique(metadata$color_pop)
names(popcols) <- unique(metadata$population)
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C")
names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC")

# load sqanti results
sqanti <- fread("../../04_evaluation/02_sqanti/data/downsampling_espresso_gencode/downsampling_espresso_classification.txt")

#### LOAD DATA
metadata <- fread("../../../00_metadata/data/pantranscriptome_samples_metadata.tsv") # metadata
annot <- fread("../../../../novelannotations/merge_espresso_downsample/merged_parsednames.gtf") # merged annotation
metadata <- metadata[merged_run_mode==TRUE]
metadata <-metadata[mixed_samples==FALSE]
#### PARSE ANNOTATION
colnames(annot) <- c("contig", "platform", "feature", "start", "end", "unk1", "strand", "unk2", "info")

annot <- annot[feature=="transcript",]
gc()
annot[, transcript := tstrsplit(info, "\"")[[2]]]
annot[, samples := tstrsplit(info, "\"")[[4]]]

# Split the tool_sample_pairs into individual pairs
annot[, samples_split := strsplit(samples, ",")]
# Transform the annotation to a long format to be able to compute sharing across samples, tools and populations
samples_split <- unique(unlist(annot$samples_split))
gc()
# Create a data.table with the transcripts and their corresponding tool_sample_pairs
annot_expanded <- annot[, .(transcript = transcript, sample = unlist(samples_split)), by = 1:nrow(annot)]
gc()


# add sqanti to annot_expanded
newannot <- sqanti[, .(isoform, structural_category)][annot_expanded, on=c("isoform"="transcript")]
newannot[, cell_line_id := tstrsplit(sample, "_")[[3]]]
newannotmeta <- metadata[, .(cell_line_id, population, map_reads_assemblymap, sample)][newannot, on="cell_line_id"]
newannotmeta[, trx_per_sample:=uniqueN(isoform), by="sample"]
newannotmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]
library(ggpubr)


categories <- c("FSM", "Intergenic", "NNC", "NIC",  "Fusion", "ISM", "Genic", "Antisense")
names(categories) <- unique(newannotmeta$structural_category)
newannotmeta$structural_category <- categories[newannotmeta$structural_category]
newannotmeta[, structural_category :=factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Intergenic", "Genic", "Fusion", "Antisense"))]
p <- ggplot(unique(newannotmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC")][, trx_per_sample_cat := uniqueN(isoform), by=c("sample", "structural_category")][, .(trx_per_sample_cat, structural_category,population, eur)]), 
       aes(x = eur, y = trx_per_sample_cat, fill = eur)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.9)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.6, width = 0.2, dodge.width = 0.8, aes(col=population)) +
  mytheme +
  scale_fill_manual(values = c("#466995", "#A53860")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = 0), 
               position = position_dodge(0.9), size=7*0.35) +
  labs(x = "", y = "# Discovered Transcripts", col="Population")+
  facet_wrap(~structural_category, nrow=1)+
  guides(fill="none")+
  scale_color_manual(values=popcols)+
  theme(legend.key.size = unit(0.2, "cm"))+
  geom_pwc(ref.group="European",
           method="t_test", label.size=7*0.35
  )+ylim(c(0,31000))
ggadjust_pvalue(
  p=p,
  p.adjust.method = "BH",
  label = "p.adj.format",
  hide.ns = NULL,
  output = c("plot"))


ggsave("../../../10_figures/01_plots/main/fig_04/violin_EURvsNonEUR_DiscoveredTranscripts_PerSqantiCategoryFSMtoNNC.pdf", dpi=700, width = 5.91, height = 2.7,  units = "in")



## NOW CHECK THE POPULATION SPECIFIC TRANSCRIPTS
newannotmeta <-newannotmeta[!sample%in%c("ITU1", "CEU1", "LWK1", "YRI1", "YRI2", "AJI1", "AJI2", "PEL1", "PEL2", "HAC2", "HAC2")]
newannotmeta[, sample_detected:=1][, sample_per_pop_detected:=sum(sample_detected), by=c("isoform", "population")]
newannotmeta[, i.sample:=NULL][, cell_line_id:=NULL][, nrow:=NULL]
newannotmeta[, populations_detected := uniqueN(population), by = "isoform"]

colsqanti <- c("#61814B",  "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353","#8EDE95")
names(colsqanti) <- unique(newannotmeta$structural_category)[c(1,4,3,2,7,5,8,6)]

newannotmeta[, novelty:=factor(fifelse(structural_category%in%c("FSM", "ISM"), "Annotated", "Novel"), levels=rev(c("Annotated", "Novel")))]
newannotmeta[, population:=factor(population, levels=c("MPC", "ITU", "CEU", "LWK", "YRI", "AJI", "PEL", "HAC"))]
ggplot(newannotmeta[populations_detected==1 & sample_per_pop_detected>=2], aes(x=population, fill=structural_category))+
  geom_bar(position="fill")+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))
ggplot(newannotmeta[populations_detected==1 & sample_per_pop_detected>=2 & structural_category%in%c("FSM", "ISM", "NIC", "NNC")], aes(x=population, fill=structural_category))+
  geom_bar(position="fill")+
  mytheme+
  scale_fill_manual(values=colsqanti)+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_fill(vjust=0.5))+
  labs(x="", y="Proportion of Population Specific Transcripts", fill="")
ggsave("figures/.pdf", dpi=700, width = 25, height = 14,  units = "cm")


popsp <-newannotmeta[populations_detected==1 & sample_per_pop_detected>=2 & structural_category%in%c("FSM", "ISM", "NIC", "NNC")]
isoform_counts <- popsp[novelty == "Annotated", uniqueN(isoform), by = population]

# Step 2: Reorder `population` based on the counts calculated in Step 1
# (Higher counts appear first)
popsp$population <- factor(popsp$population, 
                                  levels = isoform_counts[order(-V1)]$population)
ggplot(popsp, 
       aes(x=reorder(population,count ), fill=novelty))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=rev(c("#5E8531","#A2331D")))+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5))+
  labs(x="", y="Proportion of Population Specific Transcripts", fill="Intron Chain")+
  theme(legend.position="top")
ggsave("../../../10_figures/downsampling/barplotStack.annotatedVSnovelIntronChain_POPspecificTrx.pdf", dpi=700, width = 17, height = 14,  units = "cm")



