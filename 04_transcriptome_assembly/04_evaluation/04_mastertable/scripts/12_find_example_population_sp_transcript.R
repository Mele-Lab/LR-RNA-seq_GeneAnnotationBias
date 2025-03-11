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
metadata <- fread("00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- metadata$color_pop
names(popcol) <- metadata$population

# Find an example
counts <- fread("../novelannotations/quantifications/kallisto_quant/matrix.abundance.tpm.tsv")
counts <- melt(counts, value.name = "TPM", id.vars = "transcript_id", variable.name = "cell_line_id_unparsed")
counts[, cell_line_id:=gsub("_1", "", cell_line_id_unparsed)]
counts <- metadata[,.(cell_line_id, sample, population)][counts, on="cell_line_id"]



ggplot(counts[transcript_id == "transcript_148736"][order(sample), thresh := TPM < 0.1][TPM>0],
       aes(x = "", y = TPM, color = population, alpha = thresh)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#393939") +
  ggbeeswarm::geom_quasirandom(size=3) +  # Adjust width to control separation between bars
  mytheme +
  scale_color_manual(values = popcol, name="") +  # Set all bars to the same fill color
  labs(fill = "Population", x = "transcript_148736", y="TPM") +
  scale_alpha_manual(values = c(1, 0.3), name="", labels=c("TRUE"="TPM >= 0.1", "FALSE"="TPM < 0.1")) +
  guides(fill = "none")+
  scale_y_continuous(trans="log10")+
  annotation_logticks(sides = "l")
ggsave(paste0("../../../10_figures/01_plots/main/fig_03/trxplot.popSpecific_transcript_example_expression.pdf"), dpi=700, width = 2, height = 2.25,  units = "in")



library(ggtranscript)
MYTRX <- "transcript_148736"
MYGEN<- "ENSG00000113387.13"

# original gtf in CRGcluster /users/rg/fdegalez/data/annotation/241018_v47_poder_merge.placeholder_gene_name.withCDS.trBiotype.gtf
gtf <- fread("../novelannotations/ENSG00000113387.13.gtf")
# gtf <- fread("../novelannotations/241018_v47_poder_merge.placeholder_gene_name.withCDS.trBiotype.gtf")
colnames(gtf) <- c("chr", "source", "feature", "start", "end", "uk1", "strand" , "uk2", "info")
gtf <- gtf[feature!="gene"]

gtf <- gtf[grepl("ENST", info)| grepl("transcript_148736", info)]
gtf[, transcript := sub('.*transcript_id "([^"]+)".*', '\\1', info)]

gtf[, seqnames:=chr][, chr:=NULL]
gtf[, type:=feature][, feature:=NULL]
gtf[, source:=fifelse(source=="HAVANA", "GENCODE", "PODER")]



# prepare data for plotting
subgtf <- gtf
subgtf[, newtype:= type]
subgtf[type=="CDS", type:="exon"]
subgtf[newtype=="CDS", transcript:=paste0("mock", transcript)]
subgtf_rescaled  <- shorten_gaps(
  exons = subgtf[type=="exon"], 
  introns = to_intron(subgtf[type=="exon"], "transcript"), 
  group_var = "transcript"
)

# subgtf[transcript%in%c("ENST00000515355.5", "transcript_148736")]

setDT(subgtf_rescaled)
subgtf_rescaled[, transcript:=gsub("mock", "", transcript)]
ggplot(subgtf_rescaled[type=="exon" & newtype!="CDS"],
            aes(xstart = start,
                xend = end,
                y = transcript)) +
  geom_rect(aes(xmin=528, xmax=7576, ymin=13.5, ymax=14.5), fill="#ed9693", alpha=0.1)+
  geom_rect(data = subgtf_rescaled[type=="exon" & transcript==MYTRX],
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5)+
  geom_range(aes(fill = source, color=source), linewidth = 0, height=0.25) +
  geom_intron(data = unique(subgtf_rescaled[type=="intron" & newtype!="CDS"][, .(strand, start, end, transcript)]),
              aes(strand = strand), linewidth=0.25,
              arrow.min.intron.length = 300,
              arrow = grid::arrow(ends = "last", length = grid::unit(0.05, "inches")))+
  mytheme+
  labs(y="", title="SUB1", x="Distance from most upstream TSS")+
  scale_color_manual(values= c("#7F675B","#A2AD59"), name="Source")+
  scale_fill_manual(values= c("#7F675B","#A2AD59"), name="Source")+
  theme(legend.position=c(0.88, 0.8),
        plot.margin = margin(1, 1, 1, 1))+
  geom_range(
    data = subgtf_rescaled[type=="exon" & newtype == "CDS"],
    aes(fill=source),
    color = NA, 
    linewidth = 1, height=0.85
  )+
  theme(plot.title = element_text(face = "bold.italic"))+
  annotate(geom="text", label="PEL specific", y=14, vjust=0.5, hjust=0, x=7700, size=6*0.35, color="darkred", fontface="bold")
ggsave("10_figures/01_plots/main/fig_03/trxplot.popSpecific_transcript_example.pdf", dpi=700, width = 3, height = 2,  units = "in")


ggplot(subgtf_rescaled[type=="exon" & newtype!="CDS"],
       aes(xstart = start,
           xend = end,
           y = transcript)) +
  geom_rect(aes(xmin=528, xmax=7576, ymin=13.5, ymax=14.5), fill="#ed9693", alpha=0.1)+
  geom_rect(data = subgtf_rescaled[type=="exon" & transcript==MYTRX],
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5)+
  geom_range(aes(fill = source, color=source), linewidth = 0, height=0.25) +
  geom_intron(data = unique(subgtf_rescaled[type=="intron" & newtype!="CDS"][, .(strand, start, end, transcript)]),
              aes(strand = strand), linewidth=0.25,
              arrow.min.intron.length = 300,
              arrow = grid::arrow(ends = "last", length = grid::unit(0.05, "inches")))+
  mytheme+
  labs(y="", title="Novel\nSplice\nSite", x="")+
  scale_color_manual(values= c("#7F675B","#A2AD59"), name="Source")+
  scale_fill_manual(values= c("#7F675B","#A2AD59"), name="Source")+
  theme(legend.position=c(0.88, 0.8),
        plot.margin = margin(1, 1, 1, 1))+
  geom_range(
    data = subgtf_rescaled[type=="exon" & newtype == "CDS"],
    aes(fill=source),
    color = NA, 
    linewidth = 1, height=0.85
  )+
  theme(plot.title = element_text(face = "bold"))+
  annotate(geom="text", label="PEL specific", y=14, vjust=0.5, hjust=0, x=7700, size=6*0.35, color="darkred", fontface="bold")+
  coord_cartesian(xlim = c(2220, 2240), ylim = c(11, 14.25))+
  theme(legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  geom_vline(xintercept=2230, linetype="dashed")
ggsave("10_figures/01_plots/main/fig_03/trxplot.popSpecific_transcript_example_zoom.pdf", dpi=700, width = 0.5, height = 2,  units = "in")

             