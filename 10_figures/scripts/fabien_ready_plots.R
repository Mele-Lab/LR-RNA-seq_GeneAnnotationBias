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
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}

#### MAPPED READS

## personal assembly part
data <- fread("10_figures/data/primaryAlignment_forBarplot.txt")
data <- data[type!="total"]
data[, newtype:=ifelse(type%in%c("bothPG", "PG", "PG_T2T", "T2T_bothPG"), "Only\nPersonal Assembly",
                       ifelse(type%in%c("oPG_bothRef", "oPG_hg38", "bothRef", "hg38" ), "Only\nGRCh38",
                              ifelse(type%in%c("oPG", "T2T","unmapped", "oPG_T2T"), "None",
                                     ifelse(type%in%c("all", "PG_bothRef"),"Both\nAssemblies",
                                            ifelse(type%in%c("hg38_bothPG", "PG_hg38"), "Both\nAssemblies","MISSING")))))]
data[, total:=sum(freq), by="genome"]
data[, newfreq:=sum(freq), by=c("genome", "newtype")]
data[, newperc:=paste0(round(newfreq/total*100, digits=2), "%")]


mydata <- unique(data[, meancount :=mean(newfreq), by="newtype"][, sdcount:=sd(newfreq), by="newtype"][, .(meancount, sdcount, newtype)])
mydata[, percentcount:=paste0(round(meancount/sum(meancount)*100, digits=3), " %")]


ggplot(mydata, 
       aes(x = newtype, fill = newtype)) +
  geom_col( alpha = 0.8, aes(y=meancount)) +
  geom_errorbar(width=0.3, aes(ymin=meancount-sdcount, ymax=meancount+sdcount)) +
  
  ggbeeswarm::geom_quasirandom(data=unique(data[, .(genome, newfreq, newtype, newperc)]), size=0.5,alpha = 0.4, aes( col=genome,x=newtype,  y = newfreq)) +
  labs(y = "# Reads", x = "Reads Mapping to") +
  scale_fill_manual(values = c("#4A9F37","darkgrey", "#368282", "darkred", "#4A9F37",  "#4A9F37")) +
  guides(alpha = "none", fill="none", color="none") +
  mytheme+
  scale_y_continuous(trans="log10")+
  annotation_logticks(sides = "l", outside=F, linewidth=0.25)+
  scale_color_manual(values=rep("black", 6))+
  geom_text(aes(label = percentcount, y=meancount+sdcount*10),  vjust=0, size=6*0.35)+
  geom_text(aes(label = round(meancount, digits=0), y=meancount-sdcount),  vjust=3, size=6*0.35)+
  coord_cartesian(ylim=c(1, 100e6), xlim=c(0.85,4))+
  scale_x_discrete(limits=c("None", "Both\nAssemblies", "Only\nGRCh38","Only\nPersonal Assembly"),
                   labels=c("None"="Unmapping", 
                            "All Assemblies" ="All\nAssemblies",
                            "Only\nGRCh38/T2T"= "Only\nGRCh38\nand/or T2T", 
                            "PG and\nGRCh38"="GRCh38+\nPersonal\nAssembly", 
                            "PG and\nT2T"="T2T+\nPersonal\nAssembly", "Only\nPersonal Assembly"="Only\nPersonal\nAssembly"))
ggsave("10_figures/01_plots/main/fig_05/barplot_MappedReads.pdf", dpi=700, width = 2.5, height = 2.25,  units = "in")

##### FULL UNCOLLAPSED PLOT FOR SUPP
data[, newtype:=ifelse(type%in%c("bothPG", "PG"), "Only\nPersonal Assembly",
                       ifelse(type%in%c("oPG_bothRef", "oPG_hg38", "oPG_T2T", "bothRef"), "Only\nGRCh38+T2T",
                              ifelse(type%in%c("oPG"), "None",
                                     ifelse(type=="unmapped", "None", 
                                            ifelse(type%in%c("all", "PG_bothRef"),"All Assemblies",
                                                   ifelse(type%in%c("hg38_bothPG", "PG_hg38"), "PG and\nGRCh38",
                                                          ifelse(type%in%c("T2T_bothPG", "PG_T2T"),"PG and\nT2T",
                                                                 ifelse(type%in%c("hg38"), "Only GRCh38",
                                                                        ifelse(type%in%c("T2T"), "Only T2T", "MISSING")))))))))]

data[, total:=sum(freq), by="genome"]
data[, newfreq:=sum(freq), by=c("genome", "newtype")]
data[, newperc:=paste0(round(newfreq/total*100, digits=2), "%")]


mydata <- unique(data[, meancount :=mean(newfreq), by="newtype"][, sdcount:=sd(newfreq), by="newtype"][, .(meancount, sdcount, newtype)])
mydata[, percentcount:=paste0(round(meancount/sum(meancount)*100, digits=3), " %")]
mydata[, newcolor:=ifelse(newtype%in%c("Only T2T", "Only GRCh38"), "1 Reference Assembly",
                          ifelse(newtype=="None", "Unmapping",
                                 ifelse(newtype=="Only\nPersonal Assembly", "1 Personal Assembly", "At least 2 Assemblies")))]

ggplot(mydata, 
       aes(x = reorder(newtype, -meancount), fill = newtype)) +
  geom_col( alpha = 0.8, aes(y=meancount, fill=newcolor)) +
  geom_errorbar(width=0.3, aes(ymin=meancount-sdcount, ymax=meancount+sdcount)) +
  
  ggbeeswarm::geom_quasirandom(data=unique(data[, .(genome, newfreq, newtype, newperc)]), size=0.5,alpha = 0.4, aes( col=genome,x=newtype,  y = newfreq)) +
  labs(y = "# Reads", x = "Reads Mapping to", fill="") +
  guides(alpha = "none",  color="none") +
  mytheme+
  scale_y_continuous(trans="log10")+
  annotation_logticks(sides = "l", outside=F)+
  scale_color_manual(values=rep("black", 6))+
  geom_text(aes(label = percentcount, y=meancount+sdcount*10),  vjust=0, size=6*0.35)+
  geom_text(aes(label = round(meancount, digits=0), y=meancount-sdcount),  vjust=3, size=6*0.35)+
  coord_cartesian(ylim=c(1, 100e6),xlim=c(0.75,8))+
  scale_x_discrete(labels=c("None"="None\n(Unmapping)", 
                            "All Assemblies" ="All\nAssemblies",
                            "PG and\nGRCh38"="GRCh38 +\nPersonal\nAssembly", 
                            "PG and\nT2T"="T2T+\nPersonal\nAssembly", 
                            "Only\nPersonal Assembly"="Only\nPersonal\nAssembly",
                            "Only\nGRCh38+T2T"="GRCh38\n+ T2T"))+
  scale_fill_manual(values = c("At least 2 Assemblies"="#4A9F37","Unmapping"="darkgrey", "1 Personal Assembly"="darkred","1 Reference Assembly"="#368282"))+
  theme(legend.position=c(0.85, 0.9))

ggsave("10_figures/01_plots/supp/38_map_pg/barplot_MappedReads.pdf", dpi=700, width = 4.5, height = 2.75,  units = "in")


# plot gene density
data <- fread("10_figures/data/fabien_data/table1_geneDensity.tsv")
colnames(data)[grepl("pop", colnames(data))] <- "fastq"
data[, newUniqueRegionsOf:=fifelse(!is.na(uniqueRegionsOf) & !uniqueRegionsOf%in%c("GRCh38", "T2T"), "Personal\nAssembly", uniqueRegionsOf )]
data[is.na(newUniqueRegionsOf), newUniqueRegionsOf:="Shared"]
data[, newfasta:=ifelse(!FASTA%in%c("T2T", "GRCh38"), "Personal\nAssembly", FASTA)]
data[, category:=paste0("found in ", newfasta, " ", typeRegion, " with ", newUniqueRegionsOf  )]

median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(median(x), 2)))
}

data[, typeRegion:=factor(typeRegion, levels=c("all","shared", "notshared"))]
ggplot(data[FASTA != "T2T" & typeRegion != "all" & uniqueRegionsOf != "T2T"], 
       aes(x = newfasta, y = density, fill = typeRegion)) +
  geom_violin(alpha = 0.5, position = position_dodge(0.75)) + 
  ggbeeswarm::geom_quasirandom(aes(color = fastq, group = interaction(newfasta, typeRegion)), 
                               size = 0.5, alpha = 0.5, dodge.width = 0.75, show.legend = FALSE) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.75), 
               color = "black", show.legend = FALSE) + 
  mytheme +
  labs(y = "Density of Discovered Genes\n(Genes/Mb)", x = "Genome Assembly", fill = "Genes sitting in") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA))) +
  scale_color_manual(values = rep("black", 8)) +
  scale_fill_manual(values = rev(c("darkred", "#4A9F37")), labels = c("notshared" = "Unique Region", "shared" = "Shared Region")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1), size = 5 * 0.35, position = position_dodge(0.75)) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 0.2), size = 6 * 0.35, position = position_dodge(0.75)) +
  
  geom_text(data = unique(data[FASTA != "T2T" & typeRegion != "all" & uniqueRegionsOf != "T2T"][, mediann := round(median(nbGenes), 0), by = .(newfasta, typeRegion)][, .(newfasta, typeRegion, mediann)]),
            aes(label = mediann, y = -0.5), size = 6 * 0.35, position = position_dodge(0.75)) +
  
  annotate(geom = "text", label = "density =", size = 6 * 0.35, y = 0.2, x = 0.6, hjust = 1) +
  annotate(geom = "text", label = "number =", size = 6 * 0.35, y = -0.5, x = 0.6, hjust = 1) +
  coord_cartesian(xlim = c(0.75, 2))+
  theme(legend.position = "top")+
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 3)  # Increase legend square size
    )
  )
ggsave("10_figures/01_plots/main/fig_05/violin_genes_discovered.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")
# plot transcript density-------------------------------------------------------------
data <- fread("10_figures/data/fabien_data/table1_transcriptDensity.txt")
colnames(data)[grepl("pop", colnames(data))] <- "fastq"
data[, newUniqueRegionsOf:=fifelse(!is.na(uniqueRegionsOf) & !uniqueRegionsOf%in%c("GRCh38", "T2T"), "Personal\nAssembly", uniqueRegionsOf )]
data[is.na(newUniqueRegionsOf), newUniqueRegionsOf:="Shared"]
data[, newfasta:=ifelse(!FASTA%in%c("T2T", "GRCh38"), "Personal\nAssembly", FASTA)]
data[, category:=paste0("found in ", newfasta, " ", typeRegion, " with ", newUniqueRegionsOf  )]

median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(median(x), 2)))
}

data[, typeRegion:=factor(typeRegion, levels=c("all","shared", "notshared"))]
ggplot(data[FASTA != "T2T" & typeRegion != "all" & uniqueRegionsOf != "T2T"], 
       aes(x = newfasta, y = density, fill = typeRegion)) +
  geom_violin(alpha = 0.5, position = position_dodge(0.75)) + 
  ggbeeswarm::geom_quasirandom(aes(color = fastq, group = interaction(newfasta, typeRegion)), 
                               size = 0.5, alpha = 0.5, dodge.width = 0.75, show.legend = FALSE) + 
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.75), 
               color = "black", show.legend = FALSE) + 
  mytheme +
  labs(y = "Density of Discovered Transcripts\n(Transcripts/Mb)", x = "Genome Assembly", fill = "Transcripts in") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, color = NA))) +
  scale_color_manual(values = rep("black", 8)) +
  scale_fill_manual(values = rev(c("darkred", "#4A9F37")), labels = c("notshared" = "Unique Region", "shared" = "Shared Region")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1), size = 5 * 0.35, position = position_dodge(0.75)) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 0.2), size = 6 * 0.35, position = position_dodge(0.75)) +
  
  geom_text(data = unique(data[FASTA != "T2T" & typeRegion != "all" & uniqueRegionsOf != "T2T"][, mediann := round(median(nbTr), 0), by = .(newfasta, typeRegion)][, .(newfasta, typeRegion, mediann)]),
            aes(label = mediann, y = -0.5), size = 6 * 0.35, position = position_dodge(0.75)) +
  
  annotate(geom = "text", label = "density =", size = 6 * 0.35, y = 0.2, x = 0.6, hjust = 1) +
  annotate(geom = "text", label = "number =", size = 6 * 0.35, y = -0.5, x = 0.6, hjust = 1) +
  coord_cartesian(xlim = c(0.55, 2))+
  theme(legend.position = "top")+
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 3)  # Increase legend square size
    )
  )
ggsave("10_figures/01_plots/main/fig_05/violin_trx_discovered.pdf", dpi=700, width = 2.25, height = 2.25,  units = "in")

# all
data <- data[newfasta != "GRCh38" & typeRegion != "all" & uniqueRegionsOf != "GRCh38"]
data[, newfasta:=factor(newfasta, levels=c("T2T", "Personal\nAssembly"))]
ggplot(data, 
       aes(x = newfasta, y = density, fill = typeRegion)) +
  geom_violin(alpha = 0.5,  position = position_dodge(0.75)) + # Remove violin border lines
  ggbeeswarm::geom_quasirandom(aes(color = fastq, group = interaction(newfasta, typeRegion)), 
                               size = 0.5, alpha = 0.5, dodge.width = 0.75, show.legend = FALSE) + # Exclude from legend
  geom_boxplot(outliers = FALSE, width = 0.15, position = position_dodge(0.75), 
               color = "black", show.legend = FALSE) + # Remove boxplot fill from legend
  mytheme +
  labs(y = "Density of Discovered Transcripts\n(Transcripts/Mb)", x = "Genome Assembly", fill = "Transcripts\nsitting in") +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, color = NA)) # Override legend to display plain squares
  ) +
  scale_color_manual(values = rep("black", 8)) +
  scale_fill_manual(values = rev(c("darkred", "#4A9F37")), labels = c("all"="All Genomes Shared Region","notshared" = "Unique Region", "shared" = "Both Assemblies")) +
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1), size = 5 * 0.35, position = position_dodge(0.75)) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 0.2), size = 6 * 0.35, position = position_dodge(0.75)) +
  geom_text(data = unique(data[, mediann := round(median(nbTr),0), by = .(newfasta, typeRegion)][, .(newfasta, typeRegion, mediann)]),
            aes(label = mediann, y = -0.5), size = 6 * 0.35, position = position_dodge(0.75))+
  annotate(geom="text", label="density =", size=6*0.35, y=0.2, x=0.6, hjust=1)+
  annotate(geom="text", label="number =", size=6*0.35, y=-0.5, x=0.6, hjust=1)+
  coord_cartesian(xlim = c(0.75,2)) +# Add extra space by setting x-axis limits
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 3)  # Increase legend square size
    )
  )+
  theme(legend.position = "top")
ggsave("10_figures/01_plots/supp/38_map_pg/violin_trx_discovered_all.pdf", dpi=700, width = 2.5, height = 2.25,  units = "in")

# repetitiveness-----------------------------------------------------------------------------

# pop contains: region from left not found in right

random <- fread("10_figures/data/fabien_data/shannonEntropy_random.tsv")
colnames(random) <- c("pop", "contigs", "sumKmers", "shannonEntropy")
random <- random[pop!="hg38_hg38"]
random[, type:="Background\nMatching\nRegions"]
data1 <- fread("10_figures/data/fabien_data/shannonEntropy_specRegion1.tsv")
data2 <- fread("10_figures/data/fabien_data/shannonEntropy_specRegion2.tsv")
data <- rbind.data.frame(data1, data2)
data <- data[pop!="hg38_hg38"]
data[, type:="Assembly\nSpecific\nRegions"]

data <- rbind.data.frame(data, random)


median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(median(x), digits=2)))
}

data[, ref:=tstrsplit(pop, "_")[[1]]]
data[, comparison:=tstrsplit(pop, "_")[[2]]]
data[, type_comparison:=fifelse(comparison%in%c("hg38", "T2T") | ref%in%c("hg38", "T2T"), "containsref", "pgaonly")]
data <- data[type_comparison=="containsref"]
data[, reftype:=fifelse(ref=="hg38", "GRCh38",
                        fifelse(ref=="T2T", "T2T", "Personal Genome\nAssemblies"))]
data[, comptype:=fifelse(comparison=="hg38", "Compared to GRCh38",
                         fifelse(comparison=="T2T", "Compared to T2T", "Personal Genome\nAssemblies"))]

ggplot(data[reftype== "Personal Genome\nAssemblies"], aes(x=type, y=as.numeric(shannonEntropy), fill=type))+
  geom_boxplot(outliers = T, alpha=0.5)+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -1), size=5*0.35) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = 12.5), size=6*0.35)+
  ggpubr::stat_compare_means(ref.group = "FSM", method.args = list(alternative="less"), label.y=110)+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  guides(fill="none")+
  facet_wrap(~comptype)+
  labs(x="", y="Shannon Entropy\nof Region k-mers")+
  ggpubr::stat_compare_means(comparisons=list(c("Background\nMatching\nRegions", "Assembly\nSpecific\nRegions")),
                             method = "wilcox.test", 
                             method.args = list(alternative = "two.sided"),size=7*0.35, label.y=14)+
  scale_x_discrete(limits=c("Background\nMatching\nRegions", "Assembly\nSpecific\nRegions"))+
  ylim(c(-1.5, 15.5))
ggsave("10_figures/01_plots/supp/37_nonref_seq/violin_shannon_entropy.pdf", dpi=700, width = 5, height = 2.5,  units = "in")




#### regions length ----------------------------
data <-fread("10_figures/data/fabien_data/sizeSpecRegions.txt")
data[, type_comparison:=fifelse(genome%in%c("hg38", "T2T") | ref%in%c("hg38", "T2T"), "containsref", "pgaonly")]
data <- data[genome!=ref]


ggplot(data, aes(x=size/1000, col=testedAssemblies))+
  geom_density(alpha=0.5)+
  mytheme+
  scale_x_continuous(trans="log10")+
  guides(col="none")+
  labs(x="Length (Kb) of Sequences not Found in Other Genomes", y="Density")+
  annotation_logticks(sides = "b")+
  theme(plot.margin = margin(t = 10, r = 15, b = 10, l = 10))

# Compute density
density_data <- data[, {
  d <- density(size.log10, n = 512)
  data.table(x = d$x, y = d$y)
}, by=testedAssemblies]

# Define a common x grid
common_x <- seq(min(density_data$x), max(density_data$x), length.out = 1000)

# Interpolate y values for the common grid
interpolated_density <- density_data[, {
  valid_x <- common_x[common_x >= min(x) & common_x <= max(x)] # Restrict common_x
  y_interp <- approx(x, y, xout = valid_x)$y
  data.table(x = valid_x, y = y_interp)
}, by = testedAssemblies]

# Summarize mean and standard deviation on the common grid
summary_density <- interpolated_density[, .(
  mean_y = mean(y, na.rm = TRUE),
  sd_y = sd(y, na.rm = TRUE)
), by = x]

# Add ymin and ymax for the shaded area
summary_density[, `:=`(
  ymin = mean_y - sd_y,
  ymax = mean_y + sd_y
)]

# Plot
ggplot(summary_density, aes(x = x)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "black", alpha = 0.2) +
  geom_line(aes(y = mean_y), color = "black", size = 1*0.35) +
  mytheme +
  labs(
    x = "log10(Length of Sequences not Found in Other Assemblies)",
    y = "Density"
  )+
  annotation_logticks(sides="b")+
  xlim(c(0,9))+
  geom_vline(xintercept=log10(median(data$size)), linetype="dashed", color="darkgrey")+
  annotate(geom="text", label=paste0("Median size = ",median(data$size)), x=log10(median(data$size)), hjust=-0.2, y=0.75, size=7*0.35, )
ggsave("10_figures/01_plots/supp/37_nonref_seq/density_FragmentsLenght.pdf", dpi=700, width = 2.75, height = 2.75,  units = "in")

##### non-reference regions

data <- fread("10_figures/data/fabien_data/genomicComparisonWithN.txt")
ns<-data[testedAssemblies=="hg38_hg38", totalSize]
data[, totalSize:=ifelse(ref=="hg38", totalSize-ns, totalSize)]
data$genome_size <- as.numeric(data$genome_size)
data[, genome_size:=ifelse(ref=="hg38", genome_size-ns, genome_size)]
data[, totalSize.prct:=round(totalSize/genome_size*100, digits=2)]
data$ref <- factor(data$ref, levels = c("hg38","T2T", setdiff(unique(data$ref), c("hg38","T2T"))))
data$genome <- factor(data$genome, levels = c("hg38","T2T", setdiff(unique(data$ref), c("hg38","T2T"))))

ggplot(data, aes(y=ref, x=genome))+
  geom_tile(color = "white", aes(fill=totalSize/1000000)) +
  scale_fill_gradient2(low = "white", high = "#A0000D", na.value = "#9AC7CB") +
  theme_minimal() +
  labs(x = "Against Target Genome",
       y="Map Query Genome",
       fill = "Unmapping\nRegions (Mb)") +
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title=element_text(face="bold", hjust=0.5)) +
  geom_text( aes(label=paste0(totalSize.prct, " %"), size=abs(totalSize.prct)),  na.rm =T)+
  # coord_fixed()+
  scale_size_continuous(range=c(4,6)*0.35)+
  guides(size="none")+
  scale_x_discrete(labels=c("hg38"="GRCh38"))+
  scale_y_discrete(labels=c("hg38"="GRCh38"))+
  ggside::geom_ysidecol(data=unique(data[, .(ref, genome_size)]), aes(y=ref, x=genome_size/1000000), fill="darkgrey", col="darkgrey")+
  ggside::geom_ysidetext(data=unique(data[, .(ref, genome_size)]), aes(y=ref, label=paste0(round(genome_size/1000000), " Mb"), x=genome_size/1000000), 
                         col="black", size=6*0.35, hjust=-0.05)+
  ggside::scale_ysidex_continuous(guide = guide_axis(angle = 45), limits = c(0, 7250), guide_axis("TEST"))+
  theme(ggside.panel.scale = .3)  
ggsave("10_figures/01_plots/supp/37_nonref_seq/heatmap_missingregions.pdf", dpi=700, width = 4.5, height = 3,  units = "in")



ggplot(data[genome %in% c("hg38", "T2T") & ref != "hg38" & ref != "T2T"], aes(x = ref, y = totalSize / 1e6)) +
  geom_col(fill = "#9D1616") +
  mytheme +
  labs(x = "", y = "Total Size of Not Found Sequences (Mb)") +
  facet_wrap(~genome,
             labeller = labeller(genome = c(
               "hg38" = "Regions not found in GRCh38",
               "T2T" = "Regions not found in T2T"
             ))) +
  geom_hline(data = data[genome %in% c("hg38", "T2T") & ref != "hg38" & ref != "T2T", 
                         .(median_totalSize = median(totalSize / 1e6)), by = genome],
             aes(yintercept = median_totalSize),
             linetype = "dashed", col = "darkgrey") +
  geom_text(data = data[genome %in% c("hg38", "T2T") & ref != "hg38" & ref != "T2T", 
                        .(median_totalSize = median(totalSize / 1e6)), by = genome],
            aes(x = 3, y = median_totalSize, 
                label = paste0("median = ", round(median_totalSize, digits = 1), " Mb")),
            vjust = -0.5, size = 6 * 0.35)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("10_figures/01_plots/supp/37_nonref_seq/barplot_missingregions.pdf", dpi=700, width = 5, height = 2.5,  units = "in")

### PLOT NUMBER OF GENES PARTIALLY INTERSECTING NON REF REGIONS
data <- fread("10_figures/data/fabien_data/table1_geneTypeRegion.tsv")
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  paste0("median = ", round(median(x), digits=0))))
}
ggplot(data[uniqueRegionsOf=="hg38"], aes(x=type, y=freq, fill=type))+
  geom_violin(alpha=0.75)+
  geom_boxplot(outliers=F, width=0.05)+
  ggbeeswarm::geom_quasirandom(size=1)+
  scale_x_discrete(labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  labs(x="Relative Location of Gene\nto non-Reference Region", y="# Discovered Gene\nin non-Reference Regions")+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -4), size=5*0.35) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y =5), size=6*0.35)+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"))+
  guides(fill="none")
ggsave("10_figures/01_plots/main/fig_05/violin_types_genes_nonref.pdf", dpi=700, width = 2.75, height = 2.75,  units = "in")


##### PLOT NUMBER OF TRX PARTIALLY INTERSECTING NON REF REGIONS

# Meaning of the columns:
# *isIncluded=T if the transcript is fully included in the region (=included in category)
# *nb_ex: number of exons in the tr
# *nb_ex_intersected: nb of exons intersecting a region
# *nb_ex_included: nb of exons fully includin in a region	
# *propExIncl: proportion of exons included in the region
# *category: intronic = the region is fully included in an intronic part of the tr / included = the transcript is fully included in the region / intersect
# *pop: sample testes
# *unique_region: origin of the region (T2T or hg38)

data <- fread("10_figures/data/fabien_data/trCategoriesIntersectNonRefRegion.txt")
data[, counts:=uniqueN(tr_id), by=c("pop", "category", "unique_region")]
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  paste0(round(median(x), digits=2))))
}
ggplot(unique(data[unique_region=="hg38", .(category, pop, unique_region, counts)]), aes(x=category, y=counts, fill=category))+
  geom_violin(alpha=0.75)+
  geom_boxplot(outliers=F, width=0.05)+
  ggbeeswarm::geom_quasirandom(size=1)+
  scale_x_discrete(labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  labs(x="Relative Location of Transcript\nto non-Reference Region", y="# Discovered Transcripts\nin non-Reference Regions")+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -4), size=5*0.35) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y =5), size=6*0.35)+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"))+
  guides(fill="none")
ggsave("10_figures/01_plots/main/fig_05/violin_types_transcripts_nonref.pdf", dpi=700, width = 2, height = 2.25,  units = "in")
ggplot(unique(data[unique_region=="T2T", .(category, pop, unique_region, counts)]), aes(x=category, y=counts, fill=category))+
  geom_violin(alpha=0.75)+
  geom_boxplot(outliers=F, width=0.05)+
  ggbeeswarm::geom_quasirandom()+
  scale_x_discrete(labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  labs(x="Relative Location of Transcript\nto non-Reference Region", y="# Discovered Transcripts\nin non-Reference Regions")+
  mytheme+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -4), size=5*0.35) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y =5), size=6*0.35)+
  facet_wrap(
    ~unique_region,
    labeller = labeller(unique_region = c(
      "hg38" = "Regions not found in GRCh38",
      "T2T" = "Regions not found in T2T"
    )))+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"))+
  guides(fill="none")
ggsave("10_figures/01_plots/supp/38_map_pg/violin_types_transcripts_nonref.pdf", dpi=700, width = 2, height = 2.25,  units = "in")


data[is.na(nb_ex_intersected), nb_ex_intersected:=0]
data[is.na(nb_ex_included), nb_ex_included:=0]
data <- data[nb_ex_intersected!=0]
data[, nb_trx_with_exons_intersected:=uniqueN(tr_id), by=c("pop", "unique_region","category", "nb_ex_intersected")]
data[, `:=`(sd_nb_trx_with_exons_intersected=sd(nb_trx_with_exons_intersected)), by=c("unique_region", "nb_ex_intersected","category")]
data[, `:=`(mean_nb_trx_with_exons_intersected=mean(nb_trx_with_exons_intersected)), by=c("unique_region", "nb_ex_intersected","category")]

ggplot(unique(data[, .(mean_nb_trx_with_exons_intersected, unique_region, sd_nb_trx_with_exons_intersected,nb_ex_intersected, category)]), 
       aes(x=nb_ex_intersected, y=mean_nb_trx_with_exons_intersected, fill=category))+
  geom_col( position=position_dodge(preserve = "single"))+
  geom_errorbar(aes(ymin = mean_nb_trx_with_exons_intersected-sd_nb_trx_with_exons_intersected, ymax=mean_nb_trx_with_exons_intersected+sd_nb_trx_with_exons_intersected),
                position=position_dodge( preserve = "single"))+
  mytheme+
  labs(x="# Exons per Transcript Intersecting non-Reference Regions",
       y="Mean # of Transcripts\nIntersecting non-Reference Regions",
       fill="Relative location\nof Transcript to\n non-reference Region")+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
                    labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  facet_wrap(~unique_region,
             labeller = labeller(unique_region = c(
               "hg38" = "Regions not found in GRCh38",
               "T2T" = "Regions not found in T2T"
             )))
ggsave("10_figures/01_plots/supp/38_map_pg/barplot_num_exons_intersecting_nonref.pdf", dpi=700, width = 4.5, height = 2.75,  units = "in")

# ggplot(unique(data[unique_region=="hg38", .(mean_nb_trx_with_exons_intersected, unique_region, sd_nb_trx_with_exons_intersected,nb_ex_intersected, category)]), 
#        aes(x=nb_ex_intersected, y=mean_nb_trx_with_exons_intersected, fill=category))+
#   geom_col( position=position_dodge(0.75,preserve = "single"))+
#   geom_errorbar(aes(ymin = mean_nb_trx_with_exons_intersected-sd_nb_trx_with_exons_intersected, ymax=mean_nb_trx_with_exons_intersected+sd_nb_trx_with_exons_intersected),
#                 position=position_dodge(0.75, preserve = "single"),
#                 linewidth=0.1)+
#   mytheme+
#   labs(x="# Exons per Transcript Intersecting\nnon-Reference Regions",
#        y="Mean # of Transcripts\nIntersecting non-Reference Regions",
#        fill="Relative location\nof Transcript to\nnon-reference Region")+
#   scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
#                     labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
#   theme(legend.position = c(0.7, 0.7))
# ggsave("10_figures/01_plots/main/fig_05/barplot_num_exons_intersecting_nonref.pdf", dpi=700, width = 2.5, height = 2.75,  units = "in")

data[, newnb_ex_intersected:=factor(fifelse(nb_ex_intersected>=10, ">=10", as.character(nb_ex_intersected)), levels=c("1","2","3","4","5","6","7","8","9",">=10"))]

ggplot(data[unique_region=="hg38"], 
       aes(x=newnb_ex_intersected, y=after_stat(count)/6, fill=category), stat="count")+
  geom_bar()+
  mytheme+
  labs(x="# Exons per Transcript Intersecting\nnon-Reference Regions",
       y="Mean # of Transcripts\nIntersecting non-Reference Regions",
       fill="Relative location\nof Transcript to\nnon-reference Region")+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
                    labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  theme(legend.position = c(0.7, 0.7))
ggsave("10_figures/01_plots/main/fig_05/barplot_num_exons_intersecting_nonref.pdf", dpi=700, width = 2.5, height = 2.75,  units = "in")

### READS used to create a transcript

data <- fread("10_figures/data/fabien_data/table2_readQuantifFollowing.tsv")
colnames(data)[grepl("pop", colnames(data))] <- "fastq"

data[, "Reads Mapping to Reference":=nbReads.ref/nbReads.PG*100 ]
data[, "Reads Supporting Transcripts\nin Reference":=nbReads.refTr/nbReads.PG*100 ]
datalong <- melt(data, measure.vars = c("Reads Mapping to Reference",  "Reads Supporting Transcripts\nin Reference"))

datalong[, newvalue:=cut(value, 
                         breaks = c(0,25,50,75,100,101),  # Add 101 to include 100 in the last range                         include.lowest = TRUE, 
                         right = FALSE,  # Left-closed intervals [a, b)
                         labels = c(paste0("[", c(0,25,50,75), "-", c(25,50,75,100), ")"), "[100]"))]
# datalong[, newvalue:=cut(value, 
#                          breaks = c(seq(0, 90, by = 10), 100,100 + 1),  # Add 101 to include 100 in the last range                         include.lowest = TRUE, 
#                          right = FALSE,  # Left-closed intervals [a, b)
#                          labels = c(paste0("[", seq(0, 80, by = 10), "-", seq(10, 90, by = 10), ")"), "[90-100)", "[100]"))]

         
datalong[, count:=.N, by=.(newvalue, variable, tr_category.PG, fastq, ref)]
datalong[, total:=.N, by=.(variable, tr_category.PG, fastq, ref)]

datalong[, percent:=count/total*100, by=.(newvalue, variable, tr_category.PG, fastq, ref)]
datalong[, sdcount:=sd(count), by=.(newvalue, variable, tr_category.PG, ref)]
datalong[, meancount:=mean(count), by=.(newvalue, variable, tr_category.PG, ref)]
datalong[, sd:=sd(percent), by=.(newvalue, variable, tr_category.PG, ref)]
datalong[, mean:=mean(percent), by=.(newvalue, variable, tr_category.PG, ref)]

# How many transcripts do not have a 100% of reads mapped GRCh38
nonredundant <-unique(datalong[ref=="hg38", .(newvalue, meancount, ref,variable,tr_category.PG)])
nonredundant[, totalmeanreads:=sum(meancount),by="variable"]
nonredundant[, not100percentage:=sum(meancount[newvalue!="[100]"]),by="variable"]
nonredundant[, not100percentagepercentage:=not100percentage/totalmeanreads*100]


ggplot(unique(datalong[ref=="hg38", .(newvalue, mean, ref,variable,tr_category.PG, sd )]), aes(x=newvalue, y=mean, fill=tr_category.PG))+
  geom_col(position=position_dodge(0.85, preserve = "single"))+
  facet_wrap(~variable)+
  mytheme+
  labs(x="% of Reads Supporting a Transcript in non-Reference Regions from Personal Assemblies", y="% Transcripts from each Category", fill="Relative location\nof transcript\nto GRCh38")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
                    labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 0.1)  # Increase legend square size
    )
  )+
  geom_errorbar(width=0.6, aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(0.85, preserve = "single"))
ggsave("10_figures/01_plots/main/fig_05/barplot_reads_supoorting_intersecting_nonref_inpercent.pdf", dpi=700, width = 5, height = 2.75,  units = "in")

ggplot(unique(datalong[, .(newvalue, mean, ref,variable,tr_category.PG, sd )]), aes(x=newvalue, y=mean, fill=tr_category.PG))+
  geom_col(position=position_dodge(0.85, preserve = "single"))+
  facet_grid(ref~variable, labeller = labeller(ref=c("hg38"="GRCh38")))+
  mytheme+
  labs(x="% of Reads Supporting a Transcript in non-Reference Regions from Personal Assemblies",  y="% Transcripts from each Category", fill="Relative location\nof transcript\nto Reference")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
                    labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 0.1)  # Increase legend square size
    )
  )+
  geom_errorbar(width=0.6, aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(0.85, preserve = "single"))
ggsave("10_figures/01_plots/supp/38_map_pg/barplot_reads_supoorting_intersecting_nonref_inpercentage.pdf", dpi=700, width = 5, height = 4,  units = "in")


#### Now t2t
# How many transcripts do not have a 100% of reads mapped GRCh38
nonredundant <-unique(datalong[ref=="T2T", .(newvalue, meancount, ref,variable,tr_category.PG)])
nonredundant[, totalmeanreads:=sum(meancount),by="variable"]
nonredundant[, not100percentage:=sum(meancount[newvalue!="[100]"]),by="variable"]
nonredundant[, not100percentagepercentage:=not100percentage/totalmeanreads*100]


ggplot(unique(datalong[ref=="T2T", .(newvalue, mean, ref,variable,tr_category.PG, sd )]), aes(x=newvalue, y=mean, fill=tr_category.PG))+
  geom_col(position=position_dodge(0.85, preserve = "single"))+
  facet_wrap(~variable, labeller = labeller(variable=c("Reads Mapping to Reference"="Reads Mapping to T2T",
                                                       "Reads Supporting Transcripts\nin Reference"="Reads Supporting Transcripts\nin T2T")))+
  mytheme+
  labs(x="% of Reads Supporting a Transcript in non-Reference Regions from Personal Assemblies", y="% Transcripts from each Category", fill="Relative location\nof transcript\nto T2T")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_fill_manual(values=c("#9D1616","#9A4B4B", "#D18A8A"), 
                    labels=c("included"="Fully\nIncluded", "intersect"="Exonic\nIntersect", "intronic"="Only Intronic\nIntersect"))+
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA, size = 0.1)  # Increase legend square size
    )
  )+
  geom_errorbar(width=0.6, aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(0.85, preserve = "single"))
ggsave("10_figures/01_plots/supp/38_map_pg/barplot_reads_supoorting_intersecting_nonref_inpercent.pdf", dpi=700, width = 5, height = 2.75,  units = "in")




#### 
data <- fread("10_figures/data/fabien_data/blastIsoquant.txt")
data <- data[, .SD, .SDcols = c("genome", "ref",grep("quant",colnames(data), value=T))]
datalong <- melt(data, id.vars=c("genome", "ref"))
ggplot(datalong[ref=="hg38"], aes(x=variable, y=value))+
  geom_boxplot()+
  labs(y="Sequence Identity Threshold(blast)\nbetween Transcript discovered in PG and GRCh38", x="% of Transcripts")+
  mythemen+
  coord_cartesian(ylim=c(0,100))+
  scale_x_discrete(labels=c("quant_0"="100%","quant_0.001"="99.9%", "quant_0.01"="99%",  "quant_0.1"="90%" ), 
                   limits=rev(c("quant_0","quant_0.001", "quant_0.01",  "quant_0.1" )))
