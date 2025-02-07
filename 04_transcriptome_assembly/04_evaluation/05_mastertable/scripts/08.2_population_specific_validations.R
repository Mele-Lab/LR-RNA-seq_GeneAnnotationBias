## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation"
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


## DATA PREPARATION
# load data
data <- fread("05_mastertable/data/29102024_PODER_mastertable.tsv")
data <- data[structural_category%in%c("FSM","NIC" ,"NNC")]
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- unique(metadata$color_pop)
names(popcol) <- unique(metadata$population)

poptrx <- data[, popsp:=fifelse(sample_sharing>1 & population_sharing==1, "popsp", "nonPopsp")]


refseq <- fread("../04_evaluation/02_sqanti/data/poder_evaluatedBy_RefSeq/poder_evaluatedBy_RefSeq_classification.txt")
refseq <- refseq[, .(isoform, structural_category_refseq)]
categories <-c("Antisense","Intergenic", "NNC", "FSM", "NIC", "ISM", "Genic", "Fusion" )
names(categories) <- unique(refseq$structural_category)
refseq[, structural_category_refseq:=categories[structural_category]]










poptrx <- poptrx[, .(structural_category, isoform, AJI, CEU, MPC, YRI, LWK, PEL, HAC, ITU, popsp)]
poptrxlong <- melt(poptrx, id.vars=c("structural_category", "isoform", "popsp"), variable.name = "population", value.name = "nbsamples")
poptrxlong <- poptrxlong[nbsamples>=1]

poptrxlong <-refseq[, .(isoform, structural_category_refseq)][poptrxlong, on="isoform"]


poptrxlong[, trx_category:=fifelse(structural_category_refseq%in%c("FSM", "ISM"), "known", "novel")]


### COMPUTE ODDS RATIO
counts <- poptrxlong[, .N, by=.(trx_category, population,popsp)]

# Run Fisher's exact test by population
results <- lapply(unlist(unique(counts[,population])), function(pop) {
  # Filter data for the population
  pop_data <- counts[population == pop]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[popsp == "nonPopsp" & trx_category == "known", N],
      pop_data[popsp == "nonPopsp" & trx_category == "novel", N],
      pop_data[popsp == "popsp" & trx_category == "known", N],
      pop_data[popsp == "popsp" & trx_category == "novel", N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("known", "novel"),
      rev(c("popsp", "nonPopsp"))
    )
  )
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, conf.int = T)
  
  list(population = pop, p.value = test_result$p.value, odds.ratio = test_result$estimate, ci=list(test_result$conf.int))
})

resultss <- rbindlist(results)
# Extract lower and upper CI values into new columns
resultss[, `:=`(ci_lower = sapply(ci, `[`, 1),  # Extract first value of CI
                ci_upper = sapply(ci, `[`, 2))] 

resultss[, eur:=fifelse(population%in%c("AJI", "CEU"), "EUR", "nonEUR")]
resultss[, fdr:=p.adjust(p.value, method = "BH")]
resultss[, FDR:=factor(fifelse(fdr<0.05, "FDR<0.05", "FDR>=0.05"))]
resultss[, population:=factor(population, levels=levels(poptrxlong$population))]

order_samples <-counts[popsp=="popsp", total:=sum(N), by=.(population)][, known_per:=N/total][trx_category=="known" & popsp=="popsp"][order(known_per), population]

p3 <- ggplot(resultss, aes(x=odds.ratio, y=population, color=eur))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5, show.legend=F)+
  geom_vline(xintercept=1, linetype="dashed")+
  annotate(geom="rect", 
           xmin=min(resultss[eur=="nonEUR", ci_lower]),
           xmax=max(resultss[eur=="nonEUR", ci_upper]),
           ymin=0.5,
           ymax=8.5,
           fill="#A53860",
           alpha=0.2)+
  annotate(geom="rect", 
           xmin=min(resultss[eur=="EUR", ci_lower]),
           xmax=max(resultss[eur=="EUR", ci_upper]),
           ymin=0.5,
           ymax=8.5,
           fill="#466995",
           alpha=0.2)+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_alpha_manual(values=rev(c(0.45, 1)))+
  mytheme+
  labs(x="Odds Ratio\n(log10)", y="", color="", alpha="")+
  theme(legend.position = "top")+
  scale_x_continuous(trans="log10", breaks=c(0.5, 1, 2, 3))+
  scale_y_discrete(limits=order_samples)+
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha=guide_legend(override.aes = list(size = 3), nrow=2))+
  annotation_logticks(sides="b")+
  guides(color="none")
# ggsave("../../10_figures/01_plots/main/fig_04/barplot_PODER_popSpecificTrx_enrichment.bycategory.pdf", dpi=700, width = 3, height = 2.25,  units = "in")

dt1 <- unique(counts[popsp=="popsp", total:=sum(N), by=.(population)][, known_per:=N/total][, eur:=ifelse(population%in%c("CEU", "AJI"), "European", "non-European")][popsp=="popsp"])

p1 <- ggplot(dt1, 
             aes(x=population, fill=eur, y=known_per))+
  geom_col(position="fill", aes(alpha=trx_category))+
  geom_text(data=dt1,
            aes(label = paste0(round(known_per*100, 1), "%"), group=trx_category),
            position = position_fill(vjust = 0.5),
            size=6*0.35 )+
  scale_x_discrete(limits=order_samples)+
  labs(x="", y="Proportion of\nPopulation-Specific\nTranscripts", fill="")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  scale_alpha_manual(values=c( 0.25,0.75),name="",labels = c("known"="Known", "novel"="Novel"))+
  scale_color_manual(values=rep("black", 2))+
  guides(color="none")+
  coord_flip()+
  mytheme+
  theme(legend.position = "top",
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(margin = margin(r = -1, unit = "pt")),
        legend.box.margin=margin(-10,-10,-10,-45),
        legend.box = "horizontal",               # Ensure horizontal alignment
        legend.box.just = "center",              # Center the legend box
        legend.box.spacing = unit(0.5, "cm"),    # Adjust spacing between rows
        legend.direction = "horizontal",         # Horizontal direction for multiple rows
        legend.justification = "left",
        legend.spacing.y = unit(0, "cm"),
        plot.margin = margin(0, 3, 5, -5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 3),order = 2,nrow = 2),
    alpha=guide_legend(order = 1, nrow = 2, override.aes = list(size = 3))
  )+
  scale_y_continuous(expand =  c(0, 0, 0, 0))

p2 <- ggplot(unique(poptrxlong[popsp=="popsp", .(isoform,trx_category, population)])[, eur:=fifelse(population%in%c("AJI", "CEU"), "European", "Non-European")], aes(x=population, fill=eur))+
  geom_bar()+
  mytheme+
  coord_flip()+
  xlab("")+
  guides(fill="none")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  ylab("# Transcripts")+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),  plot.margin = margin(0, 0, 5, -10),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  scale_y_continuous(expand =  c(0, 0, 0.1, 0), n.breaks=3)+
  scale_x_discrete(limits=order_samples)

library(patchwork)
# Ensure both p1 and p2 have the same theme
p1 <- p1 + theme(legend.position = "top")

# Combine the plots using patchwork
p3+p1+p2 +
  plot_layout(guides="collect",  widths=c(3,5.5, 0.75))&
  theme(legend.position='top')
ggsave("../../10_figures/01_plots/main/fig_03/barplot_PODER_popSpecificTrx_trx.bycategory.pdf", dpi=700, width = 3.25, height = 2.75,  units = "in")
