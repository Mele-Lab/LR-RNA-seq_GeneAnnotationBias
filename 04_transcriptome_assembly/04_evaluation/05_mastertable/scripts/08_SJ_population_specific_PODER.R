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
sj <- fread("02_sqanti/data/poder/poder_junctions.txt")
metadata <- fread("../../00_metadata/data/pantranscriptome_samples_metadata.tsv")
metadata <- metadata[mixed_samples==F]
metadata <- metadata[merged_run_mode==T]
popcol <- unique(metadata$color_pop)
names(popcol) <- unique(metadata$population)

# merge data
sj <- data[sj, on="isoform"]

sj[, junction:=paste(chrom, strand, genomic_start_coord,genomic_end_coord, sep="_")]
sj[, sj_category := fifelse(junction_category == "known", "Known SJ",
                            fifelse(start_site_category == "known" & end_site_category == "known", "Novel SJ: Known SS",
                                    fifelse(xor(start_site_category == "known", end_site_category == "known"), 
                                            "Novel SJ: 1 Novel SS", 
                                            "Novel SJ: 2 Novel SS")))]


## EXPLORE POP SP SJ
# Identify columns showing in which samples a transcript is found
pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
pattern <- paste(pops, collapse = "|")
samplecols <- colnames(sj)[grepl(pattern, colnames(sj))]
sj2 <- sj[, .SD, .SDcols=c(samplecols, "sj_category", "isoform", "junction")]

sjlong <- melt(sj2, measure.vars = samplecols, variable.name = "sample", value.name = "detected")

## Find pop sp SJ
# In which samples is a junction detected
sjlong[, has_detected := any(detected > 0), by = .(junction, sample)]
sjlong <- sjlong[has_detected==TRUE]
sjlong[, population := gsub(".$", "", sample)]
sjlong2 <- unique(sjlong[, .(junction, sample, population, has_detected, sj_category)])
# in how many pops is a junction detected
sjlong2[, nb.pops:=uniqueN(population), by=.(junction)]
sjlong3 <- sjlong2[nb.pops==1]
# In how many samples is a junction detected
sjlong3[, samples_per_pop:=uniqueN(sample), by=.(junction, population)]
sjlongpopsp <- unique(sjlong3[samples_per_pop>1][, .(junction, sj_category, population)])


#Calculate the proportion of "Known SJ" per population
known_proportion <- sjlongpopsp[, .(known_percentage = sum(sj_category == "Known SJ")/.N*100), by = population]
sjlongpopsp[, juncperpop := uniqueN(junction), by = population]
sjlongpopsp[, juncperpoppercat := uniqueN(junction), by = .(population, sj_category)]
sjlongpopsp[, proportion := round(juncperpoppercat/juncperpop, digits=3)*100]

# Reorder population by known_percentage
sjlongpopsp[, population := factor(population, levels = known_proportion[order(-known_percentage)]$population)]
sjlongpopsp[, sj_category:=factor(sj_category, levels = c("Known SJ", "Novel SJ: Known SS", "Novel SJ: 1 Novel SS", "Novel SJ: 2 Novel SS"))]
sjlongpopsp[, eur:=fifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]

p1<-ggplot(unique(sjlongpopsp[, .(junction,sj_category, population, eur,proportion)]), aes(x=population, fill=eur))+
  geom_bar(position="fill", aes(alpha=sj_category))+
  labs(x="", y="Proportion of Population-Specific\nSplice Junctions", fill="")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  scale_alpha_manual(values=rev(c(1,0.65, 0.4, 0.25)),name="",
                    labels = c("Known SJ"="Known SJ", "Novel SJ: Known SS"="Novel SJ:\nKnown SS", "Novel SJ: 1 Novel SS"="Novel SJ:\n1 Novel SS", "Novel SJ: 2 Novel SS"="Novel SJ:\n2 Novel SS"))+
  geom_text(data=sjlongpopsp[!sj_category%in%c( "Novel SJ: Known SS","Novel SJ: 1 Novel SS","Novel SJ: 2 Novel SS")],
    aes(label = paste0(proportion, "%"), color=sj_category),
    stat = "count",
    position = position_fill(vjust = 0.5),
    size=6*0.35 )+
  scale_color_manual(values=rep("black", 4))+
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
    alpha = guide_legend(order = 1),   # Put alpha first
    fill = guide_legend(order = 2,nrow = 2)
  )+
  scale_y_continuous(expand =  c(0, 0, 0, 0))

p2<-ggplot(unique(sjlongpopsp[, .(junction,sj_category, population, eur,proportion)]), aes(x=population, fill=eur))+
  geom_bar()+
  mytheme+
  coord_flip()+
  xlab("")+
  guides(fill="none")+
  scale_fill_manual(values=c("#466995", "#A53860"), name="")+
  ylab("# SJ")+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),  plot.margin = margin(0, 0, 5, -10),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  scale_y_continuous(expand =  c(0, 0, 0.1, 0), n.breaks=3)


### COMPUTE ODDS RATIO
sjlongpopspsub <- sjlongpopsp[, .(junction, population)][, popsp:="popsp"]
sjlong
sjlongextra <- sjlongpopspsub[unique(sjlong[, .(sj_category, junction, population)]), on=c("junction", "population")]
sjlongextra[is.na(popsp), popsp:="nonPopsp"]
sjlongextra[, new_sj_category:=fifelse(sj_category=="Known SJ", "Known SJ", "Novel SJ")]
counts <- sjlongextra[, .N, by=.(new_sj_category, population,popsp)]

# Run Fisher's exact test by population
results <- lapply(unique(counts$population), function(pop) {
  # Filter data for the population
  pop_data <- counts[population == pop]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[popsp == "nonPopsp" & new_sj_category == "Known SJ", N],
      pop_data[popsp == "nonPopsp" & new_sj_category == "Novel SJ", N],
      pop_data[popsp == "popsp" & new_sj_category == "Known SJ", N],
      pop_data[popsp == "popsp" & new_sj_category == "Novel SJ", N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("Known SJ", "Novel SJ"),
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
resultss[, population:=factor(population, levels=levels(sjlongpopsp$population))]
p3 <-ggplot(resultss, aes(x=odds.ratio, y=population, color=eur))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5)+
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
  scale_color_manual(values=c("#466995", "#A53860"), name="")+
  scale_alpha_manual(values=rev(c(0.45, 1)), name="")+
  mytheme+
  labs(x="Odds Ratio", y="")+
  theme(legend.position = "top",  plot.margin = margin(0, 0, 5, 0))+
  guides(color="none")

#### PLOT
library(patchwork)
# Ensure both p1 and p2 have the same theme
p1 <- p1 + theme(legend.position = "top")

# Combine the plots using patchwork
p3+p1+p2 +
  plot_layout(guides="collect",  widths=c(3,5.5, 0.75))&
  theme(legend.position='top')
ggsave("../../10_figures/01_plots/main/fig_03/barplot_PODER_popSpecificTrx_SJ.bycategory.pdf", dpi=700, width = 3, height = 2.25,  units = "in")


## now only the odds ratio in another order
resultss[, population:=factor(population, levels=rev(c("CEU", "AJI", "LWK", "ITU", "YRI", "PEL", "MPC","HAC")))]
ggplot(resultss, aes(x=odds.ratio, y=population, color=eur))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5)+
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
  scale_color_manual(values=c("#466995", "#A53860"), name="")+
  scale_alpha_manual(values=rev(c(0.45, 1)), name="")+
  mytheme+
  labs(x="Odds Ratio\n(Splice Junctions)", y="")+
  theme(legend.position = "top")+
  guides(color="none", alpha="none")
ggsave("../../10_figures/01_plots/main/fig_03/barplot_PODER_popSpecificTrx_SJ.bycategory_onlyOR.pdf", dpi=700, width = 1.25, height = 1.75,  units = "in")






#### now check enrichment at transcript level
n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  round(median(x), 0)))
}
poptrx <- data[, popsp:=fifelse(sample_sharing>1 & population_sharing==1, "popsp", "nonPopsp")]
ggplot(unique(poptrx[, .(isoform, length, popsp)]), aes(y=length, x=popsp, fill=popsp))+
  geom_violin(alpha=.5)+
  ggpubr::stat_compare_means(comparisons = list(c("nonPopsp", "popsp")), method = "wilcox.test", method.args = list(alternative="two.sided"), size=5*0.35)+
  geom_boxplot(outliers = F, show.legend = F, width=0.1)+
  mytheme+
  labs(x="Transcript", y="Length (nt)")+
  guides(fill="none")+
  # stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = -5000), size=5*0.35) +
  # stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = -8000), size=6*0.35)+
  scale_fill_manual(values=c("#c44536", "#457b9d"))+
  scale_x_discrete(labels=c("nonPopsp"="Population\nShared", "popsp"="Population\nSpecific"))+
  scale_y_continuous(trans="log10")
ggsave("../../10_figures/01_plots/supp/16_popsp_description/violin_popsp_length.pdf", dpi=700, width = 3, height = 2.25,  units = "in")


poptrx <- poptrx[, .(structural_category, isoform, AJI, CEU, MPC, YRI, LWK, PEL, HAC, ITU, popsp)]
poptrxlong <- melt(poptrx, id.vars=c("structural_category", "isoform", "popsp"), variable.name = "population", value.name = "nbsamples")
poptrxlong <- poptrxlong[nbsamples>=1]
poptrxlong[, trx_category:=fifelse(structural_category=="FSM", "known", "novel")]

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


#### Are genes having population specific transcripts enriched in biological functions?
library(org.Hs.eg.db)
mygenes <- unique(gsub("\\..*", "",poptrx[popsp=="popsp", geneid.v]))
allgenes <-  unique(gsub("\\..*", "",poptrx[, geneid.v]))

mygenes_converted <- clusterProfiler::bitr(mygenes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
allgenes_converted <- clusterProfiler::bitr(allgenes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


res <-clusterProfiler::enrichGO(gene = mygenes,
                          universe = allgenes,
                          keyType = "ENSEMBL",
                          ont="CC",
                          qvalueCutoff = 1,
                          OrgDb = "org.Hs.eg.db")

res <-clusterProfiler::enrichGO(gene = mygenes_converted$ENTREZID,
                                universe = allgenes_converted$ENTREZID,
                                keyType = "ENTREZID",
                                ont="BP",
                                qvalueCutoff = 1,
                                OrgDb = "org.Hs.eg.db")
################ PERMUTATION TEST --------------------------------------------------------------------------------------------------------------------------
# pops <-c("AJI.", "CEU.", "ITU.", "HAC.", "PEL.", "MPC.", "YRI.", "LWK.")
# pattern <- paste(pops, collapse = "|")
# samplecols <- colnames(data)[grepl(pattern, colnames(data))]

samplingvec <- gsub(".$", "", unique(sjlong$sample))


iterations <- 600
start_time <- Sys.time()
mylist <- list()  # Pre-allocate vector for efficiency
samplingvec_len <- length(samplingvec)

for(iteration in 1:iterations){
  print(iteration)
  # Efficient sampling and mapping population names using vectorized operations
  randompop <- setNames(sample(samplingvec, 43, replace = FALSE), unique(sjlong$sample))
  
  # Modify population directly
  sjlong[, population := randompop[sample]]
  
  ## Find pop sp SJ
  # In which samples is a junction detected
  sjlong[, has_detected := any(detected > 0), by = .(junction, sample)]
  sjlong <- sjlong[has_detected==TRUE]
  sjlong2 <- unique(sjlong[, .(junction, sample, population, has_detected, sj_category)])
  # in how many pops is a junction detected
  sjlong2[, nb.pops:=uniqueN(population), by=.(junction)]
  sjlong3 <- sjlong2[nb.pops==1]
  # In how many samples is a junction detected
  sjlong3[, samples_per_pop:=uniqueN(sample), by=.(junction, population)]
  sjlongpopsp <- unique(sjlong3[samples_per_pop>1][, .(junction, sj_category, population)])
  known_proportion <- sjlongpopsp[, .(known_percentage = sum(sj_category == "Known SJ")/.N*100), by = population]
  
  
  # Direct assignment
  mylist <- append(mylist, list(known_proportion))
}
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

fwrite(rbindlist(mylist), "05_mastertable/data/PODER_sj_permutation.propKnownSJ_in_popSpSJ.all_populations.tsv")
bootstrap <- fread("05_mastertable/data/PODER_sj_permutation.propKnownSJ_in_popSpSJ.all_populations.tsv")
### Check results
bootstrap[, iteration:=factor(rep(1:600, 8)[order(rep(1:600, 8))])]
colnames(known_proportion) <- c("population", "real")
colnames(bootstrap) <- c("population", "bootstrap","iteration")
bootstrap <- known_proportion[bootstrap, on="population"]
bootstrap[, EUR:=fifelse(population%in%c("AJI", "CEU"), "EUR", "nonEUR")]

bootstrap[, boot_eur_mean:=mean(bootstrap), by=.(EUR, iteration)]
bootstrap[, real_eur_mean:=mean(real), by=EUR]
bootstrap[, eur_mean_higher:=real_eur_mean>boot_eur_mean]

prop.table(table(bootstrap$population, bootstrap$eur_mean_higher), margin=1)









fwrite(as.data.frame(unlist(resvecpos)), "05_mastertable/data/PODER_sj_permutation.propKnownSJ_in_popSpSJ.tsv")
emppval <- sum(as.numeric(unlist(resvecpos))>0.949)/1000
emppval <- sum(as.numeric(unlist(resvecpos))>0.963)/1000
ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram(bins=30, color="black", fill="darkgrey")+
  geom_vline(xintercept=0.949, linetype="dashed", color="darkred")+
  mytheme+
  annotate(geom="text",label=paste0("Empirical      \nP-value= ", emppval), x=0.949, y=75, hjust=-0.1)+
  labs(x="Proportion of Known SJ in pseudoPopulation-Specific SJ", y="Permutations Count")
ggplot(data.frame("permutations"=as.numeric(unlist(resvecpos))), aes(x=permutations))+
  geom_histogram(bins=30, color="black", fill="darkgrey")+
  geom_vline(xintercept=0.963, linetype="dashed", color="darkred")+
  mytheme+
  annotate(geom="text",label=paste0("Empirical      \nP-value= ", emppval), x=0.963, y=75, hjust=1.1)+
  labs(x="Proportion of Known SJ in pseudoPopulation-Specific SJ", y="Permutations Count")
ggsave("../../10_figures/suppfig/barplot.permutation_knownSJpopspecific.pdf", dpi=700, width = 14, height = 11,  units = "cm")





# POPULATION SPECIFIC TRANSCRIPTS
popsp <-data[sample_sharing>=2 & population_sharing==1, ]
cols_to_keep <- colnames(popsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(popsp))]
subpopsp <- popsp[, c(cols_to_keep, "sample_sharing", "structural_category", "isoform"), with = FALSE]
subpopsplong <- melt(subpopsp, measure.vars= colnames(subpopsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(subpopsp))], variable.name = "sample", value.name = "detected")
subpopsplong <- subpopsplong[!is.na(detected)]
subpopsplong <- subpopsplong[detected==1]

subpopsplongmeta <- metadata[, .(sample, map_reads_assemblymap,population )][, total_throughput:=sum(map_reads_assemblymap), by="population"][subpopsplong, on="sample"]
subpopsplongmeta[, trx_per_cat_per_pop := uniqueN(isoform), by=c("population", "structural_category")]
subpopsplongmeta[, trx_per_cat_per_pop_norm :=trx_per_cat_per_pop/total_throughput*10^6]
subpopsplongmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]
fwrite(subpopsplongmeta, "data/241128_PODER_pop_specific_transcripts.tsv", sep="\t", quote = F, row.names = F)
# biotypes of pop specific
subpopsplongmeta_biotypes <- data[, .(isoform, associated_gene_biotype_sub)][subpopsplongmeta, on="isoform"]


ggplot(unique(subpopsplongmeta_biotypes[, .(isoform, associated_gene_biotype_sub, structural_category, population)]),
       aes(x=structural_category, fill=associated_gene_biotype_sub))+
  geom_bar()+
  scale_fill_manual(values=c("#4F4F4F","#ee9b00","darkgrey", "#4C8C36"),
                    labels=c("protein_coding"="Protein Coding", "novel/ambiguous gene"="Novel/Ambiguous Gene"))+
  labs(x="", y="# PODER Population-Specific Transcripts", fill="Gene Biotype")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  mytheme+
  geom_text(aes(label=after_stat(count)), stat="count",position = position_stack(vjust = 0.5), size=6*0.35)+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 3, -10, -7),
        legend.position = c(0.75, 0.75))
ggsave("../../../10_figures/01_plots/supp/16_popsp_val_gencode/barplot_PODER_PerSqantiCategory_popSpecificTrx_biotypes.pdf", dpi=700, width = 3, height = 2.25,  units = "in")

premergepopsp <- unique(subpopsplongmeta_biotypes[, .(isoform, associated_gene_biotype_sub)])[, popsp:=TRUE]

datapopsp <- premergepopsp[data, on="isoform"]
datapopsp[is.na(popsp), popsp:=FALSE]
table(datapopsp[, .(isoform, popsp, associated_gene_biotype)]$popsp, datapopsp[, .(isoform, popsp, associated_gene_biotype)]$associated_gene_biotype)
chires <- chisq.test(table(datapopsp[, .(isoform, popsp, associated_gene_biotype)]$popsp, datapopsp[, .(isoform, popsp, associated_gene_biotype)]$associated_gene_biotype))



# LINE PLOT FOR 2 SAMPLES SHARING FACETED
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(data=unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC")&eur=="Non-European", 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
            linewidth=1, lty="11", color="darkgrey")+
  geom_line(aes(col=structural_category),linewidth=1)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("FSM"),
                                            .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
               color = "darkgrey") +
  geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"),
                                            .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop - 25, trx_per_cat_per_pop + 25)),
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("FSM"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop + 35, trx_per_cat_per_pop - 35),
                label = population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop - 35, trx_per_cat_per_pop + 35),
                label = population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.90))+
  xlim(c(52, 105))+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(4,1.5))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(alpha="", size="", color="")
ggsave("../../../10_figures/01_plots/main/fig_03/line_PODER_popSpecific_transcripts.2samplesSharing.faceted.pdf", dpi=700, width = 4.5, height = 3,  units = "in")


# LINE PLOT FOR 2 SAMPLES SHARING FACETED with ----- REFSEQ ----------
refseq <- fread("../02_sqanti/data/poder_evaluatedBy_RefSeq/poder_evaluatedBy_RefSeq_classification.txt")
refseq <- refseq[, .(isoform, structural_category)]
categories <-c("Antisense","Intergenic", "NNC", "FSM", "NIC", "ISM", "Genic", "Fusion" )
names(categories) <- unique(refseq$structural_category)
refseq[, structural_category_refseq:=categories[structural_category]]
subpopsplongmeta_refseq <- refseq[, .(isoform, structural_category_refseq)][subpopsplongmeta, on="isoform"]
subpopsplongmeta_refseq[, trx_per_cat_per_pop:=uniqueN(isoform), by=c("structural_category_refseq", "population")]
ggplot(unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM", "NIC", "NNC"),][,.(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM", "NIC", "NNC")&eur=="Non-European", 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]), 
            linewidth=1, lty="11", color="darkgrey")+
  geom_line(aes(col=structural_category_refseq),linewidth=1)+
  mytheme+
  scale_color_manual(values=c(colsqanti, "ISM"="#8EDE95"))+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  facet_wrap(~structural_category_refseq)+
  labs(col="", alpha="", size="")+
  geom_segment(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM"),
                                                   .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
               color = "darkgrey") +
  geom_segment(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("NIC", "NNC"),
                                                   .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop - 25, trx_per_cat_per_pop + 25)),
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("FSM", "ISM"), 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop + 35, trx_per_cat_per_pop - 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  geom_text(data=unique(subpopsplongmeta_refseq[structural_category_refseq%in%c("NIC", "NNC"), 
                                                .(population, eur, structural_category_refseq, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop - 35, trx_per_cat_per_pop + 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.90))+
  xlim(c(52, 105))+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(4,1.5))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(alpha="", size="", color="")
ggsave("../../../10_figures/01_plots/supp/18_popsp_val_refseq/line_PODER_popSpecific_transcripts.2samplesSharing.faceted_refseq_PODER.pdf", dpi=700, width = 7, height = 4,  units = "in")






# LINE PLOT FOR 2 SAMPLES SHARING
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(col="", alpha="", size="")+
  geom_segment(data = unique(subpopsplongmeta[structural_category == "FSM" ,.(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data = unique(subpopsplongmeta[structural_category == "FSM", .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 30, trx_per_cat_per_pop - 30),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 4.5,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/line_PODER_popSpecific_transcripts.2samplesSharing.pdf", dpi=700, width = 16, height = 14,  units = "cm")



# LINE PLOT FOR 3 SAMPLES SHARING 
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(aes(col=structural_category),linewidth=1.5)+
  mytheme+
  scale_color_manual(values=colsqanti)+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(6.5,3))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(col="", alpha="", size="")+
  geom_segment(data = unique(subpopsplongmeta[structural_category == "FSM" ,.(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 1, trx_per_cat_per_pop - 1)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data = unique(subpopsplongmeta[structural_category == "FSM", .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),  # Only FSM data
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European" & structural_category=="FSM", trx_per_cat_per_pop + 1.2, trx_per_cat_per_pop - 1.2),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 4.5,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = "top")
ggsave("../../../10_figures/suppfig/line_PODER_popSpecific_transcripts.3samplesSharing.pdf", dpi=700, width = 16, height = 14,  units = "cm")


#########################33-----------------------------------------------------------


myunique_popsp_trx <- unique(subpopsplongmeta[structural_category%in%c("FSM", "NNC", "NIC"), .(isoform, population, structural_category)])
myunique_popsp_trx <- unique(counts[, `:=`(mean_expression=mean(TPM, na.rm=TRUE),
                                           max_expression=max(TPM, na.rm=TRUE)), by=.(transcript_id, population)][, .(transcript_id, mean_expression, max_expression, population)])[myunique_popsp_trx, on=c(transcript_id="isoform", population="population")]


ggplot(myunique_popsp_trx, aes(x=population, y=mean_expression, fill=population))+
  geom_boxplot(outliers=F)+
  facet_wrap(~structural_category)+
  mytheme+
  scale_fill_manual(values=popcol)+
  labs(x="", y="Mean Expression of Population Specific Transcrip\nacross Pop samples")+
  guides(fill="none")
ggplot(myunique_popsp_trx, aes(x=population, y=max_expression, fill=population))+
  geom_boxplot(outliers=F)+
  facet_wrap(~structural_category)+
  mytheme+
  scale_fill_manual(values=popcol)+
  labs(x="", y="Max Expression of Population Specific Transcrip\nacross Pop samples")+
  guides(fill="none")






############## FIND EUROPEAN SPECIFIC TRANSCIRPTS

# POPULATION SPECIFIC TRANSCRIPTS

prepopsp <-data[sample_sharing>=4, ]
prepopsp[, number_eur_discovered := rowSums(.SD), .SDcols = c("AJI", "CEU")][, number_noneur_discovered := rowSums(.SD), .SDcols = c("ITU", "HAC", "PEL", "MPC", "YRI", "LWK")]
popsp <- prepopsp[, eurspecificity :=fifelse(number_eur_discovered>=4 & number_noneur_discovered==0, "European-Specific",
                                             fifelse(number_eur_discovered==0 & number_noneur_discovered>=11, "non-European Specific", "Shared"))]
popsp <- popsp[eurspecificity!="Shared"]
popsp[, noveltrx:=ifelse(structural_category=="FSM", "known", "novel")]
table(popsp$noveltrx, popsp$eurspecificity)
fisher.test(table(popsp$noveltrx, popsp$eurspecificity))


cols_to_keep <- colnames(popsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(popsp))]
subpopsp <- popsp[, c(cols_to_keep, "sample_sharing", "structural_category", "isoform"), with = FALSE]
subpopsplong <- melt(subpopsp, measure.vars= colnames(subpopsp)[grepl("^[a-zA-Z]{3}[0-9]$", colnames(subpopsp))], variable.name = "sample", value.name = "detected")
subpopsplong <- subpopsplong[!is.na(detected)]
subpopsplong <- subpopsplong[detected==1]

subpopsplongmeta <- metadata[, .(sample, map_reads_assemblymap,population )][, total_throughput:=sum(map_reads_assemblymap), by="population"][subpopsplong, on="sample"]
subpopsplongmeta[, trx_per_cat_per_pop := uniqueN(isoform), by=c("eur", "structural_category")]
subpopsplongmeta[, trx_per_cat_per_pop_norm :=trx_per_cat_per_pop/total_throughput*10^6]
subpopsplongmeta[, eur := ifelse(population%in%c("CEU", "AJI"), "European", "Non-European")]

# LINE PLOT FOR 2 SAMPLES SHARING FACETED
ggplot(unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC"), 
                               .( eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
       aes(x=total_throughput/10^6, y=trx_per_cat_per_pop))+
  geom_line(data=unique(subpopsplongmeta[structural_category%in%c("FSM", "ISM", "NIC", "NNC")&eur=="Non-European", 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]), 
            linewidth=1, lty="11", color="darkgrey")+
  geom_line(aes(col=structural_category),linewidth=1)+
  mytheme+
  scale_color_manual(values=c(colsqanti))+
  labs(x="Total Mapped Reads per Population (M)", y="# Population Specific PODER Transcripts")+
  guides(color="none")+
  facet_wrap(~structural_category)+
  labs(col="", alpha="", size="")+
  geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("FSM"), 
                                            .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop + 25, trx_per_cat_per_pop - 25)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  geom_segment(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"), 
                                            .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
               aes(xend = total_throughput/10^6, # Adjust for arrow length
                   yend = ifelse(eur == "European", trx_per_cat_per_pop - 25, trx_per_cat_per_pop + 25)),
               arrow = arrow(length = unit(0.2, "cm")), # Arrow size
               color = "darkgrey") +
  ggnewscale::new_scale_color()+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("FSM"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop + 35, trx_per_cat_per_pop - 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  geom_text(data=unique(subpopsplongmeta[structural_category%in%c("NIC", "NNC"), 
                                         .(population, eur, structural_category, total_throughput,trx_per_cat_per_pop)]),
            aes(x = total_throughput/10^6,  # Position text away from point
                y = ifelse(eur == "European", trx_per_cat_per_pop - 35, trx_per_cat_per_pop + 35),
                label = population,
                color=population),
            fontface="bold",
            family="Helvetica",
            size = 8*0.35,
            alpha=0.6)+
  scale_color_manual(values=popcol)+
  guides(color="none")+
  theme(legend.position = c(0.9, 0.90))+
  xlim(c(52, 105))+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=eur, size=eur, alpha=eur))+
  scale_color_manual(values=c("#466995", "#A53860"))+
  scale_size_manual(values=c(4,1.5))+
  scale_alpha_manual(values=c(0.75, 0.9))+
  labs(alpha="", size="", color="")
