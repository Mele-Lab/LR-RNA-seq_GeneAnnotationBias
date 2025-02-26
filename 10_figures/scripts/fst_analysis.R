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



data <- fread("../novelannotations/analysis_tables/250213_snp_exon_60bp_intersect_fst.tsv")
data <- data[, .(novelty, pop1_det, pop2_det, pop1, pop2, fst)]
data[, detection:=factor(fifelse(pop1_det & pop2_det, "both", "one"))]


result <- data[, .(
  wilcox_p_value =  wilcox.test(fst[detection == "one"], fst[detection == "both"], alternative = "greater")$p.value,
  wilcox_statistic =  wilcox.test(fst[detection == "one"], fst[detection == "both"], alternative = "greater")$statistic
), by = .(pop1, pop2)]


result2 <- data[,   mean_fst := mean(fst, na.rm=T), by=.(pop1, pop2, detection)]
result2 <-unique(result2[, .(pop1, pop2, mean_fst, detection)])
result2wide <- dcast(result2, pop1+pop2~detection, value.var = "mean_fst")
result2wide[, diff:=one-both]


result <- result2wide[, .(pop1, pop2, diff)][result, on=.(pop1, pop2)]

result[, fdr:=p.adjust(wilcox_p_value, method="BH")]
result[, sig := fifelse(fdr <= 0.001, "***",
                        fifelse(fdr <= 0.01, "**",
                                fifelse(fdr <= 0.05, "*", "")))]
result <- rbind.data.frame(result, c(rep("CEU", 2), NA, NA, NA,NA,NA))
result <- rbind.data.frame(result, c(rep("LWK", 2), NA, NA, NA,NA,NA))
result <- rbind.data.frame(result, c(rep("YRI", 2), NA, NA, NA,NA,NA))
result <- rbind.data.frame(result, c(rep("PEL", 2), NA, NA, NA,NA,NA))
result <- rbind.data.frame(result, c(rep("HAC", 2), NA, NA, NA,NA,NA))
result <- rbind.data.frame(result, c(rep("ITU", 2), NA, NA, NA,NA,NA))
result[pop1=="HAC" & pop2=="CEU", `:=`(pop1="CEU", pop2="HAC")]
result[pop1=="LWK" & pop2=="PEL", `:=`(pop1="PEL", pop2="LWK")]
result[pop1=="HAC" & pop2=="ITU", `:=`(pop1="ITU", pop2="HAC")]


# Corrected syntax
result[, c("wilcox_p_value", "wilcox_statistic", "fdr", "diff") := 
         lapply(.SD, as.numeric), 
       .SDcols = c("wilcox_p_value", "wilcox_statistic", "fdr", "diff")]
result[, c("pop1", "pop2") := 
         lapply(.SD, as.factor), 
       .SDcols = c("pop1", "pop2")]
result[, rmv:=factor(ifelse(pop1==pop2, TRUE, FALSE))]

ggplot(result, aes(x=pop2, y=pop1))+
  geom_tile(aes(fill=diff,alpha=rmv))+
  geom_text( aes(label=sig), size=6*0.35 )+
  scale_x_discrete(limits=c("CEU", "ITU", "HAC", "PEL", "LWK", "YRI"))+
  scale_y_discrete(limits=c("CEU", "ITU", "HAC", "PEL", "LWK", "YRI"))+
  scale_fill_continuous(low="#e6f2c2", high="#61814B", na.value="white")+
  labs(x="", y="", fill="Difference\nMean FST", title="")+
  geom_vline(xintercept=4.5, linetype="dashed")+
  annotate(geom="text", x=4.5,y=6, label="African", hjust=-0.2, size=6*0.35)+
  annotate(geom="text", x=4.5,y=6, label="OOA", hjust=1.2,  size=6*0.35)+
  mytheme+
  theme(legend.position = c(0.4,0.7))+
  scale_alpha_manual(values=c(1,0))+
  guides(alpha="none")+
  ggnewscale::new_scale_fill()+
  geom_tile(data=unique(result[, .(pop1, pop2)]),aes(y = 0.375, x=pop2,fill = pop2),  # Position tiles below y-axis labels
            height = 0.2)+
  geom_tile(aes(x = 0.375, y=pop1,fill = pop1),  # Position tiles below y-axis labels
            height = 1, width=0.2)+
  scale_fill_manual(values = popcol, na.value = "grey90")+
  guides(fill="none")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/heatmap_WilcoxonSignificantDifferencesbetweenPOPS_FST_60bp.pdf", dpi=700, width = 1.8, height = 2.25,  units = "in")



### perfomr enrichments
fst_threshold <- 0.2

data[, divergence :=fifelse(fst>=fst_threshold, "high", "low")]
data <- data[pop1=="CEU" | pop2=="CEU"]
data[, pop:=fifelse(pop1=="CEU", pop2,pop1)]
data <- data[, .(novelty, pop, divergence, detection)]






##### FISHER TO SEE IF POPSP EXONS ARE ENRICHED IN HIGHFST
# Run Fisher's exact test by population
results <- lapply(unlist(unique(data$pop)), function(pops) {
  # Filter data for the population
  pop_data <- data[pop == pops,]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[divergence == "low" & detection == "both", .N],
      pop_data[divergence == "low" & detection == "one", .N],
      pop_data[divergence == "high" & detection == "both", .N],
      pop_data[divergence == "high" & detection == "one", .N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("both", "one"),
      rev(c("high", "low"))
    )
  )
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, conf.int = T)
  
  list(population = pops, p.value = test_result$p.value, odds.ratio = test_result$estimate, ci=list(test_result$conf.int))
})

resultss <- rbindlist(results)
# Extract lower and upper CI values into new columns
resultss[, `:=`(ci_lower = sapply(ci, `[`, 1),  # Extract first value of CI
                ci_upper = sapply(ci, `[`, 2))] 

resultss[, fdr:=p.adjust(p.value, method = "BH")]
resultss[, FDR:=factor(fifelse(fdr<0.05, "FDR<0.05", "FDR>=0.05"))]
ggplot(resultss, aes(x=odds.ratio, y=reorder(population, odds.ratio)))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5, show.legend=F)+
  geom_vline(xintercept=1, linetype="dashed")+
  scale_color_manual(values="black")+
  scale_alpha_manual(values=rev(c(0.25, 1)))+
  mytheme+
  labs(x="Odds Ratio", y="", color="", alpha="")+
  theme(legend.position = "top")+
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha=guide_legend(override.aes = list(size = 3), nrow=2))+
  annotation_logticks(sides="b")+
  guides(color="none")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/dotplot_enrichment.popspExons_in_highfst_CEUvsallPOPS_FST_60bp.pdf", dpi=700, width = 1.8, height = 2.25,  units = "in")





##### FISHER TO SEE IF NOVEL EXONS ARE ENRICHED IN HIGHFST
# Run Fisher's exact test by population
results <- lapply(unlist(unique(data$pop)), function(pops) {
  # Filter data for the population
  pop_data <- data[pop == pops,]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[divergence == "low" & novelty == "Known", .N],
      pop_data[divergence == "low" & novelty == "Novel", .N],
      pop_data[divergence == "high" & novelty == "Known", .N],
      pop_data[divergence == "high" & novelty == "Novel", .N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("known", "novel"),
      rev(c("high", "low"))
    )
  )
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, conf.int = T)
  
  list(population = pops, p.value = test_result$p.value, odds.ratio = test_result$estimate, ci=list(test_result$conf.int))
})

resultss <- rbindlist(results)
# Extract lower and upper CI values into new columns
resultss[, `:=`(ci_lower = sapply(ci, `[`, 1),  # Extract first value of CI
                ci_upper = sapply(ci, `[`, 2))] 

resultss[, fdr:=p.adjust(p.value, method = "BH")]
resultss[, FDR:=factor(fifelse(fdr<0.05, "FDR<0.05", "FDR>=0.05"))]
ggplot(resultss, aes(x=odds.ratio, y=reorder(population, odds.ratio)))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5, show.legend=F)+
  geom_vline(xintercept=1, linetype="dashed")+
  scale_color_manual(values="black")+
  scale_alpha_manual(values=rev(c(0.25)))+
  mytheme+
  labs(x="Odds Ratio", y="", color="", alpha="")+
  theme(legend.position = "top")+
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha=guide_legend(override.aes = list(size = 3), nrow=2))+
  guides(color="none")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/dotplot_enrichment.novel_in_highfst_CEUvsallPOPS_FST_60bp.pdf", dpi=700, width = 1.8, height = 2.25,  units = "in")



##### FISHER TO SEE IF Novel 5'/3' EXONS ARE ENRICHED IN HIGHFST
# Run Fisher's exact test by population
results <- lapply(unlist(unique(data$pop)), function(pops) {
  # Filter data for the population
  pop_data <- data[pop == pops,]
  
  # Create contingency table
  contingency_table <- matrix(
    c(pop_data[divergence == "low" & novelty == "Known", .N],
      pop_data[divergence == "low" & novelty == "Novel 5'/3'", .N],
      pop_data[divergence == "high" & novelty == "Known", .N],
      pop_data[divergence == "high" & novelty == "Novel 5'/3'", .N]),
    nrow = 2, byrow = F,
    dimnames = list(
      c("known", "Novel 5'/3'"),
      rev(c("high", "low"))
    )
  )
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, conf.int = T)
  
  list(population = pops, p.value = test_result$p.value, odds.ratio = test_result$estimate, ci=list(test_result$conf.int))
})

resultss <- rbindlist(results)
# Extract lower and upper CI values into new columns
resultss[, `:=`(ci_lower = sapply(ci, `[`, 1),  # Extract first value of CI
                ci_upper = sapply(ci, `[`, 2))] 

resultss[, fdr:=p.adjust(p.value, method = "BH")]
resultss[, FDR:=factor(fifelse(fdr<0.05, "FDR<0.05", "FDR>=0.05"))]
ggplot(resultss, aes(x=odds.ratio, y=reorder(population, odds.ratio)))+
  geom_point(size=1.5, aes(alpha=FDR))+
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper, alpha=FDR), width = 0.25, linewidth=0.5, show.legend=F)+
  geom_vline(xintercept=1, linetype="dashed")+
  scale_color_manual(values="black")+
  scale_alpha_manual(values=rev(c(0.25, 1)))+
  mytheme+
  labs(x="Odds Ratio", y="", color="", alpha="")+
  theme(legend.position = "top")+
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha=guide_legend(override.aes = list(size = 3), nrow=2))+
  guides(color="none")
ggsave("10_figures/01_plots/supp/41_personalized_GRCh38/dotplot_enrichment.novel35prime_in_highfst_CEUvsallPOPS_FST_60bp.pdf", dpi=700, width = 1.8, height = 2.25,  units = "in")
