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

data <- fread("../novelannotations/analysis_tables/241108_n_t_structural_category_protein_category_subcategory.tsv")
ourcat <- fread("../novelannotations/analysis_tables/241113_struct_cat_aa_cat.tsv")
# orfs <- fread("../novelannotations/analysis_tables/241108_n_t_structural_category_protein_category_annot_aa.tsv")
# orfs[, V1:=NULL]
library(ggalluvial)
library(ggplot2)
library(ggalluvial)

categories <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense")
names(categories) <- unique(data$protein_splice_category)[c(2,5,7,8,6,4,3,1)]
colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense")

# prepare data

expanded_table <- ourcat %>%
  uncount(weights = n_t)
setDT(expanded_table)
expanded_table[, transcript:=paste0("mock_", .I)][, V1:=NULL][, perc:=NULL][, n_total_t:=NULL]
mydata <- merge(ourcat[, .(aa_seq_novelty, structural_category, n_t)], expanded_table, by = c("structural_category", "aa_seq_novelty"))
newdata <-melt(mydata[, n_total_t:=NULL], id.vars = c("transcript"))

target <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense", "TRUE", "FALSE")
target_ordered <- newdata %>% arrange(factor(value, levels = target))
target_ordered[, value :=fifelse(value=="TRUE", "Known\nORF", 
                                 fifelse(value=="FALSE", "Novel\nORF", value))]

newdata[, value:=fifelse(value=="Truncation", "Known\nTruncated", value)]
ggplot(newdata[variable%in%c("structural_category", "aa_seq_novelty")],
       aes(x=variable, stratum=value, fill=value, label=value, alluvium=transcript)) +
  geom_flow(aes(order=value)) +
  geom_stratum() +
  geom_text(stat = "stratum", size = 6*0.35) +
  mytheme+
  theme(
    axis.title.x = element_blank(),        # Removes x-axis title
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),  # Keeps y-axis line
    panel.background = element_blank(),    # Removes plot background
    panel.border = element_blank(),        # Removes plot border
    plot.background = element_blank(),     # Removes background outside panel
    panel.grid = element_blank()           # Removes grid lines
  )+
  labs(x = "", y = "# Transcripts - ORF pairs")+
  scale_x_discrete(labels = c("structural_category"="Transcript", "aa_seq_novelty"="Predicted\nORF"), expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  guides(fill = "none")+
  scale_fill_manual(values = c(colsqanti, "Known"="#637831", "Known\nTruncated"="#94C595", "Novel"="#C2623C"))
ggsave("10_figures/fig_02/alluvial_fromTrxtoOurORFcategory.pdf", dpi=700, width = 2, height = 2.5,  units = "in")



# data[, protein_splice_category:=categories[protein_splice_category]]
# data[, protein_splice_category:=factor(protein_splice_category, levels= c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense"))]
# data[, V1:=NULL]
# orfs[, protein_splice_category:=categories[protein_splice_category]]
# orfs[,  protein_splice_category:=factor(protein_splice_category, levels= c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense"))]
# orfs$protein_splice_category <- factor(orfs$protein_splice_category,  levels= c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense"))
# 
# orfs[, annot_aa:=factor(annot_aa, levels=c(TRUE, FALSE))]
# # ggplot(data, aes(axis1 = structural_category, axis2 = protein_splice_category, y = n_transcripts)) +
# #   geom_alluvium(aes(fill = structural_category)) +
# #   geom_stratum(aes(fill = protein_splice_category), color = "black", show.legend = FALSE) +
# #   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6*0.35, color = "black", family="Helvetica") +  # Adds labels to stratum
# #   labs(x = "", y = "# Transcripts - ORF pairs") +  # Customize axis labels
# #   scale_fill_manual(values = colsqanti) +
# #   scale_x_discrete(limits = c("Transcript", "Protein"),expand=c(0,0)) +
# #   scale_y_continuous(limits=c(0,sum(data$n_transcripts)), expand=c(0,0))+
# #   mytheme+
# #   theme(
# #     axis.title.x = element_blank(),        # Removes x-axis title
# #     axis.line.x = element_blank(),
# #     axis.line.y=element_line(color="black"),# Removes x-axis line
# #     panel.background = element_blank(),    # Removes plot background
# #     panel.border = element_blank(),        # Removes plot border
# #     plot.background = element_blank(),     # Removes background outside panel
# #     panel.grid = element_blank()          # Removes grid lines
# #   ) +
# #   guides(fill = "none")
# # ggsave("10_figures/suppfig/alluvial_fromTrxtoProtein.pdf", dpi=700, width = 3, height = 3,  units = "in")
# 
# 
# 
# 
# ggplot(orfs, aes(axis1 = structural_category, axis2 = protein_splice_category, axis3 = annot_aa, y = n_transcripts)) +
#   geom_alluvium(aes(fill = structural_category, order=vectororder)) +
#   scale_fill_manual(values = colsqanti) +
#   mytheme+
#   theme(
#     axis.title.x = element_blank(),        # Removes x-axis title
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),  # Keeps y-axis line
#     panel.background = element_blank(),    # Removes plot background
#     panel.border = element_blank(),        # Removes plot border
#     plot.background = element_blank(),     # Removes background outside panel
#     panel.grid = element_blank()           # Removes grid lines
#   )+
#   labs(x = "", y = "# Transcripts - ORF pairs") +  # Customize axis labels
#   scale_x_discrete(limits = c("Transcript", "Protein", "Annotation"), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, sum(orfs$n_transcripts)), expand = c(0, 0)) +
#   guides(fill = "none")+
#   geom_stratum(aes(fill = structural_category), color = "black", show.legend = FALSE) + 
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6 * 0.35, color = "black", family = "Helvetica")
# 
# 
# # prepare data
# expanded_table <- orfs %>%
#   uncount(weights = n_transcripts)
# setDT(expanded_table)
# expanded_table[, transcript:=paste0("mock_", .I)]
# mydata <- merge(orfs, expanded_table, by = c("structural_category", "protein_splice_category", "annot_aa"))
# newdata <-melt(mydata[, n_transcripts:=NULL], id.vars = c("transcript"))
# 
# target <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense", "TRUE", "FALSE")
# target_ordered <- newdata %>% arrange(factor(value, levels = target))
# target_ordered[, value :=fifelse(value=="TRUE", "Known\nORF", 
#                                  fifelse(value=="FALSE", "Novel\nORF", value))]
# # plot
# ggplot(target_ordered,
#        aes(x=variable, stratum=value, fill=value, label=value, alluvium=transcript)) +
#   geom_flow(aes(order=value)) +
#   geom_stratum() +
#   geom_text(data=target_ordered[!value%in%c("Genic", "Intergenic", "Antisense", "Fusion")],stat = "stratum", size = 6*0.35) +
#   mytheme+
#   theme(
#     axis.title.x = element_blank(),        # Removes x-axis title
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),  # Keeps y-axis line
#     panel.background = element_blank(),    # Removes plot background
#     panel.border = element_blank(),        # Removes plot border
#     plot.background = element_blank(),     # Removes background outside panel
#     panel.grid = element_blank()           # Removes grid lines
#   )+
#   labs(x = "", y = "# Transcripts - ORF pairs") +  # Customize axis labels
#   scale_x_discrete(labels = c("structural_category"="Transcript", "protein_splice_category"="Protein", "annot_aa"="Annotated\nORF"), expand = c(0, 0)) +
#   scale_y_continuous( expand = c(0, 0)) +
#   guides(fill = "none")+
#   scale_fill_manual(values = c(colsqanti, "Known\nORF"="#457b9d", "Novel\nORF"="#c44536"))
# ggsave("10_figures/suppfig/alluvial_fromTrxtoProtein.pdf", dpi=700, width = 3.5, height = 3,  units = "in")
# 


