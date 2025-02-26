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

data <- fread("../novelannotations/analysis_tables/250224_protein_table_for_alluvial.tsv")
# data <- fread("../novelannotations/analysis_tables/241108_n_t_structural_category_protein_category_subcategory.tsv")
# ourcat <- fread("../novelannotations/analysis_tables/241124_long_struct_cat_aa_cat.tsv")
# # orfs <- fread("../novelannotations/analysis_tables/241108_n_t_structural_category_protein_category_annot_aa.tsv")
# # orfs[, V1:=NULL]
# library(ggalluvial)
# 
# categories <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense")
# names(categories) <- unique(data$protein_splice_category)[c(2,5,7,8,6,4,3,1)]
# colsqanti <- c("#61814B", "#8EDE95", "#356CA1", "#C8773C",  "darkred", "#B5B5B5", "#4F4F4F", "#6E5353")
# names(colsqanti) <- c("FSM", "ISM", "NIC", "NNC", "Intergenic","Genic", "Fusion", "Antisense")
# ourcat[, structural_category_basic:=fifelse(structural_category%in%c("FSM", "FSM w/o CDS"), "FSM", structural_category)]
# ourcat[, structural_category_basic_novel:=fifelse(structural_category%in%c("NIC", "NNC"), "Novel", structural_category_basic)]
# 
# # Compute %
# prop.table(table(ourcat$structural_category_basic, ourcat$aa_seq_novelty), margin = 1)
# prop.table(table(ourcat$structural_category_basic_novel, ourcat$aa_seq_novelty), margin = 1)
# 
# # prepare data
# newdata <- melt(ourcat[, .(isoform, structural_category, aa_seq_novelty)], id.vars = "isoform")
# 
# newdata[, value:=fifelse(value=="Truncation", "Known\nTruncated", value)]
# newdata[, value:=fifelse(value=="FSM w/o CDS", "FSM\nwithout\nCDS", value)]
# newdata[, variable:=factor(variable, levels=rev(c("structural_category", "aa_seq_novelty" )))]
# ggplot(newdata,
#        aes(x=variable, stratum=value, fill=value, label=value, alluvium=isoform)) +
#   geom_flow(aes(order=value),width = 0.5) +
#   geom_stratum(color=NA,width = 0.5) +
#   geom_text(stat = "stratum", size = 6*0.35, family="Helvetica") +
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
#   labs(x = "", y = "# Transcripts - ORF pairs")+
#   scale_x_discrete(labels = c("structural_category"="Transcript", "aa_seq_novelty"="Predicted\nORF"), expand = c(0, 0)) +
#   scale_y_continuous( expand = c(0, 0)) +
#   guides(fill = "none")+
#   scale_fill_manual(values = c(colsqanti, "FSM\nwithout\nCDS"="#82A76A", "Known"="#61814B", "Known\nTruncated"="#94C595", "Novel"="#BC7EBF"))+
#   coord_flip()
# ggsave("10_figures/01_plots/main/fig_02/alluvial_fromTrxtoOurORFcategory.pdf", dpi=700, width = 3, height = 2,  units = "in")


#####------------------------
# **Expand the dataset to create one row per instance**
df_expanded <- data %>%
  uncount(isoform)  # Expands each row by its count
df_expanded$isoform <- paste0("iso", 1:nrow(df_expanded))
# Reshape into long format for plotting
df_long <- df_expanded %>%
  pivot_longer(cols = c(structural_category, aa_novelty_2), 
               names_to = "variable", values_to = "value")

setDT(df_long)
df_long[, value:=fifelse(value=="Known truncation", "Known\nTruncated",
                         fifelse(value=="Known elongation", "Known\nElongated", value))]
df_long[, value:=factor(value, levels=c("NNC", "NIC","Novel", "Known\nElongated", "Known", "Known\nTruncated", "NMD"))]

df_long <- df_long %>%
  group_by(variable, value) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(percentage = count / sum(count) * 100) # prepare colors
colsqanti <- c("#C8773C","#356CA1", "#BC7EBF","#738062", "#61814B", "#8EDE95", "darkgrey")
names(colsqanti) <- c("NNC", "NIC","Novel", "Known\nElongated", "Known", "Known\nTruncated", "NMD")

data[, total:=sum(isoform), by="structural_category"]
data[, percentage:=isoform/total*100, by="structural_category"]
# Create the plot using geom_flow
ggplot(df_long, 
       aes(x = variable, stratum = value, alluvium = isoform, fill = value, label = value)) +
  geom_flow(aes(order = value), width = 0.5, alpha = 0.7) +  # Flows colored by stratum values
  geom_stratum(color = NA, width = 0.5) +  # Stratum without borders
  geom_text(stat = "stratum", size = 6*0.35, family = "Helvetica") +  # Add text labels
  labs(x = "", y = "# Transcripts - ORF pairs")+
  scale_x_discrete(labels = c("structural_category"="Transcript", "aa_novelty_2"="Predicted\nORF"), expand = c(0, 0)) +
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
  scale_y_continuous( expand = c(0, 0))+
  guides(fill = "none")+
  scale_fill_manual(values = colsqanti)+
  coord_flip()
ggsave("10_figures/01_plots/main/fig_02/alluvial_fromTrxtoOurORFcategory.pdf", dpi=700, width = 3, height = 2,  units = "in")



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


