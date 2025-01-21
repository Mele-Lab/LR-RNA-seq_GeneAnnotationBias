## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell/07_differential_expressions"
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
# load metadata
metadata <- fread("../00_metadata/data/pantranscriptome_samples_metadata.tsv")

assemblies <-c(list.files("../../novelannotations/quantifications/personal_kallisto/"))
annotations <- c("gencode", "pantrx")

quantifications <-data.frame(expand_grid(assemblies, annotations))
quantifications$annotnamecode <- ifelse(quantifications$annotations=="gencode", "GENCODEv47",
                                        "PODER")
quantifications$name <- paste(quantifications$assemblies, quantifications$annotnamecode, sep=":")
# load data
quantilist <- list()
quantilist <- append(quantilist, list(fread("data/02_DEGres_gencode.tsv")))
quantilist <- append(quantilist, list(fread("data/02_DEGres_pantrx.tsv")))
namelist <- c("GRCh38:GENCODEv47", "GRCh38:PODER")

for(i in 1:nrow(quantifications)){
  quantilist <- append(quantilist, list(fread(paste0("data/other_assemblies/02_DEGres_", quantifications[i, "assemblies"],"_", quantifications[i, "annotations"], ".tsv"))))
  namelist <- c(namelist, quantifications[i, "name"])
  }

names(quantilist) <- namelist
data <- rbindlist(quantilist, idcol="quantification")



data[, popincontrast := factor(ifelse((reference == "CEU" & contrast == "MPC") | (reference == "MPC" & contrast == "CEU"), 
                               "CEU & MPC in Contrast",
                               ifelse(reference == "CEU" | contrast == "CEU", "CEU in Contrast",
                                      ifelse(reference == "MPC" | contrast == "MPC", "MPC in Contrast",
                                             "Other Populations"))), levels=c("CEU in Contrast", "MPC in Contrast","CEU & MPC in Contrast", "Other Populations" ))]


# Plot logFC distributions
ggplot(data[quantification=="GRCh38:GENCODEv47"], aes(y=contrast_name, x=logFC, fill=popincontrast))+
  ggridges::geom_density_ridges()+
  mytheme+
  scale_fill_manual(values=c(unique(metadata[population=="CEU", color_pop]),
                             unique(metadata[population=="MPC", color_pop]),
                             "darkgreen",
                             "darkgrey"))
ggplot(data, aes(y=contrast_name, x=logFC, fill=popincontrast))+
  ggridges::geom_density_ridges()+
  mytheme+
  facet_wrap(~quantification, nrow = 3)+
  scale_fill_manual(values=c(unique(metadata[population=="CEU", color_pop]),
                             unique(metadata[population=="MPC", color_pop]),
                             "darkgreen",
                             "darkgrey"))
# Plot intersection of DEG between quantifications
# we want to intersect each contrast across quantifications or in other words: in how many quantifications is each gene found, by contrast
datasig <-data[fdr<0.05][,quantification_sharing:=.N, by=c("contrast_name", "geneid.v")]
ggplot(datasig, aes(x=quantification_sharing, fill=contrast_name))+
  geom_bar()+
  mytheme+
  labs(y="# DEGs", x="Quantification Sharing (Assembly-Annotation pairs)", fill="")
ggplot(datasig, aes(x=quantification_sharing, fill=popincontrast))+
  geom_bar()+
  mytheme+
  scale_fill_manual(values=c(unique(metadata[population=="CEU", color_pop]),
                             unique(metadata[population=="MPC", color_pop]),
                             "darkgreen",
                             "darkgrey"))+
  labs(y="# DEGs", x="DEG Sharing (Assembly-Annotation pairs)", fill="")+
  geom_text(aes(label=after_stat(count)), stat="count", position=position_stack(vjust=0.5), size=6*0.35)
ggsave(paste0("../10_figures/01_plots/supp/28_deg_assemblies/heatmap_DGE_sharing.pdf"), dpi=700, width = 4, height = 2.75,  units = "in")
