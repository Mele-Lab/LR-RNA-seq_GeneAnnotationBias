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

# load data
genexpression <- fread("data/unmasked_genomic_wrong_index/gencode_v47/gene_tpm_matrix.tsv")
annot <- fread("/home/pclavell/mounts/mn5/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcript_parsed.tsv")
annot <- annot[, c("transcriptid", "transcriptid.v", "transcript_biotype", "transcript_name", "start", "end"):=NULL]

genelong <- melt(genexpression, value.name="tpm", variable.name="sample")[geneid.v!="__unassigned"]
uniqgene <- annot[unique(genelong[, meantpm:= mean(tpm, na.rm=T), by=geneid.v][, .(geneid.v, meantpm)]), on="geneid.v"]
uniqgene <- unique(uniqgene)


uniqgene[, pseudogene:= ifelse(grepl("processed", gene_biotype), ifelse(grepl("unprocessed", gene_biotype), "unprocessed", "pseudogene"), "not_pseudogene")]
uniqgene <- unique(uniqgene[!is.na(geneid)])
# # plot
ggplot(unique(uniqgene[meantpm>0.1, .(geneid, gene_biotype, meantpm, pseudogene)]), aes(x=gene_biotype, y=log10(meantpm), col=pseudogene))+
  geom_jitter(alpha=0.5)+
  mytheme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values=c("grey", "red", "purple"))

## look for parent-pseudogene families
# obtain parent names
pseudogenes <- uniqgene[pseudogene!="not_pseudogene" & !grepl("^ENSG", gene_name),
                        .(gene_name)][, family:=gsub("P.?.?$", "",gene_name)]

# keep only families whose parent is a gene and in particular a protein_coding gene
families_vec <- unique(pseudogenes$family)[unique(pseudogenes$family)%in%annot$gene_name]
families_vec <- unique(annot[data.table(gene_name=families_vec), on="gene_name"][gene_biotype=="protein_coding", gene_name])

# add row with parent info
protcod_rows <- data.table(gene_name= families_vec, family=families_vec)
pseudogenes <- rbind(pseudogenes, protcod_rows)
pseudogenes <- unique(pseudogenes[family %in% families_vec])
pseudofam <- unique(uniqgene[pseudogenes, on="gene_name"])
pseudofam$meantpm <- replace_na(pseudofam$meantpm, 0)
# remove lnRNA that are double annotated as pseudogens and lnRNA
pseudofam <- unique(pseudofam[gene_biotype!="lncRNA" & gene_biotype!="TEC"])

pseudofam <-unique(pseudofam[, .(family, meantpm, gene_biotype)])

# keep families with total expression higher tha 0.1 tpm
pseudofam <-pseudofam[, totaltpm:=sum(meantpm), by=family][totaltpm>0.1]
prot <- unique(pseudofam[gene_biotype=="protein_coding",  .(meantpm, family)])
prot$gene_biotype <- "protein_coding"
colnames(prot)[1:2] <- c("protpm", "family")
prot[, totalprotpm := sum(protpm), by=family][, protpm:=NULL]
# keep families where pseudogenes acount for at least 1% of the total expression
pseudofam2 <- unique(prot)[pseudofam, on=c("family", "gene_biotype")]
pseudofam2 <-pseudofam2[, protpercent := totalprotpm/totaltpm]
pseudofam2[, pseudopercent := 1-protpercent]


melted_dt <- melt(pseudofam2, id.vars = c("family", "totaltpm"), 
                  measure.vars = c("protpercent", "pseudopercent"),
                  variable.name = "percent_type", value.name = "percent_value")

melted_dt$percent_value <-replace_na(melted_dt$percent_value, 0)
melted_dt <- unique(melted_dt)



melted_dt2 <- melt(pseudofam2[protpercent<0.95,], id.vars = c("family", "totaltpm"), 
                  measure.vars = c("protpercent", "pseudopercent"),
                  variable.name = "percent_type", value.name = "percent_value")

melted_dt2$percent_value <-replace_na(melted_dt2$percent_value, 0)
melted_dt2 <- unique(melted_dt2)



ggplot(melted_dt2, aes(x=reorder(family, totaltpm), y=percent_value, fill=percent_type))+
  geom_col(pos="fill")+
  mytheme+
  scale_fill_manual(values=c("grey", "red"), 
                    labels=c("protein coding", "pseudogenes"))+
  ggnewscale::new_scale_color()+
  geom_point(data=melted_dt2,
             aes(y=-0.05, x=reorder(family, totaltpm), col=log10(totaltpm)), shape=15, size=4)+
  scale_color_gradientn(colours=c("white", "yellow", "red"))+
  ggnewscale::new_scale_color()+
  geom_point(data=melted_dt2[, RP:=ifelse(grepl("^RP", family), "RP", "notRP")],
             aes(y=-0.1, x=reorder(family, totaltpm), col=RP), shape=15, size=4)+
  scale_color_manual(values=c("white", "green"))+
  xlab("")+
  ylab("% expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(fill="")

melted_dt2[, blacklist :=family]
fwrite(melted_dt, "../00_metadata/blacklist_pseudogenes_problematic_genes.tsv", sep="\t", quote = F)
