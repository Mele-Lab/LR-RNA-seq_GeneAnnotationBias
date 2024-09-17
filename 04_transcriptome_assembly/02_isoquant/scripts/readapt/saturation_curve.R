# Count of isoforms discovered by Isoquant with /home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/02_isoquant/scripts/01_isoquant_discovery_quantification.sh
discisoforms <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/02_isoquant/stats/count_discovered_transcripts.tsv")
colnames(discisoforms) <- c("sample", "isoforms")
discisoforms$sample <- unlist(lapply(strsplit(discisoforms$sample, "/"), \(x) x[[2]]))
mapq10reads <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/03_mapping/data/count_mapped_reads.tsv")
colnames(mapq10reads) <- c("sample", "fail", "genomemapped", "sirvmapped", "fastq")

discisoforms <- mapq10reads[discisoforms, on=.(sample)]

ggplot(discisoforms, aes(x=genomemapped/10^6, y=isoforms))+
  stat_poly_line(show.legend = F) +
  stat_poly_eq(use_label(c("eq", "adj.R2")), label.y=0.9)  +
  stat_poly_eq(use_label(c( "p", "n")), label.y=0.8)  +  
  geom_point()+
  guides(line=F)+
  labs(x="Mapped reads (M)", y="Discovered isoforms")+
  geom_abline(slope=1, linetype="dashed", col="grey")+
  mytheme
