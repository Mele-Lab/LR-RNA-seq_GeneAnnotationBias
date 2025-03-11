setwd("[PATH]/fdegalez/mapping_minimap/RNAseqLR_mapping")

bed.colnames <- c("seqname.reg", "start.reg", "end.reg", "length.reg", 
                  "seqname", "start", "end", "gene_id", "tr_id", "ex_id")
genome <- "NA18906"
ref <- "hg38"

res.all <- NULL
res.collapse.all <- NULL

for (genome in c("HG002", "HG00621", "HG01928", "HG01952", "NA18906", "NA19240")){
    print(genome)
    for (ref in c("hg38", "T2T")){
        ## Intersect files
        print(ref)
        ex_intersect <- read.delim(paste0(genome,".fq_",genome,".fa/isoquant/OUT/OUT.exon_models_intersectUniqueRegions_",ref,".bed"),
                                   header = F, stringsAsFactors = F )
        ex_included <- read.delim(paste0(genome,".fq_",genome,".fa/isoquant/OUT/OUT.exon_models_includedUniqueRegions_",ref,".bed"),
                                  header = F, stringsAsFactors = F )
        tr_intersect <- read.delim(paste0(genome,".fq_",genome,".fa/isoquant/OUT/OUT.transcript_models_intersectUniqueRegions_",ref,".bed"),
                                   header = F, stringsAsFactors = F )
        tr_included <- read.delim(paste0(genome,".fq_",genome,".fa/isoquant/OUT/OUT.transcript_models_includedUniqueRegions_",ref,".bed"),
                                  header = F, stringsAsFactors = F )
        
        colnames(tr_intersect) <- bed.colnames
        colnames(tr_included) <- bed.colnames
        colnames(ex_intersect) <- bed.colnames
        colnames(ex_included) <- bed.colnames
        
        ## GTF info
        GTF_info <- read.delim(paste0(genome,".fq_",genome,".fa/isoquant/OUT/infoID.gtf.txt"),
                               header = F, stringsAsFactors = F)
        colnames(GTF_info) <- c("gene_id", "tr_id", "ex_id")
        tr_nbEx <- as.data.frame(table(GTF_info$tr_id))
        colnames(tr_nbEx) <- c("tr_id", "nb_ex")
        
        ## exon icnluded
        tr_nbEx_included <- as.data.frame(table(ex_included$tr_id))
        colnames(tr_nbEx_included) <- c("tr_id", "nb_ex_included")
        
        ## exon intersect
        tr_nbEx_intersected <- as.data.frame(table(ex_intersect$tr_id))
        colnames(tr_nbEx_intersected) <- c("tr_id", "nb_ex_intersected")
        
        res <- data.frame(tr_id = tr_intersect$tr_id, stringsAsFactors = F)
        res$isIncluded <- res$tr_id %in% tr_included$tr_id
        ##
        res <- merge(res, tr_nbEx, by = "tr_id", all.x = T)
        res <- merge(res, tr_nbEx_intersected, by = "tr_id", all.x = T)
        res <- merge(res, tr_nbEx_included, by = "tr_id", all.x = T)
        
        res$propExIncl <- res$nb_ex_included/res$nb_ex
        
        
        ## intronic part
        tr_intronic <- tr_intersect$tr_id[!(tr_intersect$tr_id %in% ex_intersect$tr_id)]
        res$category <- NULL 
        res$category[res$tr_id %in% ex_intersect$tr_id] <- "intersect"
        res$category[res$tr_id %in% tr_intronic] <- "intronic" 
        res$category[res$isIncluded == T] <- "included"
        print(table(res$category, useNA= "always"))
        res$pop <- genome
        res$unique_region <- ref
        
        res <- unique(res)
        
        
        ##### Add the size of the regions intersecting/includ. the transcript
        
        region_intersect <- tr_intersect[, c("seqname.reg", "length.reg", "tr_id")]
        region_intersect <- unique(region_intersect)
        
        ## Nb of regions Regions
        nbRegions <- aggregate(region_intersect$seqname.reg, by = list(region_intersect$tr_id), length)
        colnames(nbRegions) <- c("tr_id", "nbRegions")
        ## Name of the regions
        nameRegions <- aggregate(region_intersect$seqname.reg, by = list(region_intersect$tr_id), paste0, collapse = "/")
        colnames(nameRegions) <- c("tr_id", "nameRegions")
        ## Cum Length of the regiosn
        totalSizeRegions <- aggregate(region_intersect$length.reg, by = list(region_intersect$tr_id), sum)
        colnames(totalSizeRegions) <- c("tr_id", "totalSizeRegions")
        ## length of the regions
        sizeRegions <- aggregate(region_intersect$length.reg, by = list(region_intersect$tr_id), paste0, collapse = "/")
        colnames(sizeRegions) <- c("tr_id", "sizeRegions")
        
        
        infoRegions <- merge(nbRegions, nameRegions, by = "tr_id")
        infoRegions <- merge(infoRegions, totalSizeRegions, by = "tr_id")
        infoRegions <- merge(infoRegions, sizeRegions, by = "tr_id")
        
        res_region <- merge(res, infoRegions, by = "tr_id")
        
        
        res.collapse <- data.frame(table(res$category), stringsAsFactors = F)
        colnames(res.collapse) <- c("type", "freq")
        res.collapse$pop <- genome
        res.collapse$ref <- ref
        
        
        # write.table(res, paste0(genome,".fq_",genome,".fa/isoquant/OUT/transcriptTypeisoquant_", ref,".tsv"),
        #             quote = F, sep = "\t", col.names = T, row.names = F)
        
        
        res.all <- rbind(res.all, res_region)
        res.collapse.all <- rbind(res.collapse.all, res.collapse)
    }
}

write.table(res.all, "../isoquant_analysis/transcriptTypeisoquant.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(res.collapse.all, "../isoquant_analysis/transcriptTypeisoquant_collapse.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

