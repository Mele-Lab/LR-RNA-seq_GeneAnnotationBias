setwd("[PATH]/mapping_minimap/repetitionQuantification/genome/")

# Function to calculate Shannon Entropy for a single sequence
calculate_shannon_entropy <- function(kmer_counts) {
    # Total k-mers in the sequence
    total_kmers <- sum(kmer_counts, na.rm = T)
    # Normalize the k-mer counts to probabilities
    probabilities <- kmer_counts / total_kmers
    # Calculate Shannon entropy
    entropy <- -sum(probabilities * log2(probabilities), na.rm = TRUE)
    return(entropy)
}


dir.list <- list.dirs(".", recursive = F)
dir.list <- str_remove(dir.list, "./")

dir <- dir.list[1]

Kmers.all <- NULL

for (dir in dir.list){
    print(dir)
    file.list <- list.files(paste0(dir, "/kmer6"), pattern = ".txt", full.names = T)
    
    nuc <- c("A", "T", "C", "G")
    ## To generalize for A kmersize
    
    res <- data.frame(kmer = apply(expand.grid(nuc, nuc, nuc, nuc, nuc, nuc),
                                   1,
                                   function(x){return(paste0(x, collapse = ""))}),
                      stringsAsFactors = F)
    
    file <- file.list[1]
    for (file in file.list){
        #print(file)
        
        contig <- str_remove(rev(str_split(rev(str_split(file, "/", simplify = T))[1],"_", simplify = T))[1], ".txt")
        dta <- read.table(file, header = F)
        colnames(dta) <- c("kmer", contig)
        res[, contig] <- dta[match(res$kmer, dta$kmer), 2]
    }
    
    
    res[,2:ncol(res)][is.na(res[,2:ncol(res)])] <- 0
    
    sumKmers <- apply(res[, -1], 2, sum, na.rm = T)
    sumKmers
    
    shannonEntropy <- apply(res[, -1], 2, calculate_shannon_entropy)
    shannonEntropy
    
    
    Kmers <- as.data.frame(sumKmers)
    Kmers$contigs <- rownames(Kmers)
    Kmers$shannonEntropy <- shannonEntropy
    Kmers$pop <- dir
    Kmers.all <- rbind(Kmers.all, Kmers)
}

write.table(Kmers.all, "shannonEntropy_all.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


