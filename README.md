# Long-read transcriptomics of a diverse human cohort reveals widespread ancestry bias in gene annotation

The data processing and analyses for this manuscript were done by multiple people. We therefore provide this page to outline the locations of code used to perform each task.

## ONT preprocessing
- [**ONT basecalling:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/01_basecalling) DORADO basecalling of raw ONT signal files.
- [**ONT preprocessing:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/02_ONT_preprocessing/) Pipeline to process ONT basecalled data, including duplex read handling, adapter trimming, chimera splitting, and UMI deduplication. Starts from unaligned bam files.
- [**ONT read mapping:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/03_mapping/) Minimap2 mappings to GRCh38 (both annotation-free and annotation-guided).

## Creation of PODER

### [Transcript discovery + merging](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly)

We ran 4 different tools for transcript discovery.

- [**ESPRESSO**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/01_espresso)
- [**IsoQuant**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/02_isoquant)
- [**FLAIR**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/03_flair)
- [**Intron chain (IC) merging:**](https://github.com/fairliereese/240706_pantranscriptome_cerberus_gtf_merge/tree/594b554f0235b5c0d1f40f789e7e0ecacecbbb9c/merge_only_ics) Merge intron chains across samples and tools (ESPRESSO, LyRic, IsoQuant, FLAIR) to generate the unfiltered merged annotation (UMA).
- [**Add novel gene loci:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/novel_gene) Add novel gene IDs to intergenic transcripts found in UMA using [buildLoci](https://github.com/julienlag/buildLoci).

### UMA characterization and filtering

- [**SQANTI:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/04_evaluation/02_sqanti/) Running SQANTI to characterize transcripts in UMA.
- [**Recount3 processing**:](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/04_evaluation/03_recount3/) Processing Recount3 data and intersecting splice junction counts with UMA transcripts.
- [**Unambiguous quantification with FLAIR:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/06_quantification/02_flairquantify/) Quantification of UMA transcripts across the datasets. Used to filter on minimum counts as FLAIR only considers unambiguous transcripts and therefore is more conservative during quantification.
- [**Protein prediction:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/protein) Run protein prediction pipeline on UMA.
- [**Construction of "master table":**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/04_evaluation/04_mastertable/) Compilation of all relevant data (from earlier in this section or otherwise) into single table to facilitate analyses.

## PODER / Enhanced GENCODE creation and characterization
* [**Add gene entries to PODER GTF:**](https://github.com/fairliereese/240903_pt/blob/main/snakemake/novel_annotation_add_gene) Add gene entries to PODER GTF which only had transcript and exon entries, and a gene ID tag.
* [**Enhanced GENCODE creation:**](https://github.com/fairliereese/240903_pt/blob/main/snakemake/merge_v47_poder) Merge novel PODER transcript GTF with GENCODE v47 GTF to create the Enhanced GENCODE annotation.
* [**Protein prediction:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/poder_protein) Run protein prediction pipeline on PODER.
* [**PFAM:**](https://github.com/fairliereese/240903_pt/blob/main/snakemake/pfam) Run PFAM protein domain finder on GENCODE v47 reference annotation and PODER predicted proteins.
* [**Quantification:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/lr-kallisto) Quantification of our LR-RNA-seq dataset with lr-kallisto and PODER.



## Analyses

### ESPRESSO downsampling experiment to show that read depth does not affect our annotation bias finding
- [**ESPRESSO:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/05_downsampling) Running ESPRESSO (separately using RefSeq or GENCODE as the annotation) on downsampled data.
- [**Merge ESPRESSO ICs (GENCODE):**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/merge_espresso) Merging ESPRESSO transcripts from downsampling experiment when using GENCODE v47 as the annotation by intron chain.
- [**Merge ESPRESSO ICs (RefSeq):**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/merge_espresso_refseq) Merging ESPRESSO transcripts from downsampling experiment when using RefSeq v110 as the annotation by intron chain.

### Allele-specific analyses
- [**Allele specific calling:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/08_allele_specifics) ASE and ASTU and subsequent analyses
- [**ASTU example:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/astu_example) Systematically output LR-RNA-seq BAM files split by allele to facilitate visualization of allele-specific transcript usage (ASTU).

### Personalized GRCh38 transcript discovery
* [**Mapping, SQANTI, and variant position intersection with splice sites:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/transcript_discovery_personal): Map our LR-RNA-seq data from the 30 1000G-overlapping samples to their corresponding 2 personalized-GRCh38 haplotypes (ie GRCh38 with the SNPs from each sample incorporated). Run SQANTI on the ESPRESSO results to get splice junctions from each mapping (personalized-GRCh38s and GRCh38 alone). Intersect biallelic SNPs with splice-junction proximal exonic SNPs and splice site intronic SNPs to determine how splice junctions' discovery is affected by genetics.

### Personal assemblies
- [**Download:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/map_personal) Download of personal assemblies from pangenome matching our LR-RNA-seq samples. <!-- Also mapping but we ended up using Fabien's -->

### Other analyses
- [**Intercatalog overlap:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/ics_inter_catalog) Extract unique intron chains from each catalog (CHESS3, ENCODE4, GTEx) and run SQANTI on them.
- [**ENCODE portal LR-RNA-seq parsing:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/encode) Determine genetic ancestry (if annotated) of ENCODE LR-RNA-seq catalog.
- [**GTEx requantification:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/gtex_lr-kallisto) Re-quantify the GTEx LR-RNA-seq dataset from FASTQ using lr-kallisto and PODER.
- [**MAGE requantification:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/mage): Quantification of [MAGE](https://github.com/mccoy-lab/MAGE) RNA-seq dataset using kallisto and multiple annotations (PODER, GENCODE, Enhanced GENCODE).
- [**Novel exon FSTs:**](https://github.com/fairliereese/240903_pt/tree/main/snakemake/pop_div_exon_fsts): Find known exons and novel transcribed regions of novel exons. Compute the Fst values pairwise between each set of populations in the 1000G that overlap ours for the SNP variants that fall within these regions.
