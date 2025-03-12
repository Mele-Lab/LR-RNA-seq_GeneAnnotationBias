# Long-read transcriptomics of a diverse human cohort reveals widespread ancestry bias in gene annotation
Description of directories content
- [**01_basecalling:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/01_basecalling) DORADO basecalling by ONT
- [**02_ONT_preprocessing:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/02_ONT_preprocessing/) Pipeline to process ONT basecalled data starting from unaligned bam files
- [**03_mapping:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/03_mapping/) Minimap mappings to GRCh38 relying or not on an annotation
- [**04_transcriptome_assembly:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly)
  - [**01_espresso**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/01_espresso), [**02_isoquant**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/02_isoquant), [**03_flair:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/03_flair) transcript discovery pipelines
  - [**04_evaluation:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/04_evaluation)
    - [**04_mastertable:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/04_evaluation/04_mastertable) Preparation and explorations of Unfiltered Merged Annotation (UMA) and PODER (Population Diverse Enhanced longReads) annotation. It includes all gene annotation bias analyses.
  - [**05_downsampling:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/04_transcriptome_assembly/05_downsampling) Downsampling experiments using ESPRESSO with uniform read depths.
- [**06_quantification:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/06_quantification) transcript expression quantification used for transcript filtering
- [**08_allele_specifics:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/08_allele_specifics) ASE and ASTU and subsequent analyses
- [**fd_personal_genomes:**](https://github.com/Mele-Lab/LR-RNA-seq_GeneAnnotationBias/tree/master/fd_personal_genome) All code for "Personalized genome reference assemblies enhance novel splice junction discovery" section
- [**fairlie:**](https://github.com/fairliereese/240903_pt/tree/6fc2f746ed78187ecd6dae2398c73c90080e6dcd) All data processing and analysis code ran by Fairlie (see folder for more details on specific analyses)
