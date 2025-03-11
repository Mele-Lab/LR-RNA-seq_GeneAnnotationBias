# Long-read transcriptomics of a diverse human cohort reveals widespread ancestry bias in gene annotation
Description of directories content
- **01_basecalling:** DORADO basecalling by ONT
- **02_ONT_preprocessing:** Pipeline to process ONT basecalled data starting from unaligned bam files
- **03_mapping:** Minimap mappings to GRCh38 relying or not on an annotation
- **04_transcriptome_assembly:**
  - **01_espresso, 02_isoquant, 03_flair:** transcript discovery pipelines
  - **04_evaluation:**
    - **04_mastertable**: Preparation and explorations of Unfiltered Merged Annotation (UMA) and PODER (Population Diverse Enhanced longReads) annotation. It includes all gene annotation bias analyses.
- **06_quantification**: transcript expression quantification used for transcript filtering
- **08_allele_specifics**: ASE and ASTU and subsequent analyses
- **fd_personal_genomes**: All code for "Personalized genome reference assemblies enhance novel splice junction discovery" section
