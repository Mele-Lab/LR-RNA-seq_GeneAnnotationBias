        #!/bin/bash

        module load anaconda
        conda init
        source activate base
        conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot
        python snakemake/fastq_to_tsv.py --fastq data/fourkbam/fourkbam.fastq -o data/fourkbam/qc/test3