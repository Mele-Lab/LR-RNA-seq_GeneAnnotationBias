rule split_ONT_plus_full_linker:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load anaconda
        conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
        duplex_tools split_on_adapter --threads {resources.threads} --adapter_type ONT_sequencing_adapter+CapTrapSeqJoint {input.fastq} {output.fastq.gz} PCR
        
        """

rule split_full_linker:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load anaconda
        conda activate .
        duplex_tools split_on_adapter --threads {resources.threads} --adapter_type CapTrap_joint --n_bases_to_mask_tail 0 --n_bases_to_mask_head 0 --degenerate_bases 0 --edit_threshold 26 {input.fastq.gz} {output.fastq.gz} PCR
        
        """