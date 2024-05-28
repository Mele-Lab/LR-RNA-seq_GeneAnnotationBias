rule split_ONT_plus_full_linker:
    resources:
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load anaconda
        conda init
        source activate base
        conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
        /gpfs/projects/bsc83/utils/conda_envs/duplextools_env/bin/duplex_tools split_on_adapter \
            --threads {resources.threads} \
            --adapter_type ONT_sequencing_adapter+CapTrapSeqJoint \
            {params.fastqdir} \
            {params.outdir} \
            PCR
        touch {output.mockout}
        """

rule split_full_linker:
    resources:
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load anaconda
        conda init
        source activate base
        conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
        /gpfs/projects/bsc83/utils/conda_envs/duplextools_env/bin/duplex_tools split_on_adapter \
            --threads {resources.threads} \
            --adapter_type CapTrap_joint \
            --n_bases_to_mask_tail 0 \
            --n_bases_to_mask_head 0 \
            --degenerate_bases 0 \
            --edit_threshold 26 {params.fastqdir} \
            {params.outdir} \
            PCR
        cat {input.mockin2} > {output.mockout}
        """

rule skip_multisplits:
    resources:
        threads = 112
    shell:
        """
        bash snakemake/skipMultiSplits.sh snakemake {params.outdir} {params.splitdir} {wildcards.sample} {resources.threads} {input.mockinput}
        cat {input.mockinput} >> {output.mockout}
        """

