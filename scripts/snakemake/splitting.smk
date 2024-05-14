rule split_ONT_plus_full_linker:
    resources:
        threads = 8,
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
        threads = 8,
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
        threads = 8
    shell:
        """
        bash snakemake/skipMultiSplits.sh . {params.outdir} {params.splitdir} {params.sample} {resources.threads} {input.mockinput}
        cat {input.mockinput} >> {output.mockout}
        """

# ex [fomr pclavell]

# mkdir splitdir
# mkdir splitdir/split_1step/
# mkdir splitdir/split_2step/
# duplex_tools split_on_adapter --threads 8 --adapter_type "ONT_sequencing_adapter+CapTrapSeqJoint" 01_basecalling/data/modifications/tests/ splitdir/split_1step/20240212_HS_16_CH4_GM18772_40kreads_split PCR
# duplex_tools split_on_adapter --threads 8 --n_bases_to_mask_tail 0 --n_bases_to_mask_head 0 --degenerate_bases 0 --edit_threshold 26 --adapter_type "CapTrap_joint" splitdir/split_1step/20240212_HS_16_CH4_GM18772_40kreads_split splitdir/split_2step/20240212_HS_16_CH4_GM18772_40kreads_split PCR
# mkdir newdir
# bash ONT_preprocessing/scripts/snakemake/skipMultiSplits.sh ONT_preprocessing/scripts/snakemake newdir splitdir/split 20240212_HS_16_CH4_GM18772_40kreads_split 8
