rule split_ONT_plus_full_linker:
    resources:
        threads = 112,
        mem_gb = 32,
        runtime = 300
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
        mem_gb = 32,
        runtime = 300
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
        threads = 112,
        runtime = 300
    shell:
        """
        bash snakemake/skipMultiSplits.sh snakemake {params.outdir} {params.splitdir} {wildcards.sample} {resources.threads} {input.mockinput}
        """

# Substitute part of the content of split_on_adapter.py from duplex-tools with this
# CUSTOMIZE ADAPTERS
# HEAD_ADAPTER = 'AATGTACTTCGTTCAGTTACGTATTGCT'
# TAIL_ADAPTER = 'GCAATACGTAACTGAACGAAGT'
# FW_CAP_PRIMER = 'TCGTCGGCAGCGTC'
# RV_CAP_PRIMER = 'GTCTCGTGGGCTCGG'
# FIVE_LINKER = 'AGATGTGTATAAGAGACAGNNNNNNNNNNNNNNNNGTGGTATCAACGCAGAGTAC'
# THREE_LINKER = 'AGATGTGTATAAGAGACAGCGATGCTTTTTTTTTTTTTTTT'
# FULL_FIVE_LINKER = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNNNNNNNGTGGTATCAACGCAGAGTAC'
# FULL_THREE_LINKER = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCGATGCTTTTTTTTTTTTTTTT'

# rctrans = str.maketrans('ACGT', 'TGCA')


# def rev_comp(seq):
#     """Reverse complement a DNA sequence."""
#     return str.translate(seq, rctrans)[::-1]


# def build_targets(
#         n_bases_to_mask_head, n_bases_to_mask_tail,adapter_type, degenerate_bases=None,
#         pcr_primers=(
#             'ACTTGCCTGTCGCTCTATCTTCGGCGTCTGCTTGGGTGTTTAACC',  #default
#             'TTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACCT'),
#         head_adapter=HEAD_ADAPTER, tail_adapter=TAIL_ADAPTER,
#         n_replacement=None):
#     print(adapter_type)

#     if adapter_type=="ONT_sequencing_adapter":
#        head_adapter=(
#             HEAD_ADAPTER)
#        tail_adapter=(
#             TAIL_ADAPTER) 
#        pcr_primers=(
#             '',
#             '')
#     elif adapter_type=="ONT_sequencing_adapter+CapTrapSeqJoint":
#         head_adapter=(
#             HEAD_ADAPTER)
#         tail_adapter=(
#             TAIL_ADAPTER)
#         pcr_primers=(
#             FULL_FIVE_LINKER,
#             FULL_THREE_LINKER)
#     elif adapter_type=="CapTrap_primer":
#        head_adapter=(
#             '')
#        tail_adapter=(
#             '')
#        pcr_primers=(
#             FW_CAP_PRIMER,
#             RV_CAP_PRIMER)
#     elif adapter_type=="CapTrap_joint":
#        head_adapter=(
#             '')
#        tail_adapter=(
#             '')
#        pcr_primers=(
#             FULL_FIVE_LINKER,
#             FULL_THREE_LINKER)
#     else:
#        head_adapter=(
#             '')
#        tail_adapter=(
#             '')
#        pcr_primers=(
#             FIVE_LINKER,
#             THREE_LINKER)