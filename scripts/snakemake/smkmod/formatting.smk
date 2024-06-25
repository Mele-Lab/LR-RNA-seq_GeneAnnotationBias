rule add_duplex_status:
    resources:
        runtime = 60,
        mem_gb = 16,
        threads = 8
    run:
        append_duplex_tag_read_name(input.unbam, input.index, output.unbam, params.sep, resources.threads)

rule fq_to_fa:
    resources:
        runtime = 60,
        mem_gb = 4,
        threads = 1
    shell:
        """zcat {input.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}"""

rule extract_umi:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 16
    shell:
        """
        module load hdf5 python/3.12.1
        python snakemake/extract_UMI.py {input.align} {params.opref} {params.sep}
        """

rule fasta_get_read_ids_rm_umi:
    resources:
        runtime = 60,
        threads = 1,
        mem_gb = 16
    shell:
        """
        grep "^>" {input.fa} | cut -c 2- | sed 's/,.*//'> {output.txt}
        """

rule fastqgz_get_read_ids:
    resources:
        runtime = 60,
        threads = 1,
        mem_gb = 16
    shell:
        """
        zcat {input.fq} | awk 'NR % 4 == 1 {{print substr($1, 2)}}' > {output.txt}
        """



# this rule ALSO REMOVES SIMPLEX
rule read_id_union:
    resources:
        runtime = 60,
        threads = 1,
        mem_gb = 4
    shell:
        """
        cat {input.a} {input.b} {input.c} | sort | uniq > {params.temporalfile}
        cat {params.temporalfile}| grep -v ":-1"> {output.txt}
        echo "{params.text}" $(cat {params.temporalfile} | grep ":-1" | wc -l) > {output.count_parent_simplex}
        rm {params.temporalfile}
        """


# Finding Entries in file1.txt but Not in file2.txt:
rule read_id_diff:
    resources:
        runtime = 60,
        threads = 1,
        mem_gb = 4
    shell:
        """
        bash snakemake/diff.sh {input.a} {input.b} {output.txt}

        """

rule gtf_to_gt_map:
    resources:
        runtime = 60,
        threads = 1,
        mem_gb = 16
    shell:
        """
        awk '/\ttranscript\t/{{OFS="\t"}}{{split($0, a, /"/); print a[2], a[4]}}' {input.gtf} |\
        grep ENST | uniq > {output.gt_map}
        """
    # shell:
    #     """
    #     module load anaconda
    #     conda init
    #     source activate base
    #     conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
    #     python snakemake/gtf2gtmap.py {input.gtf} {output.gt_map}
    #     """


rule dedupe_umi:
    resources:
        runtime = 1440,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load miniconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/umi_tools
        umi_tools dedup \
            --extract-umi-method read_id \
            --umi-separator {params.sep} \
            --method adjacency \
            --edit-distance-threshold {params.edit_dist} \
            --per-contig \
            --per-gene \
            --gene-transcript-map {input.gt_map} \
            --stdin {input.align} \
            --stdout {output.align} \
            --log {output.log}
        """


rule fastqgz_filter:
    resources:
        runtime = 60,
        threads = 112,
        mem_gb = 32
    shell:
        """
        echo "test1#####################"
        zcat {input.fq} > $TMPDIR/{wildcards.sample}.fastq
        echo "test2 ##########################"
        grep -A 3 -Ff {input.read_ids} $TMPDIR/{wildcards.sample}.fastq | grep -v "^--$" > {output.cleanfastq}
        echo "DONE ##########################"

        """

rule fastqfilter:
    resources:
        runtime = 60,
        threads=112
    shell:
        """
        module load miniconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/fastqfilter #fastq-filter

        fastq-filter -q {params.phred} {input.align} -o {output.align}
        """

# rule trimadaptors:
#     resources:
#         runtime = 60,
#         threads = 80,
#         queue = "acc_bscls"
#     shell:
#         """
#         module load dorado

#         dorado trim \
#             --threads 80 \
#             --emit-fastq \
#             {input.align} | gzip > {output.align}
#         """

rule porechop:
    resources:
        runtime = 500,
        threads = 112
    shell:
        """
        module load anaconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/porechop

        porechop \
            -i {input.fastq} \
            -o {output.fastqgz} \
            --threads {resources.threads} \
            --min_split_read_size 200 \
            --extra_middle_trim_good_side 2 \
            --extra_middle_trim_bad_side 2
        """

rule gunzip:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        gunzip -c {input.ifile} > {output.ofile}
        """


# Content of adapters.py from porechop
# ADAPTERS = [Adapter('SQK-LSK114_HEAD',
#              start_sequence=('HEAD_ADAPTER', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
#              end_sequence=('HEAD_ADAPTERrv', 'AGCAATACGTAACTGAACGAAGTACATT')),
#             Adapter('SQK-LSK114_TAIL',
#              start_sequence=('TAIL_ADAPTER', 'GCAATACGTAACTGAACGAAGT'),
#              end_sequence=('TAIL_ADAPTERrv', 'ACTTCGTTCAGTTACGTATTGC')),
#             Adapter('FwCapTrapPrimer',
#              start_sequence=('FW_CAP_PRIMER', 'TCGTCGGCAGCGTC'),
#              end_sequence=('FW_CAP_PRIMERrv', 'GACGCTGCCGACGA')),
#             Adapter('RvCapTrapPrimer',
#              start_sequence=('RV_CAP_PRIMER', 'GTCTCGTGGGCTCGG'),
#              end_sequence=('RV_CAP_PRIMERrv', 'CCGAGCCCACGAGAC')),
#             Adapter('FIVE_LINKER_1',
#              start_sequence=('FIVE_LINKER_1', 'AGATGTGTATAAGAGACAG'),
#              end_sequence=('FIVE_LINKER_1rv', 'CTGTCTCTTATACACATCT')),
#             Adapter('FIVE_LINKER_2',
#              start_sequence=('FIVE_LINKER_2', 'GTGGTATCAACGCAGAGTAC'),
#              end_sequence=('FIVE_LINKER_2rv', 'GTACTCTGCGTTGATACCAC')),
#             Adapter('THREE_LINKER',
#              start_sequence=('THREE_LINKER', 'AGATGTGTATAAGAGACAGCGATGCTTTTTTTTTTTTTTTT'),
#              end_sequence=('THREE_LINKERrv', 'AAAAAAAAAAAAAAAAGCATCGCTGTCTCTTATACACATCT')),
#             Adapter('FULL_FIVE_LINKER_1',
#              start_sequence=('FULL_FIVE_LINKER_1', 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'),
#              end_sequence=('FULL_FIVE_LINKER_1rv', 'CTGTCTCTTATACACATCTGACGCTGCCGACGA')),
#             Adapter('FULL_THREE_LINKER',
#              start_sequence=('FULL_THREE_LINKER', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCGATGCTTTTTTTTTTTTTTTT'),
#              end_sequence=('FULL_THREE_LINKERrv', 'AAAAAAAAAAAAAAAAGCATCGCTGTCTCTTATACACATCTCCGAGCCCACGAGAC'))]
