rule add_duplex_status:
    resources:
        mem_gb = 16,
        threads = 8
    run:
        append_duplex_tag_read_name(input.align, output.align, resources.threads, params.sep)

rule fq_to_fa:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        """zcat {input.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}"""

rule extract_umi:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        module load hdf5 python/3.12.1
        python snakemake/extract_UMI.py {input.align} {params.opref} {params.sep}
        """

rule fasta_get_read_ids:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        grep "^>" {input.fa} | cut -c 2- > {output.txt}
        """

rule read_id_union:
    resources:
        threads = 1,
        mem_gb = 4
    run:
        a = set(pd.read_csv(input.a, header=None)[0].tolist())
        b = set(pd.read_csv(input.b, header=None)[0].tolist())
        union = list(a|b)
        df = pd.DataFrame()
        df['read_id'] = union
        df.to_csv(output.txt, header=False)

rule read_id_diff:
    resources:
        threads = 1,
        mem_gb = 4
    run:
        a = set(pd.read_csv(input.a, header=None)[0].tolist())
        b = set(pd.read_csv(input.b, header=None)[0].tolist())
        diff = list(a-b)
        df = pd.DataFrame()
        df['read_id'] = diff
        df.to_csv(output.txt, header=False)

rule gtf_to_gt_map:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        awk '/\ttranscript\t/{{OFS="\t"}}{{split($0, a, /"/); print a[2], a[4]}}' {input.gtf} |\
         grep ENST | uniq> {output.gt_map}
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
        threads = 4,
        mem_gb = 32
    shell:
        """
        module load miniconda
        source activate sqanti3-snakemake
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
