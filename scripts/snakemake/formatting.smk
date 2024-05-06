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
        module load python/3.10.2
        python3 extract_UMI.py {input.align} {params.opref} {params.sep}
        """

rule gtf_to_gt_map:
    resources:
        threads = 1,
        mem_gb = 4
    run:
        import pyranges as pr
        df = pr.read_gtf(input.gtf).df
        df = df[['gene_id', 'transcript_id']].drop_duplicates()
        df.to_csv(output.gt_map, header=None, index=False, sep='\t')

rule dedupe_umi:
    resources:
        threads = 4,
        mem_gb = 32
    shell:
        """
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
