rule add_duplex_status:
    resources:
        mem_gb = 16,
        threads = 8
    run:
        append_duplex_tag_read_name(input.align, output.align, resources.threads)

rule fq_to_fa:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        """zcat {input.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}"""

rule extract_umi:
    resources:
        threads = 1
    shell:
        """
        module load python/3.10.2
        python3 extract_UMI.py {input.align} {params.opref}
        """
