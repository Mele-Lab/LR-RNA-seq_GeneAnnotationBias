rule add_duplex_status:
    resources:
        mem_gb = 16,
        threads = 8
    shell:
        """
        module load samtools
        samtools view {input.align} | awk '{$1=$1"_"$13; print $0}' > {output.align}
        """

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
