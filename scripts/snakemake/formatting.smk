rule fq_to_fa:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        """sed -n '1~4s/^@/>/p;2~4p' {input.fq} > {output.fa}"""
