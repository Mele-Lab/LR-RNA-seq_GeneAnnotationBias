rule sam_filt_unmapped:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -f 4 {input.align} > {output.align}
        """

rule sam_sort:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools sort {input.align} > {output.align}
        """

rule sam_index:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools index {input.align}
        """
