rule bamindex:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools index -b --threads {resources.threads} {input.bam} -o {output.index}
        """

rule sam_filt_unmapped:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -f 4 {input.align} --threads {resources.threads} > {output.align}
        """

rule sam_to_fq:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools fastq --threads {resources.threads} {input.align} | gzip > {output.fq}
        """

rule sam_to_bam:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -hSb {input.align} > {output.align}
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

rule bam_index:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools index {input.align}
        """

rule align_filt_unmapped_supp:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -h -F 4 -F 2048 {input.align} > {output.align}
        """