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
        samtools view -F 4 {input.align} > {output.align}
        """

rule bam_get_read_ids:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools view {input.bam} | cut -f1 > {output.txt}
        """

rule sam_filt_for_read_ids:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -N {input.read_ids} {input.align} > {output.align}
        """

rule sam_filt_for_primary:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -F 2048 -F 256 {input.align} > {output.align}
        """

rule bam_to_fastq:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools bam2fq {input.unbam} | gzip -c> {output.fastqgz}
        """

# rule sam_to_fq:
#     resources:
#         threads = 8,
#         mem_gb = 32
#     shell:
#         """
#         module load samtools
#         samtools fastq --threads {resources.threads} {input.align} | gzip > {output.fq}
#         """

rule sam_to_bam:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -hSb {input.align} > {output.align}
        """

rule bam_sort:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools sort -O BAM {input.align} > {output.align}
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