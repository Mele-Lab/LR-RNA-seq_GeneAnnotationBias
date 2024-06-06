rule bamindex:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools index -b --threads {resources.threads} {input.bam} -o {output.index}
        """

rule sam_filt_unmapped:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -h -F 4 -@ {resources.threads} {input.align} > {output.align}
        """

rule bam_get_read_ids:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools view -@ {resources.threads} {input.bam} | cut -f1 | sort | uniq | sed 's/,.*//' > {output.txt}
        """

rule bam_get_unmapped_read_ids:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools view -f 4 -@ {resources.threads} {input.bam} | cut -f1 | sort | uniq | sed 's/,.*//'> {output.txt}
        """

rule sam_filt_for_read_ids:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -@ {resources.threads} -N {input.read_ids} {input.align} > {output.align}
        """

rule sam_filt_for_primary:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -h -F 2048 -F 256 -@ {resources.threads} {input.align} > {output.align}
        """

rule bam_to_fastq:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools bam2fq -@ {resources.threads} {input.unbam} | gzip -c> {output.fastqgz}
        """

# rule sam_to_fq:
#     resources:
        runtime = 120,
#         threads = 112,
#         mem_gb = 32
#     shell:
#         """
#         module load samtools
#         samtools fastq --threads {resources.threads} {input.align} | gzip > {output.fq}
#         """

rule sam_to_bam:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -@ {resources.threads} -hSb {input.align} > {output.align}
        """

rule bam_sort:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools sort -@ {resources.threads} -O BAM {input.align} > {output.align}
        """

rule bam_index:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools index -@ {resources.threads} {input.align}
        """

rule align_filt_unmapped_supp:
    resources:
        runtime = 120,
        threads = 112,
        mem_gb = 32
    shell:
        """
        module load samtools
        samtools view -@ {resources.threads} -h -F 4 -F 2048 {input.align} > {output.align}
        """