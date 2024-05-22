rule count_reads_with_primary_alignment:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view -F 1060 {input.sam} | wc -l) >> {output.txt}
        """

rule count_unbam_reads:
    resources:
        threads = 8
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view {input.unbam} | wc -l)  >> {output.txt}
        """

rule count_fq_reads:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(cat {params.fq}|wc -l)/4 )) >> {output.txt}
        """

rule count_fq_gz_reads:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(zcat {params.fqgz}|wc -l)/4 )) >> {output.txt}
        """

rule count_fa_reads:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(cat {input.fa}|wc -l)/2 )) >> {output.txt}
        """

rule count_fa_gz_reads:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(gzcat {input.fqgz}|wc -l)/2 )) >> {output.txt}
        """

rule count_sam_mappings:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(cat {input.sam}|wc -l) >> {output.txt}
        """

rule count_text_reads:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(cat {input.txt}|wc -l) >> {output.txt}
        """

rule nanoplot:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load anaconda
        conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot
        /gpfs/projects/bsc83/utils/conda_envs/nanoplot/bin/NanoPlot \
            -t {resources.threads} \
            -o {params.outdir} \
            -p {params.prefix} \
            --tsv_stats \
            -f png \
            --{params.filetype} {input.reads} #ubam or fastq
        """