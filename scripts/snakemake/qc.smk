rule count_reads_with_primary_alignment:
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view -F 1060 {input.sam}) >> {output.txt}
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
    shell:
        """
        echo "{params.text}" $(( $(cat {input.fq}|wc -l)/4 )) >> {ouput.txt}
        """

rule count_fq_gz_reads:
    shell:
        """
        echo "{params.text}" $(( $(gzcat {input.fqgz}|wc -l)/4 )) >> {ouput.txt}
        """

rule count_fa_reads:
    shell:
        """
        echo "{params.text}" $(( $(cat {input.fa}|wc -l)/2 )) >> {ouput.txt}
        """

rule count_fa_gz_reads:
    shell:
        """
        echo "{params.text}" $(( $(gzcat {input.fagz}|wc -l)/2 )) >> {ouput.txt}
        """

rule count_sam_mappings:
    shell:
        """
        echo "{params.text}" $(cat {input.sam}|wc -l) >> {ouput.txt}
        """

rule count_text_reads:
    shell:
        """
        echo "{params.text}" $(cat {input.txt}|wc -l) >> {ouput.txt}
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