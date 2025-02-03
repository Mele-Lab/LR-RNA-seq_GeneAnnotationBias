rule count_reads_with_primary_alignment:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view -@ {resources.threads} -F 1060 {input.sam} | wc -l) >> {output.txt}
        """

rule count_unbam_reads:
    resources:
        runtime = 60,
        threads = 8
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view -@ {resources.threads} -c {input.unbam} )  >> {output.txt}
        """

rule count_fq_reads:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(cat {input.fq}|wc -l)/4 )) >> {output.txt}
        """

rule count_fq_gz_reads_from_params:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(zcat {params.fqgz}|wc -l)/4 )) >> {output.txt}
        """

rule count_fq_gz_reads:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(zcat {input.fqgz}|wc -l)/4 )) >> {output.txt}
        """

rule count_fa_reads:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(cat {input.fa}|wc -l)/2 )) >> {output.txt}
        """

rule count_fa_gz_reads:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(( $(gzcat {input.fqgz}|wc -l)/2 )) >> {output.txt}
        """

rule count_sam_mappings:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        module load samtools
        echo "{params.text}" $(samtools view -@ {resources.threads} -c {input.sam}) >> {output.txt}
        """

rule count_text_reads:
    resources:
        runtime = 60,
        threads = 8,
        mem_gb = 32
    shell:
        """
        echo "{params.text}" $(cat {input.txt}|wc -l) >> {output.txt}
        """



# compute read length, query coverage, and alignment identity
# for each read in a given sam / bam file
rule alignment_qc_id:
    resources:
        runtime = 60,
        mem_gb = 32,
        threads = 8
    run:
        df = compute_read_cov_id(input.align, resources.threads)
        df.to_csv(output.tsv, sep='\t', index=False)


# Nanoplot custom script (done by wcouster, nanoplot developer)
rule pseudonanoplot:
    resources:
        runtime = 60,
        threads = 4
    shell:
        """
        module unload miniconda
        module load anaconda
        conda init
        source activate base
        conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot
        python snakemake/fastq_to_tsv.py --fastq {input.fastq} -o {output.tsv}
        """

# Make QC report
rule qcreport:
    resources:
        runtime = 60,
        threads = 48
    shell:
        """
        module load R/4.3.2
        Rscript snakemake/qc_on_nanoplot_data.R {wildcards.sample} {input.nanooutput} {input.readcount} {output.pdf} {output.globalstats} {output.reads} {params.path} {output.reads2}
        """

# Run fastQC

rule fastqc:
    resources:
        runtime = 60,
        threads = 48
    shell:
        """
        module load fastqc
        fastqc \
            -o {params.outdir} \
            -t {resources.threads} \
            -a {params.adapters} \
            -f fastq \
            {input.fastq}
        """
