rule add_duplex_status:
    resources:
        mem_gb = 16,
        threads = 8
    run:
        append_duplex_tag_read_name(input.unbam, input.index, output.unbam, params.sep, resources.threads)

rule fq_to_fa:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        """zcat {input.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}"""

rule extract_umi:
    resources:
        threads = 12,
        mem_gb = 16
    shell:
        """
        module load hdf5 python/3.12.1
        python snakemake/extract_UMI.py {input.align} {params.opref} {params.sep}
        """

rule fasta_get_read_ids_rm_umi:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        grep "^>" {input.fa} | cut -c 2- | sed 's/,.*//'> {output.txt}
        """

rule fastqgz_get_read_ids:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        zcat {input.fq} | awk 'NR % 4 == 1 {{print substr($1, 2)}}' > {output.txt}
        """



# this rule ALSO REMOVES SIMPLEX
rule read_id_union:
    resources:
        threads = 1,
        mem_gb = 4
    shell:
        """
        cat {input.a} {input.b} {input.c} | sort | uniq > {params.temporalfile}
        cat {params.temporalfile}| grep -v ":-1"> {output.txt}
        echo "{params.text}" $(cat {params.temporalfile} | grep ":-1" | wc -l) > {output.count_parent_simplex}
        rm {params.temporalfile}
        """

# rule read_id_diff:
#     resources:
#         threads = 1,
#         mem_gb = 4
#     run:
#         a = set(pd.read_csv(input.a, header=None)[0].tolist())
#         b = set(pd.read_csv(input.b, header=None)[0].tolist())
#         diff = list(a-b)
#         df = pd.DataFrame()
#         df['read_id'] = diff
#         df.to_csv(output.txt, header=False)

# Finding Entries in file1.txt but Not in file2.txt:
rule read_id_diff:
    resources:
        threads = 1,
        mem_gb = 4
    shell:
        """
        bash snakemake/diff.sh {input.a} {input.b} {output.txt}

        """

rule gtf_to_gt_map:
    resources:
        threads = 1,
        mem_gb = 16
    shell:
        """
        awk '/\ttranscript\t/{{OFS="\t"}}{{split($0, a, /"/); print a[2], a[4]}}' {input.gtf} |\
        grep ENST | uniq > {output.gt_map}
        """
    # shell:
    #     """
    #     module load anaconda
    #     conda init
    #     source activate base
    #     conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
    #     python snakemake/gtf2gtmap.py {input.gtf} {output.gt_map}
    #     """


rule dedupe_umi:
    resources:
        threads = 4,
        mem_gb = 32
    shell:
        """
        module load miniconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/umi_tools
        umi_tools dedup \
            --extract-umi-method read_id \
            --umi-separator {params.sep} \
            --method adjacency \
            --edit-distance-threshold {params.edit_dist} \
            --per-contig \
            --per-gene \
            --gene-transcript-map {input.gt_map} \
            --stdin {input.align} \
            --stdout {output.align} \
            --log {output.log} \
            --output-stats {params.stats}
        """

rule assess_dedupe_umi:
    resources:
        threads = 4,
        mem_gb = 32
    shell:
        """ 
        module load miniconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/umi_tools
        umi_tools group \
            --method adjacency \
            --edit-distance-threshold=2 \
            --umi-separator {params.sep} \
            --per-gene \
            --per-contig \
            --gene-transcript-map {input.gt_map} \
            -I {input.align} \
            --group-out {output.group} \
            --log {output.log}
        """



rule fastqgz_filter:
    resources:
        threads = 8,
        mem_gb = 32
    shell:
        """
        mkdir data/temp
        echo "test1#####################"
        zcat {input.fq} > "data/temp/{wildcards.sample}.fastq"
        echo "test2 ##########################"
        grep -A 3 -Ff {input.read_ids} "data/temp/{wildcards.sample}.fastq" | grep -v "^--$" > {output.cleanfastq}
        echo "test3 ##########################"
        gzip -c {output.cleanfastq} > {output.cleanfastqgz}
        echo "test4 ##########################"
        rm "data/temp/{wildcards.sample}.fastq"
        rm -r data/temp
        """

rule fastqfilter:
    resources:
        threads=8
    shell:
        """
        module load miniconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/fastqfilter #fastq-filter

        fastq-filter -q {params.phred} {input.align} -o {output.align}
        """

rule trimadaptors:
    resources:
        threads = 80,
        queue = "acc_bscls"
    shell:
        """
        module load dorado

        dorado trim \
            --threads 80 \
            --emit-fastq \
            {input.align} | gzip > {output.align}
        """

rule porechop:
    resources:
        threads = 112
    shell:
        """
        module load anaconda
        source activate base
        conda init
        conda activate /gpfs/projects/bsc83/utils/conda_envs/porechop

        porechop \
            -i {input.fastqgz} \
            -o {output.fastqgz} \
            --threads {resources.threads} \
            --min_split_read_size 200 \
            --extra_middle_trim_good_side 2 \
            --extra_middle_trim_bad_side 2
        """