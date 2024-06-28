rule ref_make_gatk_dict:
    resources:
        runtime = 120,
        threads = 112,
        nodes = 1
    shell:
        """
        module load java-openjdk/22.0.1
        module load gatk
        gatk CreateSequenceDictionary -R {input.fa} -O {output.fa_dict}
        """

rule ref_make_fa_ind:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        module load samtools
        samtools faidx {input.fa}
        """

rule call_variants_gatk:
    resources:
        runtime = 120,
        threads = 112,
        nodes = 4
    shell:
        """
        module purge
        module load java-openjdk/17.0.11+9
        module load gatk/4.5.0.0
        gatk \
            --java-options "-Xmx36g -Xms32g" \
            HaplotypeCaller \
            -I {input.bam} \
            -R {input.fa} \
            -O {output.vcf} \
            --sample-ploidy 2 \
            --max-reads-per-alignment-start 0 \
            --dont-use-soft-clipped-bases true \
            --native-pair-hmm-threads {resources.threads} \
            --max-alternate-alleles 6 \
            --max-genotype-count 1024 \
            --create-output-bam-index \
            --sample-name {wildcards.sample}
        """

rule filt_variants:
    params:
        qual_thresh = 1000
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    run:
        df = read_vcf(input.vcf)

        # remove non-ref alleles
        l1 = len(df.index)
        df = df.loc[df.ALT!='<NON_REF>']
        l2 = len(df.index)
        assert l1 != l2

        # remove reads w/ qual under thresh
        df.QUAL = df.QUAL.astype('float')
        df = df.loc[df.QUAL>=params.qual_thresh]
        write_vcf(df, output.vcf, input.vcf)

rule intersect_variants_with_bed:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule rev_intersect_variants_with_bed:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -v -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule merge_variants:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        bcftools merge \
            --use-header {input.header_vcf} \
            {input.vcfs} > {output.vcf}
        """

rule bgzip:
    resources:
        runtime = 120,
        threads = 112,
        nodes = 1
    shell:
        """
        module load htslib
        bgzip -c {input.ifile} -@ {resources.threads} > {output.gz}
        """

rule vcfgz_index:
    resources:
        runtime = 120,
        threads = 112
    shell:
        """
        module load bcftools
        bcftools index --threads {resources.threads} {input.vcfgz} --output {output.index}
        """

rule vcf_index:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        tabix -p vcf {input.vcf}
        """

# remove PL tag because it doesn't work with bcftools
# at higher ploidies
rule vcf_rm_PL:
    resources:
        runtime = 120,
        threads = 1,
        nodes = 1
    shell:
        """
        if [ -s "{input.vcf}" ]; then
            bcftools annotate \
                -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
                {input.vcf} > {output.vcf}
        else
            touch {output.vcf}
        fi
        """

rule vcf_norm:
    resources:
        runtime = 320,
        threads = 112,
        nodes = 2
    shell:
        """
        module load bcftools

        if [ -s "{input.vcf}" ]; then
            bcftools norm \
              -a \
              -m -any \
              --fasta-ref {input.fa} \
              --old-rec-tag INFO \
              --threads {resources.threads}\
              {input.vcf} > {output.vcf}
        else
            touch {output.vcf}
        fi
        """

rule gatk_split_reads:
    resources:
        runtime = 120,
        threads = 112,
        nodes = 10
    shell:
        """
        module load java-openjdk/22.0.1
        module load gatk
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads={resources.threads}" SplitNCigarReads \
          -R {input.fa} \
          -I {input.align} \
          -O {output.align}
        """

rule fix_read_flags:
    resources:
        runtime = 120,
        threads = 112,
        nodes = 10
    shell:
        """
        module load minimap2 samtools \
            gatk R/4.3.2 vcftools bcftools \
            oneapi/2024.1 htslib/1.19.1 tabixpp/1.1.2
        Rscript /gpfs/projects/bsc83/utils/lrRNAseqVariantCalling/tools/flagCorrection.r \
          {input.align} \
          {input.split_align} \
          {output.align} \
          {resources.threads}
        """





rule merge_samples:
    resources:
        runtime = 120,
        threads = 112
    shell:
        """
        module load bcftools

        bcftools merge \
            --missing-to-ref \
            --output {output.vcfgz} \
            --output-type z \
            --threads {resources.threads} \
            {params.inputvcfs}

        """

rule vcf_pruning:
    resources:
        runtime = 120,
        threads = 112
    shell:
        """
        module load plink

        plink2 \
            --vcf {input.vcf} \
            --allow-extra-chr \
            --double-id \
            --psam {input.psam} \
            --set-missing-var-ids @:# \
            --indep-pairwise 50 10 0.1 \
            --threads {resources.threads} \
            --out {params.output_prefix}

        touch {output.mock}
        """

rule vcf_pca:
    resources:
        runtime = 120,
        threads = 112
    shell:
        """
        module load plink

        plink2 \
            --vcf {input.vcf}\
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --extract {params.input_prefix}.prune.in \
            --make-bed \
            --psam {input.psam} \
            --pca \
            --threads {resources.threads} \
            --out {params.output_prefix}

        rm {input.mock}
        touch {output.mock}
        """

rule vcf_pca_no_pruning:
    resources:
        runtime = 120,
        threads = 112
    shell:
        """
        module load plink

        plink2 \
            --vcf {input.vcf}\
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --make-bed \
            --psam {input.psam} \
            --pca \
            --threads {resources.threads} \
            --out {params.output_prefix}

        touch {output.mock}
        """