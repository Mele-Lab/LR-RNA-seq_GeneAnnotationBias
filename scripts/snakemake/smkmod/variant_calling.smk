rule ref_make_gatk_dict:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load gatk
        gatk-launch CreateSequenceDictionary -R {input.fa}
        """

rule ref_make_fa_ind:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load samtools
        samtools faidx {input.fa}
        """

rule call_variants_gatk:
    resources:
        threads = 16,
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
            -ERC BP_RESOLUTION \
            --sample-ploidy 2 \
            --max-reads-per-alignment-start 0 \
            --dont-use-soft-clipped-bases true \
            --native-pair-hmm-threads {resources.threads} \
            --max-alternate-alleles 6 \
            --max-genotype-count 1024 \
            --bamout {output.bam} \
            --create-output-bam-index \
            --sample-name {wildcards.sample}
        """

rule filt_variants:
    params:
        qual_thresh = 1000
    resources:
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
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule rev_intersect_variants_with_bed:
    resources:
        threads = 1,
        nodes = 1
    shell:
        """
        module load bedtools
        bedtools intersect -v -a {input.vcf} -b {input.bed} > {output.bed}
        """

rule merge_variants:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        bcftools merge \
            --use-header {input.header_vcf} \
            {input.vcfs} > {output.vcf}
        """

rule bgzip:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'htslib'
    shell:
        """
        bgzip -c {input.ifile} > {output.gz}
        """

rule vcf_index:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    shell:
        """
        tabix -p vcf {input.vcf}
        """

# remove PL tag because it doesn't work with bcftools
# at higher ploidies
rule vcf_rm_PL:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
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
        threads = 1,
        nodes = 2
    conda:
        'base'
    shell:
        """
        if [ -s "{input.vcf}" ]; then
            bcftools norm \
              -a \
              -m -any \
              --fasta-ref {input.fa} \
              --old-rec-tag INFO \
              {input.vcf} > {output.vcf}
        else
            touch {output.vcf}
        fi
        """
