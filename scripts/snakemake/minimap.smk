rule minimap:
    resources:
        threads = 25
    shell:
        """
        module load minimap2/2.24-r1122
        minimap2 \
            -ax splice \
            --junc-bed {input.junc_bed} \
            -t {resources.threads} \
            -o {output.sam} \
            {input.ref_fa} \
            {input.fa}
        """

rule filt_minimap_reads:
    resources:
        threads = 2
    shell:
        """
        module load samtools
        """
