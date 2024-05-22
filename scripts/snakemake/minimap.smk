rule minimap:
    resources:
        threads = 112
    shell:
        """
        module load minimap2
        minimap2 \
            -ax splice \
            --junc-bed {input.junc_bed} \
            -t {resources.threads} \
            -o {output.sam} \
            {input.ref_fa} \
            {input.fa}
        """

