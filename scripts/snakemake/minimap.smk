rule minimap:
    resources:
        threads = 48
    shell:
        """
        minimap2 \
            -ax splice \
            --junc-bed {input.junc_bed} \
            -t {resources.threads} \
            -o {output.sam} \
            {output.mmi} \
            {input.fa}
        """
