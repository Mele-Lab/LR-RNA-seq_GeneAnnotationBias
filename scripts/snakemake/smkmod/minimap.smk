rule minimap:
    resources:
        runtime = 300,
        threads = 112
    shell:
        """
        module load minimap2
        minimap2 \
            -ax splice \
            --junc-bed {input.junc_bed} \
            --MD \
            -t {resources.threads} \
            -o {output.sam} \
            {input.ref_fa} \
            {input.fa}
        """
