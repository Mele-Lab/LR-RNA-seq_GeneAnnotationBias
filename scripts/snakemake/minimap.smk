rule minimap:
    resources:
        threads = 48
    shell:
        """
        module load minimap2/2.24-r1122
        minimap2 \
            -ax splice \
            --junc-bed {input.junc_bed} \
            -t {resources.threads} \
            -o {output.sam} \
            {input.fa}
        """
