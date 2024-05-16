rule make_blast_db:
    resources:
        threads = 8,
        mem_gb = 16
    shell:
        """
        module load blast
        makeblastdb -in {input.fa} \
            -input_type fasta \
            -dbtype nucl \
            -out {params.db}
        """

# TODO figure out how to iterate over parameters for optimazation
# probably make standalone
rule blast:
    resources:
        threads = 25,
        mem_gb = 64
    shell:
        """
        module load blast/2.11.0
        module load samtools
        blastn -query {input.fa} \
            -db {params.db} \
            -outfmt 5 \
            -out {output.xml} \
            -num_threads 48 \
            -strand 'both' \
            -task 'blastn' \
            -gapopen 3 \
            -gapextend 1 \
            -penalty -2 \
            -reward 1
        """

rule blast2bam:
    resources:
        threads = 8,
        mem_gb = 16
    shell:
        """
        {params.blast2bam} \
            --minAlignLength 44 \
            {input.xml} \
            {input.ref_fa} \
            {input.fa} > {output.sam}
        """
