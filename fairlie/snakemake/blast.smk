rule make_blast_db:
    shell:
        """
        module load blast
        makeblastdb -in {input.fa} \
            -input_type fasta \
            -dbtype nucl \
            -out {output.db}
        """

rule blast:
    resources:
        threads = 48,
        mem_gb = 64
    shell:
        """
        module load blast/2.11.0
        module load samtools
        blastn -query {input.fa} \
            -db {input.db} \
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
