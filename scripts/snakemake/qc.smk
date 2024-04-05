rule count_uniq_mapped_sam_reads:
    shell:
        """
        module load samtools
        samtools view -F 4 {input.align} | cut -f1 | sort | uniq -c > {output.txt}
        """

rule count_fq_reads:
    shell:
        """
        echo $(cat {input.fq}|wc -l)/4|bc > {ouput.txt}
        """

rule count_fq_gz_reads:
    shell:
        """
        echo $(gzcat {input.fq}|wc -l)/4|bc > {ouput.txt}
        """

rule count_fa_reads:
    shell:
        """
        echo $(cat {input.fq}|wc -l)/2|bc > {ouput.txt}
        """

rule count_fa_gz_reads:
    shell:
        """
        echo $(gzcat {input.fq}|wc -l)/2|bc > {ouput.txt}
        """
