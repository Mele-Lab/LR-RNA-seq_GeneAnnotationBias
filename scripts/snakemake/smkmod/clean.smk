rule cleandir:
    resources:
        threads = 1
    shell:
        """
        rm -r {params.inputdir}
        touch {output.cleanmock}
        """