rule cleandir:
    resources:
        runtime = 60,
        threads = 1
    shell:
        """
        rm -r {params.inputdir}
        touch {output.cleanmock}
        """