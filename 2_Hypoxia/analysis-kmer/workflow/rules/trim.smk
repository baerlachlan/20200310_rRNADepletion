rule trim:
    input:
        r1 = "results/raw_data/fastq/{SAMPLE}" + config["fastq_ext"],
    output:
        r1 = "results/trim_data/fastq/{SAMPLE}" + config["fastq_ext"],
        html = "results/trim_data/log/{SAMPLE}.html"
    params:
        qual = config["trim"]["phred_qual"],
        length = config["trim"]["length"],
        extra = config["trim"]["extra"],
    conda:
        "../envs/trim.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-02:00:00",
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            --qualified_quality_phred {params.qual} \
            --length_required {params.length} \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
            {params.extra}
        """