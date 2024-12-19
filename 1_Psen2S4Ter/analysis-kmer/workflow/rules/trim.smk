rule trim:
    input:
        r1 = "results/raw_data/fastq/{SAMPLE}_R1" + config["fastq_ext"],
        r2 = "results/raw_data/fastq/{SAMPLE}_R2" + config["fastq_ext"],
    output:
        r1 = "results/trim_data/fastq/{SAMPLE}_R1" + config["fastq_ext"],
        r2 = "results/trim_data/fastq/{SAMPLE}_R2" + config["fastq_ext"],
        html = "results/trim_data/log/{SAMPLE}.html"
    params:
        qual = config["trim"]["phred_qual"],
        length = config["trim"]["length"],
        extra = config["trim"]["extra"],
    conda:
        "../envs/trim.yml"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120,
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred {params.qual} \
            --length_required {params.length} \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
            {params.extra}
        """