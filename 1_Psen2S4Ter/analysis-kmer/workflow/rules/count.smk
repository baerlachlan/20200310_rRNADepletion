rule count:
    input:
        bams = expand(rules.star_align.output.bam, SAMPLE=samples),
        ref_gtf = rules.get_refs.output.ref_gtf
    output:
        counts = "results/count/counts.out",
        counts_summary = "results/count/counts.out.summary",
    conda:
        "../envs/count.yml"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120,
    shell:
        """
        featureCounts -Q 10 \
            -s 0 \
            -T {threads} \
            -p \
            --fracOverlap 1 \
            -a {input.ref_gtf} \
            -o {output.counts} {input.bams}
        """
