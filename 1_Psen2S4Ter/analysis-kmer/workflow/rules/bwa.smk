rule bwa_align:
    input:
        rules.bwa_index.output,
        r1 = rules.trim.output.r1,
        r2 = rules.trim.output.r2,
        ref_rrna = rules.get_refs.output.ref_rrna
    output:
        bam = "results/bwa/bam/{SAMPLE}.sorted.bam",
        bam_index = "results/bwa/bam/{SAMPLE}.sorted.bam.bai",
        log = "results/bwa/log/{SAMPLE}.log",
        rrna_props = "results/bwa/rrna_props/{SAMPLE}.rrna"
    conda:
        "../envs/bwa.yml"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 120,
    shell:
        """
        bwa mem -t {threads} {input.ref_rrna} {input.r1} {input.r2} \
        | samtools sort -o {output.bam} -
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.log}
        awk -F '[(| ]' 'NR == 7 {{print $6}}' {output.log} > {output.rrna_props}
        """

rule bwa_unmapped_to_fq:
    input:
        bam = rules.bwa_align.output.bam
    output:
        r1 = "results/bwa/fastq/{SAMPLE}_R1" + config["fastq_ext"],
        r2 = "results/bwa/fastq/{SAMPLE}_R2" + config["fastq_ext"]
    conda:
        "../envs/samtools.yml"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 60,
    shell:
        """
        samtools view -u -h -f 12 {input.bam} \
        | samtools sort -n \
        | samtools fastq --threads {threads} -c 6 -1 {output.r1} -2 {output.r2}
        """