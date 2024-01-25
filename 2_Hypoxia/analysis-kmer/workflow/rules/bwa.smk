rule bwa_align:
    input:
        rules.bwa_index.output,
        r1 = rules.trim.output.r1,
        ref_rrna = rules.get_refs.output.ref_rrna
    output:
        bam = "results/bwa/bam/{SAMPLE}.sorted.bam",
        bam_index = "results/bwa/bam/{SAMPLE}.sorted.bam.bai",
        log = "results/bwa/log/{SAMPLE}.log",
        rrna_props = "results/bwa/rrna_props/{SAMPLE}.rrna"
    conda:
        "../envs/bwa.yml"
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-02:00:00",
    shell:
        """
        bwa mem -t {resources.cpu} {input.ref_rrna} {input.r1} \
        | samtools sort -o {output.bam} -
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.log}
        awk -F '[(| ]' 'NR == 7 {{print $6}}' {output.log} > {output.rrna_props}
        """

rule bwa_unmapped_to_fq:
    input:
        bam = rules.bwa_align.output.bam
    output:
        r1 = "results/bwa/fastq/{SAMPLE}" + config["fastq_ext"],
    conda:
        "../envs/samtools.yml"
    resources:
        cpu = 8,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-01:00:00",
    shell:
        """
        samtools view -u -h -f 12 -F 256 {input.bam} \
        | samtools sort -n \
        | samtools fastq --threads {resources.cpu} -c 6 -1 {output.r1}
        """