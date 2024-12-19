rule fastqc_raw:
    input:
        "results/raw_data/fastq/{SAMPLE}{PAIRTAG}" + config["fastq_ext"],
    output:
        multiext("results/raw_data/FastQC/{SAMPLE}{PAIRTAG}_fastqc", ".zip", ".html"),
    params:
        outDir = "results/raw_data/FastQC",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        # cpu = 1,
        # ntasks = 1,
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """

rule fastqc_trim:
    input:
        "results/trim_data/fastq/{SAMPLE}{PAIRTAG}" + config["fastq_ext"],
    output:
        multiext("results/trim_data/FastQC/{SAMPLE}{PAIRTAG}_fastqc", ".zip", ".html"),
    params:
        outDir = "results/trim_data/FastQC",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        # cpu = 1,
        # ntasks = 1,
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """

rule fastqc_bwa_bam:
    input:
        "results/bwa/bam/{SAMPLE}.sorted.bam",
    output:
        multiext("results/bwa/FastQC_bam/{SAMPLE}.sorted_fastqc", ".zip", ".html"),
    params:
        outDir = "results/bwa/FastQC_bam",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """

rule fastqc_bwa_fastq:
    input:
        "results/bwa/fastq/{SAMPLE}{PAIRTAG}" + config["fastq_ext"],
    output:
        multiext("results/bwa/FastQC_fastq/{SAMPLE}{PAIRTAG}_fastqc", ".zip", ".html"),
    params:
        outDir = "results/bwa/FastQC_fastq",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """

rule fastqc_star_bam:
    input:
        "results/star/bam/{SAMPLE}.bam",
    output:
        multiext("results/star/FastQC_bam/{SAMPLE}_fastqc", ".zip", ".html"),
    params:
        outDir = "results/star/FastQC_bam",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """

rule fastqc_star_fastq:
    input:
        "results/star/fastq/{SAMPLE}{PAIRTAG}" + config["fastq_ext"],
    output:
        multiext("results/star/FastQC_fastq/{SAMPLE}{PAIRTAG}_fastqc", ".zip", ".html"),
    params:
        outDir = "results/star/FastQC_fastq",
    conda:
        "../envs/fastqc.yml"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 60,
    shell:
        """
        HOME=/hpcfs/users/$USER
        fastqc -t {threads} -o {params.outDir} --noextract {input}
        """