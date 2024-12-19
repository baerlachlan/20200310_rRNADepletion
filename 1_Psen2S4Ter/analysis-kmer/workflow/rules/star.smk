rule star_align:
    input:
        r1 = rules.bwa_unmapped_to_fq.output.r1,
        r2 = rules.bwa_unmapped_to_fq.output.r2,
        star_index = rules.star_index.output
    output:
        bam = "results/star/bam/{SAMPLE}.bam",
        bam_index = "results/star/bam/{SAMPLE}.bam.bai",
    params:
        base_name = "results/star/bam/{SAMPLE}",
        bam_unsorted = "results/star/bam/{SAMPLE}Aligned.out.bam",
        align_dir = "results/star"
    conda:
        "../envs/star.yml"
    threads: 16
    resources:
        mem_mb = 32000,
        runtime = 240,
    shell:
        """
        STAR \
            --genomeDir {input.star_index}\
            --runThreadN {threads} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand "gunzip -c" \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.base_name}

        samtools sort {params.bam_unsorted} > {output.bam}
        samtools index {output.bam}
        rm {params.bam_unsorted}

        mkdir -p {params.align_dir}/log
        mv {params.base_name}*out {params.align_dir}/log
        mv {params.base_name}*tab {params.align_dir}/log
        """

rule star_mapped_to_fq:
    input:
        bam = rules.star_align.output.bam
    output:
        r1 = "results/star/fastq/{SAMPLE}_R1" + config["fastq_ext"],
        r2 = "results/star/fastq/{SAMPLE}_R2" + config["fastq_ext"],
    conda:
        "../envs/samtools.yml"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 240,
    shell:
        """
        samtools fastq -F 256 --threads {threads} -c 6 -1 {output.r1} -2 {output.r2} {input.bam}
        """