rule jellyfish_k5:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k5/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k5/{SAMPLE}.jf"
    params:
        k_size = 5
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k6:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k6/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k6/{SAMPLE}.jf"
    params:
        k_size = 6
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k7:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k7/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k7/{SAMPLE}.jf"
    params:
        k_size = 7
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k8:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k8/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k8/{SAMPLE}.jf"
    params:
        k_size = 8
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k9:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k9/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k9/{SAMPLE}.jf"
    params:
        k_size = 9
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k10:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k10/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k10/{SAMPLE}.jf"
    params:
        k_size = 10
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k11:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k11/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k11/{SAMPLE}.jf"
    params:
        k_size = 11
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """

rule jellyfish_k12:
    input:
        r1 = rules.star_mapped_to_fq.output.r1,
        r2 = rules.star_mapped_to_fq.output.r2
    output:
        dumps = "results/jellyfish/k12/{SAMPLE}.dumps",
        k_counts = "results/jellyfish/k12/{SAMPLE}.jf"
    params:
        k_size = 12
    conda:
        "../envs/jellyfish.yml",
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-04:00:00"
    shell:
        """
        jellyfish count -m {params.k_size} -s 8G --bf-size 8G -C -t {resources.cpu} -o {output.k_counts} <(zcat {input.r1}) <(zcat {input.r2})
        jellyfish dump -c {output.k_counts} > {output.dumps}
        """