fa_gz = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        "dna.primary_assembly.fa.gz"
    ])
)
fa_gz_path = os.path.join("resources", fa_gz)
fa_gz_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "fasta",
    config["ref"]["species"],
    "dna",
    fa_gz
)
fa = fa_gz.rstrip(".gz")
fa_path = os.path.join("resources", fa)
gtf = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        str(config["ref"]["ensembl_release"]),
        "chr.gtf.gz"
    ])
)
gtf_path = os.path.join("resources", gtf)
gtf_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "gtf",
    config["ref"]["species"],
    gtf
)
rrna = config["ref"]["rrna"]
rrna_path = os.path.join("resources", os.path.basename(rrna))

rule get_refs:
    input:
        ref_rrna = rrna,
    output:
        ref_fa = fa_path,
        ref_gtf = gtf_path,
        ref_rrna = rrna_path
    params:
        fa_gz_url = fa_gz_url,
        fa_gz_path = fa_gz_path,
        gtf_url = gtf_url,
    threads: 1
    shell:
        """
        rsync -avP {params.fa_gz_url} {params.fa_gz_path}
        gzip -d {params.fa_gz_path}
        rsync -avP {params.gtf_url} {output.ref_gtf}
        ln -s {input.ref_rrna} {output.ref_rrna}
        """

rule star_index:
    input:
        ref_fa = rules.get_refs.output.ref_fa,
        ref_gtf = rules.get_refs.output.ref_gtf,
    output:
        directory("resources/star"),
    params:
        overhang = config["align"]["read_length"] - 1,
    conda:
        "../envs/star.yml"
    threads: 16
    resources:
        mem_mb = 32000,
        runtime = 90,
    shell:
        """
        zcat {input.ref_gtf} > temp.gtf

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_fa} \
            --sjdbGTFfile temp.gtf \
            --sjdbOverhang {params.overhang}

        rm temp.gtf
        """

rule bwa_index:
    input:
        rules.get_refs.output.ref_rrna,
    output:
        multiext(rules.get_refs.output.ref_rrna, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda:
        "../envs/bwa.yml"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 10,
    shell:
        """
        bwa index {input}
        """