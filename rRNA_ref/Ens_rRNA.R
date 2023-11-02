## This script gets rRNA sequences using biomaRt
library(biomaRt)
library(tidyverse)
library(AnnotationHub)
library(ensembldb)

## Check biotypes to use as values in getBM()
ah <- AnnotationHub() %>%
    subset(species == "Danio rerio") %>%
    subset(rdataclass == "EnsDb")
ensDb <- ah[["AH83189"]]
genes <- genes(ensDb)
genes$gene_biotype %>%
    unique() %>%
    .[str_detect(., "rRNA")]

## Get and export rRNA sequences as FASTA
mart <- useEnsembl(
    biomart = "genes",
    dataset = "drerio_gene_ensembl",
    version = 101
)
rRna <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype"),
    filters = "biotype",
    values = c("rRNA", "Mt_rRNA"),
    mart = mart
)
seq <- getSequence(
    id = rRna$ensembl_gene_id,
    type = "ensembl_gene_id",
    seqType = "cdna",
    mart = mart
)
exportFASTA(
    seq,
    "/hpcfs/users/a1647910/20200310_rRNADepletion/rRNA_ref/Ens_rRNA.fa"
)

