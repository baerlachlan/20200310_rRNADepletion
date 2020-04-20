#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Directories
PROJROOT=/fast/users/a1647910/20200310_rRNADepletion

## Modules
module load Salmon/1.1.0-foss-2016b

## Index
salmon index \
    -t ${PROJROOT}/files/gentrome.fa.gz \
    -d ${PROJROOT}/files/decoys.txt \
    -k 31 \
    -p 16 \
    -i ${PROJROOT}/files/salmon_index