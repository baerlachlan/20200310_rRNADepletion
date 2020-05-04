#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=16

## Modules
module load Salmon/1.1.0-foss-2016b

## Directories
PROJROOT=/fast/users/a1647910/20200310_rRNADepletion
ALIGNDATABWA=${PROJROOT}/2_alignedDataBwa
ALIGNDATASALMON=${PROJROOT}/4_alignedDataSalmon

## Setup for salmon
mkdir -p ${ALIGNDATASALMON}/quant

## Run salmon
for R1 in ${ALIGNDATABWA}/fastq/ERR*.fastq.gz
do

  BNAME=$(basename ${R1%.fastq.gz})

  salmon quant \
    -i ${PROJROOT}/files/salmon_index \
    -l A \
    -r ${R1} \
    --validateMappings \
    --threads ${CORES} \
    -o ${ALIGNDATASALMON}/quant/${BNAME}

done