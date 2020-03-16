#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH -o /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=16

########################################
## Full run of the data from Novogene ##
########################################

## Modules
module load FastQC/0.11.7

# ## Genomic Data Files
# REFS=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/gtf_temp/
# GTF=${REFS}Danio_rerio.GRCz11.94.chr.gtf

## Directories
PROJROOT=/fast/users/a1647910/20200310_rRNADepletion

## Directories for Initial FastQC
RAWDATA=${PROJROOT}/0_rawData
mkdir -p ${RAWDATA}/FastQC

## Setup for Trimmed data
TRIMDATA=${PROJROOT}/1_trimmedData
mkdir -p ${TRIMDATA}/fastq
mkdir -p ${TRIMDATA}/FastQC
mkdir -p ${TRIMDATA}/log

## Setup for genome alignment
ALIGNDATA=${PROJROOT}/2_alignedData
mkdir -p ${ALIGNDATA}/log
mkdir -p ${ALIGNDATA}/bam
mkdir -p ${ALIGNDATA}/FastQC
mkdir -p ${ALIGNDATA}/featureCounts


##--------------------------------------------------------------------------------------------##
## FastQC on the raw data
##--------------------------------------------------------------------------------------------##

fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/*.fastq.gz
