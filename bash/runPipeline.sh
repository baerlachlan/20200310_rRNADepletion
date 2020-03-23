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

## Cores
CORES=16

## Modules
module load FastQC/0.11.7
module load AdapterRemoval/2.2.1-foss-2016b
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.9-foss-2016b

## Ribosomal RNA Reference files
REFS=/data/biorefs/rRNA/danio_rerio/bwa/danRer11

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

## Setup for rRNA alignment
ALIGNDATA=${PROJROOT}/2_alignedData
mkdir -p ${ALIGNDATA}/bam
mkdir -p ${ALIGNDATA}/log
mkdir -p ${ALIGNDATA}/FastQC


##--------------------------------------------------------------------------------------------##
## FastQC on the raw data
##--------------------------------------------------------------------------------------------##

# fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/*.fastq.gz

##--------------------------------------------------------------------------------------------##
## Trimming the merged data
##--------------------------------------------------------------------------------------------##

# for R1 in ${RAWDATA}/fastq/*.fastq.gz
# do

#  echo -e "Currently working on ${R1}"

#  # Now create the output filenames
#  out1=${TRIMDATA}/fastq/$(basename $R1)
#  BNAME=${TRIMDATA}/fastq/$(basename ${R1%.fastq.gz})
#  echo -e "Output file will be ${out1}"
#  echo -e "Trimming:\t${BNAME}"

#  # Trim
#  AdapterRemoval \
#    --gzip \
#    --trimns \
#    --trimqualities \
#    --minquality 30 \
#    --threads ${CORES} \
#    --basename ${BNAME} \
#    --output1 ${out1} \
#    --file1 ${R1}

# done

# ## Move the log files into their own folder
# mv ${TRIMDATA}/fastq/*settings ${TRIMDATA}/log

# ## Run FastQC
# fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastq/*.fastq.gz

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to ribosomal RNA
##--------------------------------------------------------------------------------------------##

# ## Aligning and sorting
# for R1 in ${TRIMDATA}/fastq/*fastq.gz
# do

#   out=${ALIGNDATA}/bam/$(basename ${R1%.fastq.gz}).sorted.bam
#   echo -e "Output file will be ${out}"

#   bwa mem -t ${CORES} ${REFS} ${R1} \
#   | samtools view -uh \
#   | samtools sort -o ${out}

# done

## Fastqc, indexing and flagstats
for BAM in ${ALIGNDATA}/bam/SRR218*.bam
do

  out=${ALIGNDATA}/log/$(basename ${BAM%.sorted.bam})

  fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATA}/FastQC --noextract ${BAM}
  samtools index ${BAM}
  samtools stats ${BAM} > ${out}.stats
  samtools flagstat ${BAM} > ${out}.flagstat
  samtools idxstats ${BAM} > ${out}.idxstats

done