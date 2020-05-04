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
#SBATCH --mail-user=baerlachlan@gmail.com

## Cores
CORES=16

## Modules
module load FastQC/0.11.7
module load AdapterRemoval/2.2.1-foss-2016b
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.9-foss-2016b
module load STAR/2.7.0d-foss-2016b
module load Subread/1.5.2-foss-2016b

# ## Reference files zebrafish
# RRNA=/fast/users/a1647910/20200310_rRNADepletion/files/bwa/zebrafish/danRer11
# STAR=
# GTF=

## Reference files mouse
RRNA=/fast/users/a1647910/20200310_rRNADepletion/files/bwa/mouse/mm10
STAR=/data/biorefs/reference_genomes/ensembl-release-98/mus_musculus/star
GTF=/data/biorefs/reference_genomes/ensembl-release-98/mus_musculus/Mus_musculus.GRCm38.98.chr.gtf.gz

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

## Setup for BWA alignment
ALIGNDATABWA=${PROJROOT}/2_alignedDataBwa
mkdir -p ${ALIGNDATABWA}/bam
mkdir -p ${ALIGNDATABWA}/fastq
mkdir -p ${ALIGNDATABWA}/log
mkdir -p ${ALIGNDATABWA}/FastQC

## Setup for STAR alignment
ALIGNDATASTAR=${PROJROOT}/3_alignedDataStar
mkdir -p ${ALIGNDATASTAR}/bam
mkdir -p ${ALIGNDATASTAR}/log
mkdir -p ${ALIGNDATASTAR}/FastQC
mkdir -p ${ALIGNDATASTAR}/featureCounts

##--------------------------------------------------------------------------------------------##
## FastQC on the raw data
##--------------------------------------------------------------------------------------------##

# fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/ERR*.fastq.gz              ###

##--------------------------------------------------------------------------------------------##
## Trimming the merged data
##--------------------------------------------------------------------------------------------##

# for R1 in ${RAWDATA}/fastq/ERR*.fastq.gz                                                        ###
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
# fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastq/ERR*.fastq.gz            ###

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference rRNA
##--------------------------------------------------------------------------------------------##

# ## Aligning and sorting
# for R1 in ${TRIMDATA}/fastq/ERR*.fastq.gz                                                       ###
# do

#   out=${ALIGNDATABWA}/bam/$(basename ${R1%.fastq.gz})
#   echo -e "Output file will be ${out}"

#   ## Align and return reads as .bam
#   bwa mem -t ${CORES} ${RRNA} ${R1} \
#   | samtools view -u -h \
#   | samtools sort -o ${out}.sorted.bam

# done

# ## Run FastQC
# fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATABWA}/FastQC --noextract ${ALIGNDATABWA}/bam/ERR*.bam

# ## Indexing, flagstat, and conversion of unmapped reads to fastq for further alignment
# for BAM in ${ALIGNDATABWA}/bam/ERR*.bam                                                         
# do

#   outbam=${ALIGNDATABWA}/log/$(basename ${BAM%.sorted.bam})
#   outfastq=${ALIGNDATABWA}/fastq/$(basename ${BAM%.sorted.bam})
#   echo -e "Working on ${outbam}"

#   samtools index ${BAM}
#   samtools flagstat ${BAM} > ${outbam}.flagstat

#   echo -e "Now working on ${outfastq}"

#   ## Output only unmapped reads as fastq using -f 4
#   samtools fastq -f 4 -c 6 --threads ${CORES} ${BAM} > ${outfastq}.fastq

# done

# gzip ${ALIGNDATABWA}/fastq/*.fastq

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to the genome
##--------------------------------------------------------------------------------------------##

# ## Aligning, filtering and sorting
# for R1 in ${ALIGNDATABWA}/fastq/ERR*.fastq.gz
# do

#   BNAME=$(basename ${R1%.fastq.gz})
#   echo -e "STAR will align:\t${R1}"

#   STAR \
#     --runThreadN ${CORES} \
#     --genomeDir ${STAR} \
#     --readFilesIn ${R1} \
#     --readFilesCommand gunzip -c \
#     --outFileNamePrefix ${ALIGNDATASTAR}/bam/${BNAME} \
#     --outSAMtype BAM SortedByCoordinate

# done

# # Move the log files into their own folder
# mv ${ALIGNDATASTAR}/bam/*out ${ALIGNDATASTAR}/log
# mv ${ALIGNDATASTAR}/bam/*tab ${ALIGNDATASTAR}/log

# Fastqc and indexing
# for BAM in ${ALIGNDATASTAR}/bam/ERR*.bam
# do
#   fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATASTAR}/FastQC --noextract ${BAM}
#   samtools index ${BAM}
# done

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
samples=`find ${ALIGNDATASTAR}/bam -name "ERR*out.bam" | tr '\n' ' '`

## Extract gtf for featureCounts
zcat ${GTF} > temp.gtf

## Running featureCounts on the sorted bam files
featureCounts -Q 10 \
  -s 0 \
  -T ${CORES} \
  -a temp.gtf \
  -o ${ALIGNDATASTAR}/featureCounts/counts.out ${samples}

## Remove the temporary gtf
rm temp.gtf

## Storing the outputs in a single file
cut -f1,7- ${ALIGNDATASTAR}/featureCounts/counts.out \
| sed 1d > ${ALIGNDATASTAR}/featureCounts/genes.out
