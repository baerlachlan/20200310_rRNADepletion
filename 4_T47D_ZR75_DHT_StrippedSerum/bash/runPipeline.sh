#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --ntasks-per-core=2
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH -o /hpcfs/users/a1647910/20200310_rRNADepletion/4_T47D_ZR75_DHT_StrippedSerum/slurm/%x_%j.out
#SBATCH -e /hpcfs/users/a1647910/20200310_rRNADepletion/4_T47D_ZR75_DHT_StrippedSerum/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=8

## HPC modules
module load arch/haswell
module load arch/skylake

## Modules
module load FastQC/0.11.7
module load STAR/2.7.0d-foss-2016b
module load SAMtools/1.10-foss-2016b
module load cutadapt/1.14-foss-2016b-Python-2.7.13
module load Subread/1.5.2-foss-2016b
module load BWA/0.7.15-foss-2017a
module load Salmon/1.1.0-foss-2016b

## Function for checking directories
checkAndMake () {
  echo "Checking if $1 exists"
  if [[ ! -d $1 ]]
    then 
      echo "Creating $1"
      mkdir -p $1
  fi
    
  if [[ -d $1 ]]
    then
      echo "Found $1"
    else
      echo "$1 could not be created or found"
      exit 1
  fi  
}

## Directories
PROJROOT=/hpcfs/users/a1647910/20200310_rRNADepletion/4_T47D_ZR75_DHT_StrippedSerum
REFS=/hpcfs/users/a1647910/refs
if [[ ! -d ${REFS} ]]
then
  echo "Couldn't find ${REFS}"
  exit 1
fi
GTF=${REFS}/ensembl-release-101/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
if [[ ! -f ${GTF} ]]
then
  echo "Couldn't find ${GTF}"
  exit 1
fi
RRNA=${REFS}/rRNA/homo_sapiens/ensembl-release-101/bwa/rRNA_EnsSilva.fa
if [[ ! -f ${RRNA} ]]
then
  echo "Couldn't find ${RRNA}"
  exit 1
fi

# Raw Data
RAWDIR=${PROJROOT}/0_rawData
checkAndMake ${RAWDIR}
checkAndMake ${RAWDIR}/FastQC

## Trimmed Data
TRIMDIR=${PROJROOT}/1_trimmedData
checkAndMake ${TRIMDIR}/fastq
checkAndMake ${TRIMDIR}/FastQC
checkAndMake ${TRIMDIR}/log

## Aligned Data
ALIGNDIR=${PROJROOT}/2_alignedData
checkAndMake ${ALIGNDIR}
checkAndMake ${ALIGNDIR}/bam
checkAndMake ${ALIGNDIR}/FastQC
checkAndMake ${ALIGNDIR}/log
checkAndMake ${ALIGNDIR}/featureCounts

## BWA
BWA=${PROJROOT}/3_bwa
checkAndMake ${BWA}/bam
checkAndMake ${BWA}/fastq
checkAndMake ${BWA}/log
checkAndMake ${BWA}/FastQC

## STAR 2-pass
STAR2=${PROJROOT}/4_star2pass
checkAndMake ${STAR2}/bam
checkAndMake ${STAR2}/fastq
checkAndMake ${STAR2}/FastQC
checkAndMake ${STAR2}/log
checkAndMake ${STAR2}/featureCounts

## Jellyfish
JELLYFQ=${PROJROOT}/5_jellyfishFq

echo "All directories checked and created"

##----------------------------------------------------------------------------##
##                              Initial FastQC                                ##
##----------------------------------------------------------------------------##

# fastqc -t ${CORES} -o ${RAWDIR}/FastQC --noextract ${RAWDIR}/fastq/*fastq.gz

##----------------------------------------------------------------------------##
##                             Identify Adapters                              ##
##----------------------------------------------------------------------------##

## This didn't work for this dataset, but keeping it here for future reference

# ## bbmerge script
# BBMERGE=/hpcfs/users/a1647910/packages/bbmap/bbmerge.sh
# ## Adapter folder
# checkAndMake ${RAWDIR}/adapters

# for R1 in ${RAWDIR}/fastq/*R1*.fastq.gz
# do

#   R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz
#   echo -e "R1 file should be ${R1}"
#   echo -e "R2 file should be ${R2}"
#   out=${RAWDIR}/adapters/adapters.fa

#   ${BBMERGE} in1=${R1} in2=${R2} outa=${out}
#   exit

# done

##----------------------------------------------------------------------------##
##                                 Trimming                                   ##
##----------------------------------------------------------------------------##

# for R1 in ${RAWDIR}/fastq/*R1*.fastq.gz
#   do
#     R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz
#     echo -e "The R1 file should be ${R1}"
#     echo -e "The R2 file should be ${R2}"

#     ## Create output filenames
#     out1=${TRIMDIR}/fastq/$(basename $R1)
#     out2=${TRIMDIR}/fastq/$(basename $R2)
#     BNAME=${TRIMDIR}/fastq/$(basename ${R1%_R1_001.fastq.gz})
#     echo -e "Output file 1 will be\t${out1}"
#     echo -e "Output file 2 will be\t${out2}"
#     echo -e "Trimming:\t${BNAME}"

#     LOG=${TRIMDIR}/log/$(basename ${BNAME}).info
#     echo -e "Trimming info will be written to ${LOG}"

#     cutadapt \
#       -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#       -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#       -o ${out1} \
#       -p ${out2} \
#       -m 35 \
#       --trim-n \
#       --max-n=1 \
#       --nextseq-trim=30 \
#       ${R1} \
#       ${R2} > ${LOG}

#   done

# fastqc -t ${CORES} -o ${TRIMDIR}/FastQC --noextract ${TRIMDIR}/fastq/*fastq.gz

##----------------------------------------------------------------------------##
##                              STAR alignment                                ##
##----------------------------------------------------------------------------##

# ## Aligning, filtering and sorting
# for R1 in ${TRIMDIR}/fastq/*R1*.fastq.gz
# do

#   BNAME=$(basename ${R1%_R1_001.fastq.gz})
#   R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz
#   echo -e "STAR will align:\t${R1}"
#   echo -e "STAR will also align:\t${R2}"

#   STAR \
#     --runThreadN ${CORES} \
#     --genomeDir ${REFS}/ensembl-release-101/homo_sapiens/star \
#     --readFilesIn ${R1} ${R2} \
#     --readFilesCommand gunzip -c \
#     --outFileNamePrefix ${ALIGNDIR}/bam/${BNAME} \
#     --outSAMtype BAM SortedByCoordinate

# done

# ## Move the log files into their own folder
# mv ${ALIGNDIR}/bam/*out ${ALIGNDIR}/log
# mv ${ALIGNDIR}/bam/*tab ${ALIGNDIR}/log

# ## Fastqc and indexing
# for BAM in ${ALIGNDIR}/bam/*.bam
# do
#   fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDIR}/FastQC --noextract ${BAM}
#   samtools index ${BAM}
# done

##----------------------------------------------------------------------------##
##                               featureCounts                                ##
##----------------------------------------------------------------------------##

# ## Feature Counts - obtaining all sorted bam files
# sampleList=`find ${ALIGNDIR}/bam -name "*out.bam" | tr '\n' ' '`

# ## Extract gtf for featureCounts
# zcat ${GTF} > temp.gtf

# ## Running featureCounts on the sorted bam files
# featureCounts -Q 10 \
#   -s 2 \
#   -T ${CORES} \
#   -p \
#   --fracOverlap 1 \
#   -a temp.gtf \
#   -o ${ALIGNDIR}/featureCounts/counts.out ${sampleList}

# ## Remove the temporary gtf
# rm temp.gtf

#  ## Storing the output in a single file
# cut -f1,7- ${ALIGNDIR}/featureCounts/counts.out | \
#   sed 1d > ${ALIGNDIR}/featureCounts/genes.out

##----------------------------------------------------------------------------##
##                         rRNA proportions with BWA                          ##
##----------------------------------------------------------------------------##

# ## Aligning and sorting
# for R1 in ${TRIMDIR}/fastq/*R1*.fastq.gz                                                       
# do

#   BNAME=$(basename ${R1%_R1_001.fastq.gz})
#   R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz
#   out=${BWA}/bam/${BNAME}.sorted.bam
#   echo -e "bwa will align:\t${R1}"
#   echo -e "bwa will also align:\t${R2}"
#   echo -e "Output file will be:\t${out}"

#   ## Align and return reads as .bam
#   bwa mem -t ${CORES} ${RRNA} ${R1} ${R2} \
#   | samtools view -u -h \
#   | samtools sort -o ${out}

# done

# ## Run FastQC
# fastqc -t ${CORES} -f bam_mapped -o ${BWA}/FastQC --noextract ${BWA}/bam/*.bam

# ## Indexing, flagstat, and conversion of unmapped reads to fastq for further alignment
# for BAM in ${BWA}/bam/*.bam                                                         
# do

#   outlog=${BWA}/log/$(basename ${BAM%.sorted.bam})
#   outbam=${BWA}/bam/$(basename ${BAM%.sorted.bam})
#   outfastq=${BWA}/fastq/$(basename ${BAM%.sorted.bam})

#   samtools index ${BAM}
#   samtools flagstat ${BAM} > ${outlog}.flagstat 
#   awk -F '[(|%]' 'NR == 5 {print $2}' ${outlog}.flagstat > ${outlog}.mapped

#   ## Output only reads where both mate pairs were not mapped as fastq
#   ## Reads need to be sorted by name (-n) to retain pair structure
#   ## This code was found from https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd
#   samtools view -u -h -f 12 -F 256 ${BAM} \
#   | samtools sort -n \
#   | samtools fastq --threads ${CORES} -c 6 -1 ${outfastq}_R1.fastq.gz -2 ${outfastq}_R2.fastq.gz

# done

# ## Concatenate all mapping stats including the filenames and mapped proportions
# grep -H "" ${BWA}/log/*.mapped > ${BWA}/log/samples.mapped.all

##----------------------------------------------------------------------------##
##                         STAR second-pass alignment                         ##
##----------------------------------------------------------------------------##

## The reason for running STAR a second time is to fully remove any rRNA sequences 
## that may have been missed in previous steps.
## BWA was used to remove known rRNA sequences sourced from the .gtf file as well
## as SILVA.
## Running STAR a second time in the absence of reads that map to these known sequences,
## followed by extraction of only mapped reads will exclude any rRNA reads that are not 
## known from the resulting fastq files.
## Of course some non-rRNA reads that do not map will be lost too, but it will stop 
## any bias in kmer counts due to the presence of rRNA.

# ## Aligning, filtering and sorting
# for R1 in ${BWA}/fastq/*R1*.fastq.gz  
# do

#   BNAME=$(basename ${R1%_R1.fastq.gz})
#   R2=${R1%_R1.fastq.gz}_R2.fastq.gz
#   echo -e "STAR will align:\t${R1}"
#   echo -e "STAR will also align:\t${R2}"

#   STAR \
#     --runThreadN ${CORES} \
#     --genomeDir ${REFS}/ensembl-release-101/homo_sapiens/star \
#     --readFilesIn ${R1} ${R2} \
#     --readFilesCommand gunzip -c \
#     --outFileNamePrefix ${STAR2}/bam/${BNAME} \
#     --outSAMtype BAM SortedByCoordinate

# done

# ## Move the log files into their own folder
# mv ${STAR2}/bam/*out ${STAR2}/log
# mv ${STAR2}/bam/*tab ${STAR2}/log

# ## Fastqc and indexing
# for BAM in ${STAR2}/bam/*.bam
# do

#   fastqc -t ${CORES} -f bam_mapped -o ${STAR2}/FastQC --noextract ${BAM}
#   samtools index ${BAM}

#   ## Conversion of mapped reads to fastq for kmer counting
#   ## SAM Flag 256 removes any reads not in the primary alignment
#   outfastq=${STAR2}/fastq/$(basename ${BAM%Aligned.sortedByCoord.out.bam})
#   samtools fastq -F 256 --threads ${CORES} -c 6 -1 ${outfastq}_R1.fastq.gz -2 ${outfastq}_R2.fastq.gz ${BAM}

# done

##----------------------------------------------------------------------------##
##                         featureCounts second-pass                          ##
##----------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
sampleList=`find ${STAR2}/bam -name "*out.bam" | tr '\n' ' '`

## Extract gtf for featureCounts
zcat ${GTF} > temp.gtf

## Running featureCounts on the sorted bam files
featureCounts -Q 10 \
  -s 2 \
  -T ${CORES} \
  -p \
  --fracOverlap 1 \
  -a temp.gtf \
  -o ${STAR2}/featureCounts/counts.out ${sampleList}

## Remove the temporary gtf
rm temp.gtf

 ## Storing the output in a single file
cut -f1,7- ${STAR2}/featureCounts/counts.out | \
  sed 1d > ${STAR2}/featureCounts/genes.out

##----------------------------------------------------------------------------##
##                                   Salmon                                   ##
##----------------------------------------------------------------------------##

# ## Run salmon
# for R1 in ${TRIMDIR}/fastq/*R1.fastq.gz
# do

#   BNAME=$(basename ${R1%_R1.fastq.gz})
#   R2=${R1%_R1.fastq.gz}_R2.fastq.gz
#   echo -e "file 1 is ${R1}"
#   echo -e "file 2 is ${R2}"

#   salmon quant \
#     -i ${PROJROOT}/5_salmon/salmon_index_drerio98 \
#     -l A \
#     -1 ${R1} \
#     -2 ${R2} \
#     --validateMappings \
#     --threads ${CORES} \
#     -o ${SALMON}/${BNAME}

# done

##----------------------------------------------------------------------------##
##                                Jellyfish                                   ##
##----------------------------------------------------------------------------##

## Due to issues getting jellyfish to run on the new HPC1 system, jellyfish was run locally from script localJelly.sh
## Jellyfish was crashing on HPC1 with error "OUT_OF_MEMORY" even when up to 128gb was requested
## When the same script was run locally, jellyfish used no more than 4gb memory