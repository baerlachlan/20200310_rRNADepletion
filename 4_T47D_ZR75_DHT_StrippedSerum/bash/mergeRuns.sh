#!bin/bash

## This script was merges raw FASTQ files from the T47D_ZR75_DHT_StrippedSerum dataset for each sample that was split over 2 lanes

## Directory
PROJROOT=/hpcfs/users/a1647910/20200310_rRNADepletion/4_T47D_ZR75_DHT_StrippedSerum
UNMERGED=${PROJROOT}/files

## Get sample filenames
sampleList=`find ${UNMERGED}/Hiseq_1/* -name "*.fastq.gz" | tr '\n' ' '`

for i in ${sampleList}
  do

  basename=$(basename ${i})
  parentDir=$(basename $(dirname ${i}))
  file1=${UNMERGED}/Hiseq_1/${parentDir}/${basename}
  file2=${UNMERGED}/Hiseq_2/${parentDir}/${basename}
  merged=${PROJROOT}/0_rawData/fastq/${basename%.gz}

  echo -e "Merging\t${file1}"
  echo -e "with\t${file2}"
  echo -e "Merged file will be ${merged}.gz"

  echo -e "zcat ${file1} ${file2} > ${merged}
  gzip ${merged}"

  echo -e "Merge complete\n"

  done