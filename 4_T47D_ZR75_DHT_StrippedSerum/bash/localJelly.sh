#!/bin/bash

## Cores
CORES=6

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
BWA=${PROJROOT}/3_bwa

## Jellyfish
JELLYFQ=${PROJROOT}/5_jellyfishFq

##----------------------------------------------------------------------------##
##                                Jellyfish                                   ##
##----------------------------------------------------------------------------##

## This script is written to run on a local linux machine.
## Numerous attempts were made at submitting an equivalent script to the HPC1 cluster but kept failing with error OUT_OF_MEMORY.
## The script runs successfully on local machine with the same hardware requests, so it appears jellyfish may corrupt on HPC1.

for k in $(seq 5 10)
  do

    ## Kmer size
    KMER=${k}
    checkAndMake ${JELLYFQ}/k${KMER}

    for R1 in ${BWA}/fastq/*R1.fastq.gz
      do

        R2=${R1%_R1.fastq.gz}_R2.fastq.gz
        echo -e "Counting ${KMER}-mers in R1 file:\t${R1}"
        echo -e "Also counting ${KMER}-mers in R2 file:\t${R2}"

        ## Create output name
        out=${JELLYFQ}/k${KMER}/$(basename ${R1%_R1.fastq.gz})
        echo -e "Kmer counts of both reads will be written to single file ${out}_dumps.txt"

        /usr/bin/jellyfish count -m ${KMER} -s 4G --bf-size 4G -C -t ${CORES} -o ${out}_counts.jf <(zcat ${R1}) <(zcat ${R2})
        /usr/bin/jellyfish dump -c ${out}_counts.jf > ${out}_dumps.txt
        rm ${out}_counts.jf

      done

  done