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
PROJROOT=/hpcfs/users/a1647910/20200910_Hypoxia
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

    for R1 in ${BWA}/fastq/*.fastq.gz
      do

        echo -e "Counting ${k}-mers in file ${R1}"

        ## Create output name
        out=${JELLYFQ}/k${KMER}/$(basename ${R1%.fastq.gz})
        echo -e "Kmer counts will be written to file ${out}_dumps.txt"

        /usr/bin/jellyfish count -m ${KMER} -s 4G --bf-size 4G -C -t ${CORES} -o ${out}_counts.jf <(zcat ${R1})
        /usr/bin/jellyfish dump -c ${out}_counts.jf > ${out}_dumps.txt
        rm ${out}_counts.jf

      done

  done