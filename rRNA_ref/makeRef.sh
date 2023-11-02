#!/bin/bash

## This scipt was used to create an rRNA reference consisting of Ensembl rRNA annotations as well as SILVA database sequences
## Initial files:
## 	- SILVA RefNR SSU (18S) rRNA fasta
##	- SILVA RefNR LSU (28S) rRNA fasta
##	- Ensembl 5/5.8S rRNA fasta (generated from R script Ens_rRNA.R

## Define files
SSU="SILVA_138_1_SSU_RefNR.fasta"
LSU="SILVA_138_1_LSU_RefNR.fasta"
ENS="Ens_rRNA.fa"

## First we edit the header of SILVA fasta files:
##	- awk performs the function line by line
##	- first it splits each line by separator " "
##	- then, if the first component of the split has ">" at the start of the line, it prints the first component followed by ";18S_rRNA;rRNA" for example

## Run this script from the makeRef dir

## SSU
awk '{split($0,a," "); if(a[1] ~ /^>/) print a[1]";18S_rRNA;rRNA"; else print; }' ${SSU} > SSU.fa
## LSU
awk '{split($0,a," "); if(a[1] ~ /^>/) print a[1]";28S_rRNA;rRNA"; else print; }' ${LSU} > LSU.fa

## Then concatenate the SSU and LSU files:
cat SSU.fa LSU.fa > SILVA_rRNA.fa

## Then replace all Us with Ts
sed -i '/^>/! {s/U/T/g}' SILVA_rRNA.fa

## This can now be appended to the Ensembl rRNA
cat ${ENS} SILVA_rRNA.fa > rRNA_Ens_SILVA.fa

## Remove intermediary files
rm SSU.fa
rm LSU.fa
rm SILVA_rRNA.fa
