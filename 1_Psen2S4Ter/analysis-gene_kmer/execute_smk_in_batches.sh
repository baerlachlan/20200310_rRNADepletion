#!/bin/bash

for num in {1..100}
do
	snakemake --use-conda --rerun-incomplete --scheduler greedy --batch all=$num/100
done