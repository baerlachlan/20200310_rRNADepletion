#!/bin/bash

wget -i fileUrls2.txt

mv *.fastq.gz ../0_rawData/fastq
