#!/bin/bash

wget -i fileUrls.txt

mv *.fastq.gz ../0_rawData/fastq
