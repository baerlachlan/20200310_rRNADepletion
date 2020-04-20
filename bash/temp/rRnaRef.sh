#!/bin/bash/

cut -d ";" -f 1 /data/biorefs/rRNA/danio_rerio/danRer11.fa \
| sed -n "1~2p" \
| sed -e "s/>//g" \
| sed  "1iensembl_gene_id" \
> /fast/users/a1647910/20200310_rRNADepletion/files/rRnaRef.txt
