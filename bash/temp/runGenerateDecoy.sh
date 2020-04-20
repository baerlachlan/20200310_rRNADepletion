#!/bin/bash

## Directories
PROJROOT=/fast/users/a1647910/20200310_rRNADepletion
SOFTWARE=/apps/software

sbatch generateDecoyTranscriptome.sh \
	-g ${PROJROOT}/files/Danio_rerio.GRCz11.dna.primary_assembly.fa \
	-t ${PROJROOT}/files/Danio_rerio.GRCz11.cdna.primary_assembly.fa.gz \
	-a ${PROJROOT}/files/Danio_rerio.GRCz11.98.chr.gtf.gz \
	-m ${SOFTWARE}/mashmap/2.0/bin/mashmap \
	-b ${SOFTWARE}/BEDTools/2.25.0-foss-2015b/bin/bedtools \
	-o ${PROJROOT}/files