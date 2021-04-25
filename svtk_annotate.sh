#!/bin/bash/
file=${1}
REF=/data/CARD/tprojects/refs/svtk_ref/broad/v1
GENCODE=${REF}/gencode.canonical_pc.gtf.gz
NONCODING=${REF}/noncoding.sort.hg38.bed
svtk annotate --gencode ${GENCODE} --noncoding ${NONCODING} ${file} anno.${file}