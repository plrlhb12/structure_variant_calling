#!/bin/bash
# pe-test

file=${1}
svtk pe-test ${file}_short_diploidSV.standardize.vcf ${file}_discfile.gz ${file}_pe
