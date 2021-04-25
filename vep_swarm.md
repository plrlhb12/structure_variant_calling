
# preprocessing vcf
# get the common from two short svs
file1=psomagen_short_diploidSV.standardize.vcf.gz
file2=jax_short_diploidSV.standardize.vcf.gz

for file in ${file1} ${file2}
do
bcftools index  $file
done

`bcftools isec --threads 8 ${file1} ${file2} -n+2 -p 2common_sv_std`

using one of the output feeding to VEP: /2common_sv_std/0000.txt

# get the common variants from 3 kolf2.1

file1=UNHS_GT19-38445_bcfc1_2.vcf.gz
file2=Jax.KOLF2-ARID2-A2_bcfc1_2.vcf.gz
file3=qual30_pass_annotated_long_Psomagen_KOLF2.1_anno.vcf.gz

for file in ${file1} ${file2} ${file3}
do
bcftools view --threads 56 -e '%QUAL/(INFO/DP) < 2.0 | INFO/DP<10 | INFO/MQ<40' ${file} | bgzip -c > filtered_${file};
tabix -p vcf filtered_${file};
done

`bcftools isec --threads 64 filtered_$file1 filtered_$file2 filtered_$file3 -n+3 -p 3common_snp`

using one of the output to feed to VEP : 3common_snp/0000.txt 


# codes in swarm file

# use the example vcf in the VEP package
`vep -i homo_sapiens_GRCh38.vcf -o test.out --offline --cache --dir_cache $VEP_CACHEDIR --assembly GRCh38 --everything`
has output with WARNING

# test annotating one short sv and took forever, failed 
vep -i psomagen_short_diploidSV.standardize.vcf.gz -o psomagen_sv_vep.vcf --offline --cache --dir_cache $VEP_CACHEDIR --assembly GRCh38 --species human --fasta $VEP_CACHEDIR/GRCh38.fa --everything
`swarm -f run_vep.swarm -t 16 -g 50 --time=8:00:00 --module VEP`
run with swarm, no output and stuck


# annotate the short sv (common from two standardize), commom snp, and long sv
vep -i 2common_sv_std/0000.vcf -o short_sv_vep.vcf --offline --cache --dir_cache $VEP_CACHEDIR --assembly GRCh38 --species human --fasta $VEP_CACHEDIR/GRCh38.fa --force_overwrite --sift b --polyphen b --symbol --af_gnomad --canonical --gencode_basic

vep -i 3common_snp/0000.vcf -o snp_vep.vcf --offline --cache --dir_cache $VEP_CACHEDIR --assembly GRCh38 --species human --fasta $VEP_CACHEDIR/GRCh38.fa --force_overwrite --sift b --polyphen b --symbol --af_gnomad --canonical --gencode_basic

vep -i long_merged_normed_pass.vcf -o long_sv_vep.vcf --offline --cache --dir_cache $VEP_CACHEDIR --assembly GRCh38 --species human --fasta $VEP_CACHEDIR/GRCh38.fa --force_overwrite --sift b --polyphen b --symbol --af_gnomad --canonical --gencode_basic

`swarm -f run_vep_2.swarm -t 8 -g 50 --time=8:00:00 --module VEP`