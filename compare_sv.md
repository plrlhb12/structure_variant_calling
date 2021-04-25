# 16Nov2020: compare SV files

## pengl7 at biowulf in /data/CARD/tprojects/projects/vep_kolf2


### 1. compare short psomagen and short jax sv

### 1.a first compare the short files without svtk standardization
`vcf-compare GT19-38445.diploidSV.vcf.gz KOLF2-ARID2-A2.diploidSV.vcf.gz > 2short_sv_vcf_compare`

VN      3137    KOLF2-ARID2-A2.diploidSV.vcf.gz (31.4%)
VN      6352    GT19-38445.diploidSV.vcf.gz (48.1%)
VN      6860    GT19-38445.diploidSV.vcf.gz (51.9%)     KOLF2-ARID2-A2.diploidSV.vcf.gz (68.6%)

### 1.b compare the short files after standardization
`vcf-compare psomagen_short_diploidSV.standardize.vcf.gz jax_short_diploidSV.standardize.vcf.gz > 2short_sv_std_compare`

VN      1294    jax_short_diploidSV.standardize.vcf.gz (20.4%)
VN      3584    psomagen_short_diploidSV.standardize.vcf.gz (41.6%)
VN      5034    jax_short_diploidSV.standardize.vcf.gz (79.6%)  psomagen_short_diploidSV.standardize.vcf.gz (58.4%)


## 2. Compare long to short standardized sv
### to use vcf-compare, the file must be block ziped and indexed
### choose the long sv file which passeed filter"pass"

prefix=long_merged_normed_pass
bcftools view ${prefix}.vcf -Oz -o ${prefix}.vcf.gz
bcftools index ${prefix}.vcf.gz

prefix=2common_sv_std/0000
bcftools view ${prefix}.vcf -Oz -o common_2sv_std.vcf.gz
bcftools index common_2sv_std.vcf.gz

### 2.1 compare long sv to short sv_common, found no overlapping at all

`vcf-compare long_merged_normed_pass.vcf.gz common_2sv_std.vcf.gz > long_short_sv_compare`

VN      5026    common_2sv_std.vcf.gz (100.0%)
VN      9802    long_merged_normed_pass.vcf.gz (100.0%)

### 2.2 compare long sv to the individual standardized sv, no overlapping
`vcf-compare long_merged_normed_pass.vcf.gz psomagen_short_diploidSV.standardize.vcf.gz > long_psomagen_short_sv_std_compareq`
nonoverlap

`vcf-compare long_merged_normed_pass.vcf.gz jax_short_diploidSV.standardize.vcf.gz > long_jax_short_sv_std_compare`
non-overlap

### 2.3 compare long sv to the individual unstandardized sv
`vcf-compare long_merged_normed_pass.vcf.gz GT19-38445.diploidSV.vcf.gz > long_psomagen_short_sv_compare`
VN      1940    GT19-38445.diploidSV.vcf.gz (14.7%)     long_merged_normed_pass.vcf.gz (19.8%)
VN      7862    long_merged_normed_pass.vcf.gz (80.2%)
VN      11272   GT19-38445.diploidSV.vcf.gz (85.3%)

`vcf-compare long_merged_normed_pass.vcf.gz KOLF2-ARID2-A2.diploidSV.vcf.gz > long_jax_short_sv_compare`
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      1817    KOLF2-ARID2-A2.diploidSV.vcf.gz (18.2%) long_merged_normed_pass.vcf.gz (18.5%)
VN      7985    long_merged_normed_pass.vcf.gz (81.5%)
VN      8180    KOLF2-ARID2-A2.diploidSV.vcf.gz (81.8%)


## 3. compare long sv to short unstandardized sv
### re-generate the intersect between two unstandardized short sv 
# first pass "PASS", need to index again because the new generated file don't have index yet
# then isec
```
for sample in GT19-38445.diploidSV.vcf.gz KOLF2-ARID2-A2.diploidSV.vcf.gz 
do
#bcftools view -f 'PASS' ${sample} -Oz -o pass_${sample}
bcftools index pass_${sample}
done
```

`bcftools isec pass_GT19-38445.diploidSV.vcf.gz pass_KOLF2-ARID2-A2.diploidSV.vcf.gz -n+2 -p 2common_sv`

```
prefix=2common_sv/0000
bcftools view ${prefix}.vcf -Oz -o common_2sv.vcf.gz
bcftools index common_2sv.vcf.gz
```
`vcf-compare long_merged_normed_pass.vcf.gz common_2sv.vcf.gz > long_short_sv_compare`

#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      1626    common_2sv.vcf.gz (28.8%)       long_merged_normed_pass.vcf.gz (16.6%)
VN      4011    common_2sv.vcf.gz (71.2%)
VN      8176    long_merged_normed_pass.vcf.gz (83.4%)

`vcf-compare pass_GT19-38445.diploidSV.vcf.gz pass_KOLF2-ARID2-A2.diploidSV.vcf.gz > 2sv_pass_compare`

#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      2320    pass_KOLF2-ARID2-A2.diploidSV.vcf.gz (28.2%)
VN      5304    pass_GT19-38445.diploidSV.vcf.gz (47.3%)
VN      5903    pass_GT19-38445.diploidSV.vcf.gz (52.7%)        pass_KOLF2-ARID2-A2.diploidSV.vcf.gz (71.8%)



# check the NIST 
# precessing longrange sv files
# merge del file together with large sv togther in the long ranger vcf, -a allow overlapping
file1=dels.vcf.gz
file2=large_svs.vcf.gz
bcftools concat -a ${file1} ${file2} -Oz -o long_merged_nist.vcf.gz
bcftools norm -m-any long_merged_nist.vcf.gz -Oz -o long_merged_normed_nist.vcf.gz
bcftools index long_merged_normed_nist.vcf.gz

for file in long_merged_normed_nist.vcf.gz NIST-reference-sample.diploidSV.vcf.gz
do
bcftools view -f 'PASS' ${file} -Oz -o pass_${file}
bcftools index pass_${file}
done

for file in long_merged_normed_nist.vcf.gz NIST-reference-sample.diploidSV.vcf.gz
do
echo ${file}
echo $(zcat ${file} | wc -l)
echo $(zcat pass_${file} | wc -l)
done

long_merged_normed_nist.vcf.gz
9940
5349
NIST-reference-sample.diploidSV.vcf.gz
15598
13705

`vcf-compare long_merged_normed_nist.vcf.gz NIST-reference-sample.diploidSV.vcf.gz > nist_compare`

#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      1026    NIST-reference-sample.diploidSV.vcf.gz (8.5%)   long_merged_normed_nist.vcf.gz (10.4%)
VN      8810    long_merged_normed_nist.vcf.gz (89.6%)
VN      11062   NIST-reference-sample.diploidSV.vcf.gz (91.5%)


`vcf-compare pass_long_merged_normed_nist.vcf.gz pass_NIST-reference-sample.diploidSV.vcf.gz > pass_nist_compare`

#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN      1088    pass_NIST-reference-sample.diploidSV.vcf.gz (10.7%)     pass_long_merged_normed_nist.vcf.gz (20.7%)
VN      4162    pass_long_merged_normed_nist.vcf.gz (79.3%)
VN      9126    pass_NIST-reference-sample.diploidSV.vcf.gz (89.3%)
