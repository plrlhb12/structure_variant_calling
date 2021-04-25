# svtk for 4 KOLF2.1 cell lines

### Nov 10th, 2020

This the comments from Anni "I’ve attached my script and I believe all the necessary external files but let me know if I forgot any. Basically in the script I ran manta version 1.6 to call the variants but had to convert the output back to 1.4 versioning since the newer versioning didn’t handle inversions correctly from what I remember. Then following that I was able to filter and standardize".

These are comments from Kim "For our current PD analysis we actually ran an **older version** of manta with default settings, we are then running gnomad-SV with the manta output that runs multiple filtering modules within the pipeline."

### Biowulf has already installed svtk. Although I have one installed locally, it would betteer to run biowulf incase later I have to apply bcftools

## 1.  check the diff between vcf before and after convertion

```
cd compare
prefix=kolf2.1.4
bcftools view ${prefix}.vcf -Oz -o ${prefix}.vcf.gz
bcftools index ${prefix}.vcf.gz

file1=kolf2
file2=kolf2.1.4
bcftools isec -p compare_manta1.6vs1.4 ${file1}.vcf.gz ${file2}.vcf.gz
vcf-compare ${file1}.vcf.gz ${file2}.vcf.gz | grep ^VN | cut -f 2-
```

VN	134	kolf2.1.4.vcf.gz (1.4%)
VN	426	kolf2.vcf.gz (4.3%)
VN	9571	kolf2.1.4.vcf.gz (98.6%)	kolf2.vcf.gz (95.7%)


## 2. copied the required files to biowulf
cd /data/CARD/tprojects/projects/sv-manta
mkdir svtk_kolf2
WKD='/data/CARD/tprojects/projects/sv-manta/svtk_kolf2'
REF='/data/CARD/tprojects/refs/manta'

### 2.1  copy the sv callsets from google bucket
moudle load google-cloud-sdk
gcloud init

#### kolf2.1 doesn't have the long read data, while kolf2-C1 has
for file in dels.vcf.gz dels.vcf.gz.tbi large_svs.vcf.gz large_svs.vcf.gz.tbi
do 
gsutil cp -r gs://singlecellindi/WGS/2004UAHS-0254/hg38-variants/longranger/KOLF2-ARID2-A2/outs/${file} psomagen_longranger_${file}
done 

### 2.2 copy the origian vcf output from manta, also need bam and its index file for svtk's other usages
for file in diploidSV.vcf.gz diploidSV.vcf.gz.tbi
do 
gsutil cp -r gs://singlecellindi/WGS/2001UNHS-0021/structural-variants/GT19-38445.${file} psomagen_short_${file}
gsutil cp -r gs://singlecellindi/WGS/Jax/IlluminaWGS/structural-variants-after-svtk/KOLF2-ARID2-A2.${file} jax_short_${file}
done 

for file in bai bam
do
gsutil cp -r gs://singlecellindi/WGS/2001UNHS-0021/hg38/bams/GT19-38445.${file} psomagen.${file}
gsutil cp -r gs://singlecellindi/WGS/Jax/IlluminaWGS/hg38/bams/KOLF2-ARID2-A2.${file} jax.${file}
done

### 2.3 copy the vcf already standardized by svtk
for file in diploidSV.standardize.vcf
do
gsutil cp -r gs://singlecellindi/WGS/2001UNHS-0021/structural-variants/GT19-38445.${file} psomagen_short_${file}
gsutil cp -r gs://singlecellindi/WGS/Jax/IlluminaWGS/structural-variants-after-svtk/KOLF2-ARID2-A2.${file} jax_short_${file}
done

## 3. Prepreoceessing before running svtk
# I already have the svtk for short read data and can skip it

gunzip -c long_merged_normed.vcf.gz > long_merged_normed.vcf

for sample in long_merged_normed.vcf
do
(grep "^#" ${sample}; awk '$7 == "PASS"' ${sample}) > ${sample}_pass
awk '{{ gsub(/chr/, ""); print}}' ${sample}_pass > edit.${sample}
done 

# check for inversions
for sample in edit.long_merged_normed.vcf 
do
awk '$5 == "<INV>"' ${sample} | wc -l
done

3


# 4. apply svtk to my longranger output sv, but not work 
using the source of manta although mine is longranger
Usage: svtk input output source
options of sources [delly,lumpy,manta,wham,melt]
`svtk standardize ${sample}.diploidSV.1.4.edit.vcf ${sample}.diploidSV.standardize.vcf manta`
```
sample='edit.long_merged_normed.vcf' 
svtk standardize ${sample} standardized.${sample} manta
```
error: 
only finished 1k lines
File "/opt/svtk/svtk/standardize/std_manta.py", line 95, in standardize_info
 std_rec.info['STRANDS'] = strands
UnboundLocalError: local variable 'strands' referenced before assignment



# compare two standardized short sv vcf and edited long vcf
file1=edit.long_merged_normed.vcf
file2=jax_short_diploidSV.standardize.vcf
file3=psomagen_short_diploidSV.standardize.vcf

for file in $file1 $file2 $file3
do
bcftools view ${file} -Oz -o ${file}.gz
bcftools index -t ${file}.gz
done

`vcf-compare ${file2}.gz ${file3}.gz > compare_jax_psomagen`  # only 59% common
`grep ^VN compare_jax_psomagen | cut -f 2-`

######################################## below are just playing svtk
### try svtk vcfcluster
First generatee a file list to contain the two standarideed vcf, the first sample is from jax, the second is from psomagen

Intersect SV called by PE/SR-based algorithms.

Paired-end and split-read callers provide a reasonably precise estimation of
an SV breakpoint. This program identifies variant calls that fall within
the expected margin of error made by these programs and clusters them together.
The cluster distance defaults to 500 bp but it is recommended to use the
maximum individual clustering distance across the libraries being analyzed.(Generally median + 7 * MAD)

```
vim file_list
svtk vcfcluster file_list clustered.vcf
```

for file in clustered.vcf
do
bcftools sort ${file} -Oz -o sort.${file}.gz  # after cluster have to sort the file, gzip, and index again
bcftools index -t sort.${file}.gz
done

`vcf-compare ${file3}.gz sort.clustered.vcf.gz | grep ^VN | cut -f 2-`

# check statistics
svtk count-svtypes [-h] [--no-header] [--total-obs] [--total-variants] vcf [fout]

for file in $file1 $file2 $file3
do
svtk count-svtypes --no-header --total-obs --total-variants ${file} ${file}.stats
done

# check coverage using bed files
bam=
svtk bincov bam chr cov_out(bed file)
for file in $file2 $file2
svtk bincov bam chr 1 cov_out ${file}.bed

for file in psomagen jax
do
svtk bincov ${file}.bam chr 1-22 ${file}.bed
done

svtk bincov jax.bam chr1 cov_out_chr1.bed

# PE/SR analysis: Collect split read and discordant pair data from a bam alignment. Generate split file and discordance file
svtk collect-pesr [-h] [--index-dir INDEX_DIR] [-r REGION] [-z] bam sample-name-to-use-as-id splitfile discfile

for file in psomagen jax
do
svtk collect-pesr ${file}.bam ${file} ${file}_splitfile ${file}_discfile
done

# calulate enrichnent of cliped reads at sv breakpoints
svtk sr-test vcf splitfile fout

# make sure the file is sorted
for file in psomagen jax
do
sort -k1,1 -k2,2n ${file}_splitfile | bgzip > ${file}_splitfile.gz
tabix -p vcf ${file}_splitfile.gz
done

`will it be working if i use bcftools sort?`

# it took 5 hours to be finished one task
# so use swarm to do the job
for file in psomagen jax
do
svtk sr-test ${file}_short_diploidSV.standardize.vcf ${file}_splitfile.gz ${file}_sr
done

# pe-test: this is very quick
for file in psomagen jax
do
(sort -k1,1 -k2,2n ${file}_discfile | bgzip > ${file}_discfile.gz); tabix -p vcf ${file}_discfile.gz
done

svtk pe-test ${file}_short_diploidSV.standardize.vcf ${file}_discfile.gz ${file}_pe

# download the associated files from broad institute
gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/

# need to add / (represent folder instead of file) aftere v1, otherwise not work
gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/ /data/CARD/tprojects/refs/svtk_ref/broad/

# under the folder of /data/CARD/tprojects/refs/svtk_ref/broad/v1
REF=/data/CARD/tprojects/refs/svtk_ref/broad/v1
CYTOBANDS=${REF}/cytobands_hg38.bed.gz
MEI_BED=${REF}/mei_hg38.bed.gz

# first resolve complex variant of breakpoints before annoate: Resolve complex SV from inversion/translocation breakpoints and CNV intervals.
Resolve complex SV from inversion/translocation breakpoints and CNV intervals.

positional arguments:
  raw                   Filtered breakpoints and CNV intervals.
  resolved              Resolved simple and complex variants.

# 
for file in psomagen jax
do
svtk resolve --discfile ${file}_disctfile ${file}_short_diploidSV.standardize.vcf ${file}.resolved
done

svtk resolve --discfile psomagen_disctfile --mei-bed $MEI_BED --cytobands $CYTOBANDS psomagen_short_diploidSV.standardize.vcf psomagen

FileNotFoundError: [Errno 2] No such file or directory: 'vcf-sort': 'vcf-sort'

`will it be working if i use bcftools sort or vcf-sort?`

$ svtk resolve --discfile psomagen_disctfile --mei-bed $MEI_BED --cytobands $CYTOBANDS psomagen_short_diploidSV.standardize.vcf psomagen
svtk resolve @ 14:06:59: starting variant resolution.
Traceback (most recent call last):
  File "/opt/conda/envs/app/bin/svtk", line 7, in <module>
    exec(compile(f.read(), __file__, 'exec'))
  File "/opt/svtk/scripts/svtk", line 67, in <module>
    main()
  File "/opt/svtk/scripts/svtk", line 64, in main
    getattr(cli, command)(sys.argv[2:])
  File "/opt/svtk/svtk/cli/resolve.py", line 451, in main
    stdout=args.resolved)
  File "/opt/conda/envs/app/lib/python3.7/subprocess.py", line 800, in __init__
    restore_signals, start_new_session)
  File "/opt/conda/envs/app/lib/python3.7/subprocess.py", line 1551, in _execute_child
    raise child_exception_type(errno_num, err_msg, err_filename)
FileNotFoundError: [Errno 2] No such file or directory: 'vcf-sort': 'vcf-sort'

# annotation: Annotate resolved SV with genic effects and noncoding hits.
# this is the latest version of gtf I downloaded but will not be used in this case
# use the gencode and noncode files from broad
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz -O /data/CARD/tprojects/refs/gencode

REF=/data/CARD/tprojects/refs/svtk_ref/broad/v1
GENCODE=${REF}/gencode.canonical_pc.gtf.gz
NONCODING=${REF}/noncoding.sort.hg38.bed
file=psomagen
# usage: svtk annotate [-h] [--gencode GENCODE] [--noncoding NONCODING] vcf annotated_vcf

$file2 $file3
for file in $file2
do
svtk annotate --gencode ${GTF} ${file} ${file}.anno
done

file=long_merged_normed.vcf
svtk annotate --gencode ${GENCODE} --noncoding ${NONCODING} ${file) anno.${file}

svtk annotate --gencode ${GENCODE} ${file}_short_diploidSV.standardize.vcf ${file}.anno

# change the annotation of CHROM by adding chr before numbers
(grep '^#' ${file}_short_diploidSV.standardize.vcf; grep -v '^#' ${file}_short_diploidSV.standardize.vcf | sed 's/^/chr/g') > ${file}.add_chr

still not work because [W::vcf_parse] Contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)

Error message was:
pybedtools.helpers.BEDToolsError: 
Command was:

	bedtools intersect -wb -wa -b /data/CARD/tprojects/refs/svtk_ref/broad/v1/gencode.canonical_pc.gtf.gz -a /tmp/pybedtools.nb0yw07l.tmp

[W::sam_read1] Parse error at line 1
***** WARNING: File /tmp/pybedtools.eukek45u.tmp has inconsistent naming convention for record:
1	778704	778705	MantaINS_12_0_0_0_0_0	DEL	+-

