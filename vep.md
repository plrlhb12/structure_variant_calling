# use vep following hpc instruction
https://hpc.nih.gov/apps/VEP.html

# reference paths: 
# Most of the reference data is available within $VEP_CACHEDIR, but some are available in the /fdb tree
--fasta $VEP_CACHEDIR/GRCh38.fa --species human --assembly GRCh38

# Examples
./vep -i in.vcf -o out.txt -cache -everything
./filter_vep -i out.txt -o out_filtered.txt -filter "[filter_text]"
./vep -i in.vcf -o stdout -cache -check_existing | ./filter_vep -filter "not Existing_variation" -o out.txt

# basic options for vep
# default vep mayoutput one variant for several transcripts of the same gene, 
--cache -i ${file} -o output_${file}  --compress_output bgzip --force_overwrite \
--sift b \
--polygen b \
--symbol \
--vcf 
--af_gnomad \
--canonical \
--gencode_basic \
--pick \ # Pick one line or block of consequence data per variant, including transcript-specific columns.
--assembly GRCh38 2>&1 >> '$log_file'

# VEP options that I can add
--pick: each line for each transcripts
--per-gene
--stats_file [filename]
--summary : Output only a comma-separated list of all observed consequences per variant. Transcript-specific columns will be left blank.
--af_gnomad / --af / --af_1kg / --af_esp : Include allele frequency, not that it will export multiple populations
--canonical : label the varaint is canonical or not (YES/NOT)

# example
Let's say we're only interested in what is considered the canonical transcript for this gene (i.e., only output the representative transcript from mulitple matched isoforms for a varaint)

./vep -i examples/homo_sapiens_GRCh38.vcf --cache --force_overwrite --sift b --canonical --symbol --tab --fields Uploaded_variation,SYMBOL,CANONICAL,SIFT -o STDOUT | ./filter_vep --filter "CANONICAL is YES and SIFT is deleterious"


# some options that I can try later
--coding_only : Only return consequences that fall in the coding regions of transcripts
--protein 
--freq_pop : filter by freq according population gnomAD_NFE/1KG_EUR
--fields: define the columns to output



## plugins that I can try: 
CADD
LOFTEE


# post-run processing using the script of filter_vep instead of vep
http://uswest.ensembl.org/info/docs/tools/vep/script/vep_filter.html

The VEP package includes a tool, filter_vep, to filter results files on a variety of attributes.
It operates on standard, tab-delimited or VCF formatted output (NB only VCF output produced by VEP or in the same format can be used).

# SIFT calls variants either "deleterious" or "tolerated".
./filter_vep -i OUTPUT_FROM_VEP -filter "SIFT is deleterious" | grep -v "##" | head -n5

# writing filters: can use is, <, or, not, mtach, in, exists... in " " 
--ontology --filter "Consequence is coding_sequence_variant"   (mtach all coding consequnces (e.g., missense_variant, synomymous_variant))

--filter "SYMBOL"
--filter "SYMBOL in BRCA1,BRCA2"

--filter "Feature is ENST0000307301"
--filter "Feature in /data/files/motifs_list/txt"
--filter "Feature_type is Transcript"    # get only transcript consequencees

--filter "Exon > 1"
--filter "Consequence match stop"
--filter "Consequence match stream"
--filter "SIFT != tolerated"

--filter " AF < 0.01 or not AF" # finding varaints where the AF is less than 1% or absent
--filter "(AFR_AF < 0.1 or EUR_AF > 0.1) and (EAS_AF < 0.1 and SAS_AF < 0.1)"
--filter "not Exiting_variation" # filter out known variants
--filter "AFR_AF < #EUR_AF"







