echo -n KOLF2-ARID2-A2 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/KOLF2-ARID2-A2.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/KOLF2-ARID2-A2.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/KOLF2-ARID2-A2/KOLF2-ARID2-A2.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/KOLF2-ARID2-A2     --labels=pipe=manta,sample=manta,cohort=jax

echo -n KUCG3-C1 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/KUCG3-C1.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/KUCG3-C1.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/KUCG3-C1/KUCG3-C1.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/KUCG3-C1     --labels=pipe=manta,sample=manta,cohort=jax

echo -n LNGPI1-C1 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/LNGPI1-C1.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/LNGPI1-C1.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/LNGPI1-C1/LNGPI1-C1.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/LNGPI1-C1     --labels=pipe=manta,sample=manta,cohort=jax

echo -n NCRM1-C6 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/NCRM1-C6.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/NCRM1-C6.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM1-C6/NCRM1-C6.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/NCRM1-C6     --labels=pipe=manta,sample=manta,cohort=jax

echo -n NCRM5-C5 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/NCRM5-C5.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/NCRM5-C5.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/NCRM5-C5/NCRM5-C5.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/NCRM5-C5     --labels=pipe=manta,sample=manta,cohort=jax

echo -n NN0003932-C3 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/NN0003932-C3.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/NN0003932-C3.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0003932-C3/NN0003932-C3.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/NN0003932-C3     --labels=pipe=manta,sample=manta,cohort=jax

echo -n NN0004297-C1 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/NN0004297-C1.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/NN0004297-C1.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/NN0004297-C1/NN0004297-C1.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/NN0004297-C1     --labels=pipe=manta,sample=manta,cohort=jax

echo -n PGP1-C2 OPID=
gcloud alpha genomics pipelines run     --project singlecellseq      --pipeline-file Manta_v1.6.yaml     --zones us-central1-f     --inputs EXECTOOLmanta=gs://test-7cee72c0e768/Jax/tools/manta-1.6.0.centos6_x86_64.tar.gz     --inputs REFSEQ=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     --inputs REFSEQINDEX=gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai     --inputs INCRAM=gs://test-7cee72c0e768/Jax/hg38/crams/PGP1-C2.cram      --inputs INCRAMINDEX=gs://test-7cee72c0e768/Jax/hg38/crams/PGP1-C2.cram.crai     --outputs OUTSUMMARY=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.alignmentStatsSummary.txt     --outputs OUTSMALLINDELS=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.candidateSmallIndels.vcf.gz     --outputs OUTSMALLINDELSINDEX=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.candidateSmallIndels.vcf.gz.tbi     --outputs OUTSV=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.candidateSV.vcf.gz     --outputs OUTSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.candidateSV.vcf.gz.tbi     --outputs OUTDIPLOIDSV=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.diploidSV.vcf.gz     --outputs OUTDIPLOIDSVINDEX=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.diploidSV.vcf.gz.tbi     --outputs OUTSTATSTSV=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.svCandidateGenerationStats.tsv     --outputs OUTSTATSXML=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.svCandidateGenerationStats.xml     --outputs OUTGRAPHTSV=gs://test-7cee72c0e768/Jax/sv-manta/PGP1-C2/PGP1-C2.svLocusGraphStats.tsv     --logging gs://test-7cee72c0e768/Jax/logs/Manta/PGP1-C2     --labels=pipe=manta,sample=manta,cohort=jax
