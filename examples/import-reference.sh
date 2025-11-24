# NOTE: This script assumes a copy of the T2T Human genome version CHM13 is located at ./T2T/CHM13.fa
# A copy can be download from the Human Pangenomics AWS S3 bucket:
# https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/

# Imports the T2T Human genome to use as a reference with the "generate" command
synthetic-vcf-generator import-reference \
  ./T2T/CHM13.fa \
  ./T2T/imported/ \
  --included-chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
