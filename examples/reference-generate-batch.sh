# Uses 10 threads to generate 100 single sample VCF files using the T2T Human genome for REF values
synthetic-vcf-generator generate-batch \
  --num-samples 1 \
  --id-type uuid \
  --num-rows 210000 \
  --reference-dir-path ./T2T/imported/ \
  --chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  --sample-prefix T2T_ \
  --output-directory generate_batch_output \
  --vcf-prefix T2T_ \
  --output-type bgzip \
  --num-vcfs 100 \
  --num-threads 10
