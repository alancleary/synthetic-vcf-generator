# Examples

This directory contains examples of how to use synthetic-vcf-generator.
Specifically, each of the following scripts contains a different use case:
```shell
./examples
├── README.md
├── generate.sh
├── generate-batch.sh
├── import-reference.sh
├── reference-generate.sh
└── reference-generate-batch.sh
```

### `generate.sh`

This script generates a single sample VCF file for a human (`./human.vcf`) with random REF values.
This is useful if you need to generate a single VCF file that will be used by itself.
If you need real REF values, then use the `import-reference.sh` and `reference-generate.sh` scripts.
Note that omitting the `--synthetic-vcf-path` flag will cause the VCF file to written to the standard output.

### `generate-batch.sh`

This script generates 100 single sample VCF files for human with random REF values.
The filename of each file generated will be a UUID optionally prefixed by the `--vcf-prefix` flag's argument.
This is useful if you need to generate a collection of VCF files and consistent REF values is **not** necessary.
If you need real REF values, then use of the `import-reference.sh` and `reference-batch.sh` scripts.

### `import-reference.sh`

This script imports a genome from a FASTA file to use as a reference when generating REF values.
This script specifically assumes that a copy of the T2T Human genome version CHM13 is located at `./T2T/CHM13.fa`.
A copy can be download from the Human Pangenomics AWS S3 bucket: [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/).

### `reference-generate.sh`

This script generates a single sample VCF file for a human (`./T2T.vcf`) using the T2T Human genome imported by `import-reference.sh` for REF values.
This is useful if you need to generate a single VCF file but the REF values need to be accurate.
As with `generate.sh`, omitting the `--synthetic-vcf-path` flag will cause the VCF file to written to the standard output.

### `reference-generate-batch.sh`

This script generates 100 single sample VCF files for human using the T2T Human genome imported by `import-reference.sh` for REF values.
As with `generate-batch.sh`, the filename of each file generated will be a UUID optionally prefixed by the `--vcf-prefix` flag's argument.
This is useful if you need to generate a collection of VCF files and consistent REF values **is** necessary.
