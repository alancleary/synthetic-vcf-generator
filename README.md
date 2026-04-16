# synthetic-vcf-generator

synthetic-vcf-generator generates synthetic Variant Call Format (VCF) files for testing purposes.
This project is a hard fork of [fake-vcf](https://github.com/endast/fake-vcf) by [Magnus Wahlberg ](https://github.com/endast) (MIT license).
The purpose of this fork is to provide a much more extensive set of features with run-time performance sufficient for generating large collections of VCF files at scale.

## Usage

### Setup

To isntall synthetic-vcf-generator:
```shell
git clone https://github.com/alancleary/synthetic-vcf-generator.git
cd synthetic-vcf-generator
pip install .
```

### Run

synthetic-vcf-generator is under active development.
Use the `--help` flag to get up-to-date usage documentation:
```shell
synthetic-vcf-generator --help
```

synthetic-vcf-generator currently supports 3 commands: `generate`, `generate-batch`, and `import-reference`.

### `generate`

This command generates a single VCF file and writes it to the standard output or a user-defined file.
By default it generates random values for the REF column, but it can use a reference genome imported with the `import-reference` command to use real REF values instead.
For more details run:
```shell
synthetic-vcf-generator generate --help
```

#### Variant type distribution

The generator emits a mix of SNPs, MNPs, indels (INS/DEL), and SVs (DEL/INS/DUP/INV) controlled by three flags.
Each flag accepts a CSV of `key=integer` pairs whose values must sum to exactly 100.
Unspecified keys within a flag default to 0, and unknown keys are rejected.

| Flag | Default | Valid keys |
|---|---|---|
| `--type-weights` | `snp=80,mnp=5,indel=10,sv=5` | `snp`, `mnp`, `indel`, `sv` |
| `--indel-weights` | `ins=50,del=50` | `ins`, `del` |
| `--sv-weights` | `del=40,ins=20,dup=20,inv=20` | `del`, `ins`, `dup`, `inv` |

SV-specific header lines (`##ALT=<ID=...>`, `##INFO=<ID=SVTYPE,...>`, `##INFO=<ID=END,...>`, `##INFO=<ID=SVLEN,...>`) are emitted only when the `sv` weight is greater than 0.

Examples:
```shell
# SNP-only output
synthetic-vcf-generator generate --type-weights snp=100

# Skew the default mix
synthetic-vcf-generator generate --type-weights snp=50,mnp=10,indel=20,sv=20

# 100% insertions (ignores del sub-weight)
synthetic-vcf-generator generate --type-weights indel=100 --indel-weights ins=100

# 100% inversions
synthetic-vcf-generator generate --type-weights sv=100 --sv-weights inv=100
```

Size ranges are fixed in this version:
- MNP: 2–5 bp
- Indel (small INS/DEL): 1–50 bp
- SV (DEL/INS/DUP/INV): 50–10,000 bp

### `generate-batch`

This command generates a collection of VCF files and writes them to a user-defined directory.
The filename of each file generated will be a UUID with an optional user-defined prefix.
By default it generates random values for the REF column, but it can use a reference genome imported with the `import-reference` command to use real REF values instead.
The `--type-weights`, `--indel-weights`, and `--sv-weights` flags described under `generate` are accepted by this command as well.
For more details run:
```shell
synthetic-vcf-generator generate-batch --help
```

### `import-reference`

This command imports a genome from a FASTA file to use as a reference for the `generate` and `generate-batch` commands.
When the `generate` and `generate-batch` are given an imported genome, they will reference this genome when generating REF values.
For more details run:
```shell
synthetic-vcf-generator import-reference --help
```

### Examples

Scripts containing examples of common use cases are located in this repository's `./examples` directory.

## Development

### Dependencies

To install all dependencies, including development and testing, run:
```shell
pip install . --all-extras
```

To install only development dependencies, run:
```shell
pip install . --group dev
```

To install only test dependencies, run:
```shell
pip install . --group test
```

### Testing

Tests can be run with:
```shell
pytest
```

#### Code Coverage

[Coverage.py](https://coverage.readthedocs.io/) is included as a dev dependency.
To measure code coverage while running the test suite:
```shell
coverage run -m pytest
```

After the run completes, view a summary report in the terminal:
```shell
coverage report
```

Or generate a browsable HTML report in `./htmlcov`:
```shell
coverage html
```

### Linting and Formatting

Linting and formatting are automatically applied via [pre-commit](https://pre-commit.com/) hooks.
These hooks can be run manually using:
```shell
pre-commit run --all-files
```
