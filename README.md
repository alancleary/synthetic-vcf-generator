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

### `generate-batch`

This command generates a collection of VCF files and writes them to a user-defined directory.
The filename of each file generated will be a UUID with an optional user-defined prefix.
By default it generates random values for the REF column, but it can use a reference genome imported with the `import-reference` command to use real REF values instead.
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

### Linting and Formatting

Linting and formatting are automatically applied via [pre-commit](https://pre-commit.com/) hooks.
These hooks can be run manually using:
```shell
pre-commit run --all-files
```
