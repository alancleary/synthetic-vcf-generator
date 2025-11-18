# synthetic-vcf-generator

synthetic-vcf-generator generates synthetic Variant Call Format (VCF) files for testing purposes.
This project is a hard fork of [fake-vcf](https://github.com/endast/fake-vcf) by [Magnus Wahlberg ](https://github.com/endast) (MIT license).
The purpose of this fork is to provide a much more extensive set of features with run-time performance sufficient for generating large collections of VCF files at scale.

## Usage

### Setup

Currently synthetic-vcf-generator uses [PDM](https://pdm-project.org/en/latest/) for project and dependency management:
```shell
git clone https://github.com/alancleary/synthetic-vcf-generator.git
cd synthetic-vcf-generator
pdm install
pdm build
```

If you want to write bgzip files instead gzip when writing compressed gzip files then install the optional bgzip dependencies:
```shell
pdm install -Gbgzip
```

### Run

synthetic-vcf-generator is under active development.
Use the `--help` flag to get up-to-date usage documentation:
```shell
pdm run synthetic-vcf-generator --help
```

## Development

### Dependencies

To install all dependencies, including development and testing, run:
```shell
pdm install -G:all
```

To install only development dependencies, run:
```shell
pdm install -Gdev
```

To install only test dependencies, run:
```shell
pdm install -Gtest
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
