from typing import List, Literal

import os
import time
from pathlib import Path

import typer
from rich.console import Console

from synthetic_vcf_generator import version
from synthetic_vcf_generator.vcf_generator import (
    batch_synthetic_vcf_data,
    synthetic_vcf_data,
)
from synthetic_vcf_generator.vcf_reference import import_reference

app = typer.Typer(
    name="synthetic-vcf-generator",
    help="A synthetic VCF file generator",
    add_completion=False,
)
console = Console()


def version_callback(print_version: bool) -> None:
    """
    Callback function to print the version of the package.

    Args:
        print_version (bool): Flag to print the version.

    Raises:
        Exit: If the print_version flag is set.
    """
    if print_version:
        console.print(
            f"[yellow]synthetic-vcf-generator[/] version: [bold blue]{version}[/]"
        )
        raise typer.Exit()


@app.command(name="import-reference")
def vcf_reference_import(
    reference_file_path: Path = typer.Argument(
        help="Path to reference fasta file.",
    ),
    reference_storage_path: Path = typer.Argument(
        help="Where to store the references.",
    ),
    included_chromosomes: List[str] = typer.Option(
        None,
        "--included_chromosomes",
        "-c",
        help="List of chromosomes to extract from reference, if not specified all will be imported",
    ),
) -> None:
    """
    Import reference fasta file and extract specified chromosomes if provided.

    Parameters:
        reference_file_path (Path): Path to reference fasta file.
        reference_storage_path (Path): Where to store the references.
        included_chromosomes (Optional[List[str]], optional): List of chromosomes
            to extract from reference. If not specified, all will be imported.

    Example:
        To import a reference file and extract specific chromosomes:
        ```
        vcf_reference_import("path/to/reference.fasta", "output/directory", included_chromosomes=["chr1", "chr2"])
        ```

        To import a reference file without extracting specific chromosomes:
        ```
        vcf_reference_import("path/to/reference.fasta", "output/directory")
        ```
    """

    print(f"Importing reference {reference_file_path}")
    if included_chromosomes:
        print(
            f"Importing {len(included_chromosomes)} chromosomes from reference {reference_file_path}: {', '.join(included_chromosomes)}"
        )
    else:
        print(f"Importing all chromosomes from reference {reference_file_path}")

    start_time = time.time()
    import_reference(
        file_path=reference_file_path,
        output_dir=reference_storage_path,
        include_sequences=included_chromosomes,
    )
    end_time = time.time()

    print(f"Reference imported in {end_time - start_time:.2f} seconds.")


@app.command(name="generate")
def main(
    synthetic_vcf_path: Path = typer.Option(
        None,
        "--synthetic_vcf_path",
        "-o",
        help="Path to synthetic vcf file. If the path ends with .gz the file will be gzipped.",
    ),
    num_rows: int = typer.Option(
        10,
        "--num_rows",
        "-r",
        help="Number of rows (variants) to generate per chromosome",
    ),
    num_samples: int = typer.Option(
        10, "--num_samples", "-s", help="Number of sample to generate."
    ),
    chromosomes: str = typer.Option(
        "chr1",
        "--chromosomes",
        "-c",
        help="CSV list of chromosome IDs to include in the VCF",
    ),
    seed: int = typer.Option(None, "--seed", help="Random seed to use, default none."),
    sample_prefix: str = typer.Option(
        "sample_",
        "--sample_prefix",
        "-p",
        help="Sample prefix, e.g. SAM => SAM0000001 SAM0000002",
    ),
    id_type: Literal["count", "padded_count", "uuid"] = typer.Option(
        "padded_count",
        "--id_type",
        "-i",
        help="Type of unique ID to use for samples",
    ),
    phased: bool = typer.Option(default=True, help="Simulate phased"),
    large_format: bool = typer.Option(default=True, help="Write large format vcf"),
    print_version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the synthetic-vcf-generator package.",
    ),
    reference_dir: Path = typer.Option(
        None,
        "--reference-dir-path",
        "-f",
        help="Path to imported refernce directory.",
        exists=True,
    ),
) -> None:
    """
    Generate synthetic VCF data

    Args:
        synthetic_vcf_path (Path): Path to synthetic VCF file or None to write to standard output.
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (str): CSV list of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Simulate phased genotypes.
        large_format (bool): Write large format VCF.
        print_version (bool): Flag to print the version of the synthetic-vcf-generator package.
        reference_dir (Path): Path to directory containing imported reference_data.
    """
    chromosomes = chromosomes.split(",")
    synthetic_vcf_data(
        synthetic_vcf_path=synthetic_vcf_path,
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
        seed=seed,
        sample_prefix=sample_prefix,
        id_type=id_type,
        phased=phased,
        large_format=large_format,
        reference_dir_path=reference_dir,
    )


@app.command(name="generate-batch")
def generate_batch(
    synthetic_vcf_dir: Path = typer.Option(
        None,
        "--output-directory",
        "-d",
        help="Path to the directory where output synthetic VCF files should be saved.",
    ),
    num_vcfs: int = typer.Option(
        10, "--num-vcfs", "-n", help="Number of VCF files to generate."
    ),
    vcf_prefix: str = typer.Option(
        "",
        "-v",
        "--vcf_prefix",
        help="VCF file prefix, e.g.: GRCh38.p14_ => GRCh38.p14_293439a-5334-431d-b2ab-0f831463fada.vcf",
    ),
    num_rows: int = typer.Option(
        10,
        "--num_rows",
        "-r",
        help="Number of rows (variants) to generate per chromosome",
    ),
    num_samples: int = typer.Option(
        10, "--num_samples", "-s", help="Number of sample to generate."
    ),
    chromosomes: str = typer.Option(
        "chr1",
        "--chromosomes",
        "-c",
        help="CSV list of chromosome IDs to include in the VCF",
    ),
    seed: int = typer.Option(None, "--seed", help="Random seed to use, default none."),
    sample_prefix: str = typer.Option(
        "sample_",
        "--sample_prefix",
        "-p",
        help="Sample prefix, e.g. SAM => SAM0000001 SAM0000002",
    ),
    id_type: Literal["count", "padded_count", "uuid"] = typer.Option(
        "padded_count",
        "--id_type",
        "-i",
        help="Type of unique ID to use for samples",
    ),
    phased: bool = typer.Option(default=True, help="Simulate phased"),
    large_format: bool = typer.Option(default=True, help="Write large format vcf"),
    print_version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the synthetic-vcf-generator package.",
    ),
    reference_dir: Path = typer.Option(
        None,
        "--reference-dir-path",
        "-f",
        help="Path to imported refernce directory.",
        exists=True,
    ),
    num_threads: int = typer.Option(
        os.cpu_count(), "--num-threads", "-t", help="Number of threads."
    ),
) -> None:
    """
    Generate a batch of synthetic VCF data

    Args:
        synthetic_vcf_dir (Path): Path to the directory where output synthetic VCF files should be saved.
        num_vcfs (int): Number of VCF files to generate.
        vcf_prefix (str): Prefix added to the filename of every generated VCF.
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (str): CSV list of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Simulate phased genotypes.
        large_format (bool): Write large format VCF.
        print_version (bool): Flag to print the version of the synthetic-vcf-generator package.
        reference_dir (Path): Path to directory containing imported reference_data.
        num_threads (int): Number of threads.
    """
    chromosomes = chromosomes.split(",")
    batch_synthetic_vcf_data(
        synthetic_vcf_dir=synthetic_vcf_dir,
        num_vcfs=num_vcfs,
        vcf_prefix=vcf_prefix,
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
        seed=seed,
        sample_prefix=sample_prefix,
        id_type=id_type,
        phased=phased,
        large_format=large_format,
        reference_dir_path=reference_dir,
        num_threads=num_threads,
    )


if __name__ == "__main__":
    app()  # pragma: no cover
