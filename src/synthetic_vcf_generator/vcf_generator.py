import sys
from pathlib import Path

import tqdm

from synthetic_vcf_generator.virtual_vcf import VirtualVCF


def to_std_out(virtual_vcf: VirtualVCF) -> None:
    """
    Writes VirtualVCF data to standard output.

    Args:
        virtual_vcf (VirtualVCF): VirtualVCF object containing the data.
    """
    with virtual_vcf as v_vcf:
        for line in v_vcf:
            sys.stdout.write(line)


def to_vcf_file(
    virtual_vcf: VirtualVCF, synthetic_vcf_path: Path, num_rows: int
) -> None:
    """
    Writes VirtualVCF data to a VCF file.

    Args:
        virtual_vcf (VirtualVCF): VirtualVCF object containing the data.
        synthetic_vcf_path (Path): Path to the synthetic VCF file.
        num_rows (int): Number of rows.
    """
    print(f"Writing to file {synthetic_vcf_path}")

    if synthetic_vcf_path.suffix == ".gz":
        print("(Using compression)")
        try:
            from Bio import bgzf as compressor
        except ImportError:  # pragma: no cover
            print("Biopython not installed, falling back to gzip instead of bgzip")
            import gzip as compressor

        with compressor.open(synthetic_vcf_path, "wt") as gz_file, virtual_vcf as v_vcf:
            for line in tqdm.tqdm(v_vcf, total=num_rows + 1):
                gz_file.write(line)
    else:
        print("(No compression)")
        with (
            open(synthetic_vcf_path, "w", encoding="utf-8") as txt_file,
            virtual_vcf as v_vcf,
        ):
            for line in tqdm.tqdm(v_vcf, total=num_rows + 1):
                txt_file.write(line)

    print(f"Done, data written to {synthetic_vcf_path}")


def synthetic_vcf_data(
    synthetic_vcf_path,
    num_rows,
    num_samples,
    chromosome,
    seed,
    sample_prefix,
    id_type,
    phased,
    large_format,
    reference_dir_path,
):
    """
    Generates synthetic VCF data and writes it to either a file or standard output.

    Args:
        synthetic_vcf_path (Path or None): Path to the synthetic VCF file or None to write to standard output.
        num_rows (int): Number of rows.
        num_samples (int): Number of samples.
        chromosome (str): Chromosome identifier.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Phased or unphased genotypes.
        large_format (bool): Use large format VCF.
        reference_dir_path (Path or None): Path to imported reference data.
    """
    virtual_vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=num_samples,
        chromosome=chromosome,
        sample_prefix=sample_prefix,
        id_type=id_type,
        random_seed=seed,
        phased=phased,
        large_format=large_format,
        reference_dir=reference_dir_path,
    )

    if synthetic_vcf_path is None:
        to_std_out(virtual_vcf=virtual_vcf)
        return

    to_vcf_file(
        virtual_vcf=virtual_vcf,
        synthetic_vcf_path=synthetic_vcf_path,
        num_rows=num_rows,
    )
