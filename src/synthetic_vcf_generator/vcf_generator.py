import sys
import uuid
import multiprocessing
from pathlib import Path


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
    if synthetic_vcf_path.suffix == ".gz":
        try:
            from Bio import bgzf as compressor
        except ImportError:  # pragma: no cover
            import gzip as compressor

        with compressor.open(synthetic_vcf_path, "wt") as gz_file, virtual_vcf as v_vcf:
            # for line in tqdm.tqdm(v_vcf, total=num_rows + 1):
            for line in v_vcf:
                gz_file.write(line)
    else:
        with (
            open(synthetic_vcf_path, "w", encoding="utf-8") as txt_file,
            virtual_vcf as v_vcf,
        ):
            # for line in tqdm.tqdm(v_vcf, total=num_rows + 1):
            for line in v_vcf:
                txt_file.write(line)


def synthetic_vcf_data(
    synthetic_vcf_path,
    num_rows,
    num_samples,
    chromosomes,
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
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (List[str]): List of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Phased or unphased genotypes.
        large_format (bool): Use large format VCF.
        reference_dir_path (Path or None): Path to imported reference data.
    """
    # Create and start the profiler
    # profiler = cProfile.Profile()
    # profiler.enable()

    virtual_vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
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

    total_num_rows = num_rows * len(chromosomes)
    to_vcf_file(
        virtual_vcf=virtual_vcf,
        synthetic_vcf_path=synthetic_vcf_path,
        num_rows=total_num_rows,
    )

    # Disable the profiler and print stats
    # profiler.disable()

    # Format and display the results
    # s = io.StringIO()
    # stats = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    # stats.print_stats(20)  # Print top 20 functions
    # print(s.getvalue())


def batch_synthetic_vcf_data(
    synthetic_vcf_dir,
    num_vcfs,
    vcf_prefix,
    num_rows,
    num_samples,
    chromosomes,
    seed,
    sample_prefix,
    id_type,
    phased,
    large_format,
    reference_dir_path,
    num_threads,
):
    """
    Generates synthetic VCF data and writes it to either a file or standard output.

    Args:
        synthetic_vcf_dir (Path): Path to the directory where output synthetic VCF files should be saved.
        num_vcfs (int): Number of VCF files to generate.
        vcf_prefix (str): Prefix added to the filename of every generated VCF.
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (List[str]): List of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Phased or unphased genotypes.
        large_format (bool): Use large format VCF.
        reference_dir_path (Path or None): Path to imported reference data.
        num_threads (int): Number of threads.
    """

    virtual_vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
        sample_prefix=sample_prefix,
        id_type=id_type,
        random_seed=seed,
        phased=phased,
        large_format=large_format,
        reference_dir=reference_dir_path,
    )

    total_num_rows = num_rows * len(chromosomes)
    results = []
    with multiprocessing.Pool(num_threads) as pool:
        for i in range(num_vcfs):
            filename = f"{vcf_prefix}{uuid.uuid4()}.vcf"
            kwargs = {
                "virtual_vcf": virtual_vcf,
                "synthetic_vcf_path": synthetic_vcf_dir / filename,
                "num_rows": total_num_rows,
            }
            result = pool.apply_async(to_vcf_file, (), kwargs)
            results.append(result)
        [r.wait() for r in results]
