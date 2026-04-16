import importlib.util
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from synthetic_vcf_generator import version
from synthetic_vcf_generator.__main__ import app
from tests.test_virtual_vcf import NUMBER_NON_SAMPLE_COL

runner = CliRunner()

GENERATE_CMD = "generate"
GENERATE_BATCH_CMD = "generate-batch"
IMPORT_REFERENCE_CMD = "import-reference"
test_data_dir = Path(__file__).resolve().parent / "test_data"
reference_dir = test_data_dir / "reference/seq"


def is_gz_file(filepath):
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


def is_bgzip_compressed(file_path):
    with open(file_path, "rb") as file:
        # Read the first three bytes from the file
        magic_bytes = file.read(3)

    # Check if the magic bytes indicate bgzip compression
    return magic_bytes == b"\x1f\x8b\x08"


script_dir = Path(__file__).resolve().parent


@pytest.mark.generate_vcf
def test_cli_reference_import_no_input():
    result = runner.invoke(app, [IMPORT_REFERENCE_CMD])
    assert result.exit_code == 2


@pytest.mark.reference_import
def test_cli_reference_import_only_input():
    result = runner.invoke(
        app,
        [
            IMPORT_REFERENCE_CMD,
            "reference.fa",
        ],
    )
    assert result.exit_code == 2


@pytest.mark.reference_import
def test_cli_reference_import_input_output(tmp_path):
    small_reference_path = script_dir / (
        "../tests/test_data/reference/reference_small.fa"
    )

    result = runner.invoke(
        app,
        [IMPORT_REFERENCE_CMD, small_reference_path.as_posix(), tmp_path.as_posix()],
    )
    assert result.exit_code == 0


@pytest.mark.reference_import
def test_cli_reference_import_input_output_incl_chrom(tmp_path):
    small_reference_path = script_dir / (
        "../tests/test_data/reference/reference_small.fa"
    )

    result = runner.invoke(
        app,
        [
            IMPORT_REFERENCE_CMD,
            "--included-chromosomes",
            "chr1",
            small_reference_path.as_posix(),
            tmp_path.as_posix(),
        ],
    )
    assert result.exit_code == 0


@pytest.mark.generate_vcf
def test_cli_generate_no_input():
    result = runner.invoke(app, [GENERATE_CMD])
    assert result.exit_code == 0
    assert "source=VirtualVCF" in result.stdout


@pytest.mark.generate_vcf
def test_cli_generate_no_compression_output(tmp_path):
    output_file = tmp_path / "example.vcf"
    result = runner.invoke(app, [GENERATE_CMD, "-o", output_file])
    assert result.exit_code == 0
    assert output_file.exists()

    try:
        assert not is_bgzip_compressed(output_file)
    except AssertionError:
        assert not is_gz_file(output_file)


@pytest.mark.generate_vcf
def test_cli_generation_compression(tmp_path):
    output_file = tmp_path / "example.vcf.gz"
    result = runner.invoke(app, [GENERATE_CMD, "-o", output_file])
    assert result.exit_code == 0
    assert output_file.exists()

    # If biopython is installed check that we wrote a bgzip file
    try:
        assert is_bgzip_compressed(output_file)
    except AssertionError:
        assert is_gz_file(output_file)


@pytest.mark.generate_vcf
def test_cli_generation_compression_no_bgzip(tmp_path):
    output_file = tmp_path / "example.vcf.gz"
    result = runner.invoke(app, [GENERATE_CMD, "-o", output_file])
    assert result.exit_code == 0
    assert output_file.exists()

    # If biopython is installed check that we wrote a bgzip file
    try:
        assert is_bgzip_compressed(output_file)
    except AssertionError:
        assert is_gz_file(output_file)


@pytest.mark.generate_vcf
def test_cli_generation_seed_same(tmp_path):
    args = [GENERATE_CMD, "--seed", "42"]

    result_1 = runner.invoke(app, args)
    result_2 = runner.invoke(app, args)
    assert result_1.exit_code == 0
    assert result_2.exit_code == 0
    assert result_1.stdout == result_2.stdout


@pytest.mark.generate_vcf
def test_cli_generation_seed_differ(tmp_path):
    base_args = [GENERATE_CMD, "--seed"]
    result_1 = runner.invoke(app, base_args + ["42"])
    result_2 = runner.invoke(app, base_args + ["1337"])
    assert result_1.exit_code == 0
    assert result_2.exit_code == 0
    assert result_1.stdout != result_2.stdout


@pytest.mark.generate_vcf
def test_cli_generation_version(tmp_path):
    result = runner.invoke(app, [GENERATE_CMD, "--version"])
    assert result.exit_code == 0
    assert version in result.stdout


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("chr",),
    [
        *[(f"chr{c}",) for c in range(1, 22)],
    ],
)
def test_cli_generation_chr_flag(chr):
    result = runner.invoke(app, [GENERATE_CMD, "-c", chr])
    assert result.exit_code == 0
    assert chr in result.stdout


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("prefix",),
    [
        ("BAM",),
        ("SAM",),
        ("MAM",),
        ("wham",),
        ("tapir",),
    ],
)
def test_cli_generation_sample_prefix_flag(prefix):
    result = runner.invoke(app, [GENERATE_CMD, "-p", prefix])
    assert result.exit_code == 0
    assert f"{prefix}0000" in result.stdout


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("expected_rows",),
    [
        *[(r,) for r in range(1, 100, 10)],
    ],
)
def test_cli_generation_number_rows(expected_rows):
    result = runner.invoke(app, [GENERATE_CMD, "-r", f"{expected_rows}"])
    row_count = len([r for r in result.stdout.split("\n") if r.startswith("chr1")])
    assert result.exit_code == 0
    assert row_count == expected_rows


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("expected_rows",),
    [
        *[(r,) for r in range(1, 10, 1)],
    ],
)
def test_cli_generation_number_rows_with_reference(expected_rows):
    result = runner.invoke(
        app, [GENERATE_CMD, "-r", f"{expected_rows}", "-f", f"{reference_dir}"]
    )
    row_count = len([r for r in result.stdout.split("\n") if r.startswith("chr1")])
    assert result.exit_code == 0
    assert row_count == expected_rows


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("expected_sample_count",),
    [
        *[(r,) for r in range(1, 100, 10)],
    ],
)
def test_cli_generation_number_samples(expected_sample_count):
    result = runner.invoke(app, [GENERATE_CMD, "-s", f"{expected_sample_count}"])
    sample_count = len(
        [r for r in result.stdout.split("\n") if r.startswith("chr1")][0].split("\t")
    )
    assert result.exit_code == 0
    assert sample_count == expected_sample_count + NUMBER_NON_SAMPLE_COL


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("expected_sample_count",),
    [
        *[(r,) for r in range(1, 100, 10)],
    ],
)
def test_cli_generation_number_samples_with_reference(expected_sample_count):
    result = runner.invoke(
        app, [GENERATE_CMD, "-s", f"{expected_sample_count}", "-f", f"{reference_dir}"]
    )
    sample_count = len(
        [r for r in result.stdout.split("\n") if r.startswith("chr1")][0].split("\t")
    )
    assert result.exit_code == 0
    assert sample_count == expected_sample_count + NUMBER_NON_SAMPLE_COL


@pytest.mark.generate_vcf
def test_cli_generation_gzip_explicit(tmp_path):
    import gzip

    output_file = tmp_path / "example.vcf.gz"
    result = runner.invoke(
        app, [GENERATE_CMD, "-o", output_file.as_posix(), "-e", "gzip"]
    )
    assert result.exit_code == 0
    assert output_file.exists()
    assert is_gz_file(output_file)
    with gzip.open(output_file, "rb") as f:
        assert f.read().startswith(b"##fileformat=VCFv4.2")


@pytest.mark.generate_vcf
def test_cli_generation_id_type_count():
    result = runner.invoke(app, [GENERATE_CMD, "-i", "count", "-s", "3"])
    assert result.exit_code == 0
    header_line = [r for r in result.stdout.split("\n") if r.startswith("#CHROM")][0]
    samples = header_line.split("\t")[NUMBER_NON_SAMPLE_COL:]
    assert samples == ["sample_1", "sample_2", "sample_3"]


@pytest.mark.generate_vcf
def test_cli_generation_id_type_uuid():
    result = runner.invoke(app, [GENERATE_CMD, "-i", "uuid", "-s", "2"])
    assert result.exit_code == 0
    header_line = [r for r in result.stdout.split("\n") if r.startswith("#CHROM")][0]
    samples = header_line.split("\t")[NUMBER_NON_SAMPLE_COL:]
    uuid_re = re.compile(
        r"^sample_[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
    )
    assert len(samples) == 2
    for name in samples:
        assert uuid_re.match(name)


@pytest.mark.generate_vcf
def test_cli_generation_invalid_extension(tmp_path):
    output_file = tmp_path / "example.txt"
    result = runner.invoke(app, [GENERATE_CMD, "-o", output_file.as_posix()])
    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)
    assert 'extension ".vcf" or ".gz"' in str(result.exception)


@pytest.mark.generate_vcf
def test_cli_generation_gz_with_vcf_type(tmp_path):
    output_file = tmp_path / "example.vcf.gz"
    result = runner.invoke(
        app, [GENERATE_CMD, "-o", output_file.as_posix(), "-e", "vcf"]
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)
    assert 'cannot have output type "vcf"' in str(result.exception)


@pytest.mark.generate_vcf
def test_cli_generation_gz_bio_missing(tmp_path, monkeypatch):
    original_find_spec = importlib.util.find_spec

    def fake_find_spec(name, *args, **kwargs):
        if name == "Bio":
            return None
        return original_find_spec(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", fake_find_spec)

    output_file = tmp_path / "example.vcf.gz"
    result = runner.invoke(app, [GENERATE_CMD, "-o", output_file.as_posix()])
    assert result.exit_code == 0
    assert output_file.exists()
    assert is_gz_file(output_file)


@pytest.mark.generate_vcf
def test_cli_generation_bgzip_bio_missing(monkeypatch):
    original_find_spec = importlib.util.find_spec

    def fake_find_spec(name, *args, **kwargs):
        if name == "Bio":
            return None
        return original_find_spec(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", fake_find_spec)

    result = runner.invoke(app, [GENERATE_CMD, "-e", "bgzip"])
    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)
    assert "requires Biopython" in str(result.exception)


@pytest.mark.generate_batch
def test_cli_generate_batch_default(tmp_path):
    out_dir = tmp_path / "batch"
    result = runner.invoke(
        app, [GENERATE_BATCH_CMD, "-d", out_dir.as_posix(), "-t", "1"]
    )
    assert result.exit_code == 0
    files = sorted(out_dir.glob("*.vcf"))
    assert len(files) == 10


@pytest.mark.generate_batch
def test_cli_generate_batch_num_vcfs(tmp_path):
    out_dir = tmp_path / "batch"
    result = runner.invoke(
        app,
        [GENERATE_BATCH_CMD, "-d", out_dir.as_posix(), "-n", "3", "-t", "1"],
    )
    assert result.exit_code == 0
    files = sorted(out_dir.glob("*.vcf"))
    assert len(files) == 3


@pytest.mark.generate_batch
def test_cli_generate_batch_prefix(tmp_path):
    out_dir = tmp_path / "batch"
    result = runner.invoke(
        app,
        [
            GENERATE_BATCH_CMD,
            "-d",
            out_dir.as_posix(),
            "-n",
            "2",
            "-v",
            "GRCh38_",
            "-t",
            "1",
        ],
    )
    assert result.exit_code == 0
    files = sorted(out_dir.glob("*.vcf"))
    assert len(files) == 2
    for f in files:
        assert f.name.startswith("GRCh38_")


@pytest.mark.generate_batch
def test_cli_generate_batch_creates_dir(tmp_path):
    out_dir = tmp_path / "nested" / "batch"
    assert not out_dir.exists()
    result = runner.invoke(
        app,
        [GENERATE_BATCH_CMD, "-d", out_dir.as_posix(), "-n", "1", "-t", "1"],
    )
    assert result.exit_code == 0
    assert out_dir.is_dir()
    assert len(list(out_dir.glob("*.vcf"))) == 1


@pytest.mark.generate_batch
def test_cli_generate_batch_compression(tmp_path):
    out_dir = tmp_path / "batch"
    result = runner.invoke(
        app,
        [
            GENERATE_BATCH_CMD,
            "-d",
            out_dir.as_posix(),
            "-n",
            "2",
            "-e",
            "gzip",
            "-t",
            "1",
        ],
    )
    assert result.exit_code == 0
    files = sorted(out_dir.glob("*.vcf.gz"))
    assert len(files) == 2
    for f in files:
        assert is_gz_file(f)


@pytest.mark.generate_batch
def test_cli_generate_batch_version():
    result = runner.invoke(app, [GENERATE_BATCH_CMD, "--version"])
    assert result.exit_code == 0
    assert version in result.stdout


@pytest.mark.generate_batch
def test_cli_generate_batch_seed(tmp_path):
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    args_common = [
        GENERATE_BATCH_CMD,
        "-n",
        "1",
        "-v",
        "run_",
        "--seed",
        "42",
        "-t",
        "1",
    ]
    r1 = runner.invoke(app, args_common + ["-d", out_a.as_posix()])
    r2 = runner.invoke(app, args_common + ["-d", out_b.as_posix()])
    assert r1.exit_code == 0
    assert r2.exit_code == 0

    def read_body(d: Path) -> str:
        [f] = list(d.glob("*.vcf"))
        return "\n".join(
            line for line in f.read_text().splitlines() if not line.startswith("##")
        )

    assert read_body(out_a) == read_body(out_b)
