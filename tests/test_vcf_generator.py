import importlib.util
import io

import pytest

from synthetic_vcf_generator.vcf_generator import (
    batch_synthetic_vcf_data,
    synthetic_vcf_data,
    write_file,
    write_fileobj,
    write_stdout,
)
from synthetic_vcf_generator.virtual_vcf import VirtualVCF


BIO_AVAILABLE = importlib.util.find_spec("Bio") is not None


def make_virtual_vcf() -> VirtualVCF:
    return VirtualVCF(
        num_rows=2,
        num_samples=2,
        chromosomes=["chr1"],
        random_seed=42,
    )


@pytest.mark.generate_vcf
def test_write_fileobj_vcf():
    buf = io.BytesIO()
    write_fileobj(buf, make_virtual_vcf(), "vcf")
    data = buf.getvalue()
    assert data.startswith(b"##fileformat=VCFv4.2")


@pytest.mark.generate_vcf
def test_write_fileobj_gzip():
    buf = io.BytesIO()
    write_fileobj(buf, make_virtual_vcf(), "gzip")
    data = buf.getvalue()
    assert data[:2] == b"\x1f\x8b"


@pytest.mark.generate_vcf
@pytest.mark.skipif(not BIO_AVAILABLE, reason="Biopython not installed")
def test_write_fileobj_bgzip(tmp_path):
    out = tmp_path / "out.vcf.gz"
    with open(out, "wb") as f:
        write_fileobj(f, make_virtual_vcf(), "bgzip")
    data = out.read_bytes()
    assert data[:3] == b"\x1f\x8b\x08"


@pytest.mark.generate_vcf
def test_write_fileobj_invalid_type():
    buf = io.BytesIO()
    with pytest.raises(ValueError, match="not supported"):
        write_fileobj(buf, make_virtual_vcf(), "bogus")


@pytest.mark.generate_vcf
def test_write_stdout(capfdbinary):
    write_stdout(make_virtual_vcf(), "vcf")
    captured = capfdbinary.readouterr()
    assert captured.out.startswith(b"##fileformat=VCFv4.2")


@pytest.mark.generate_vcf
def test_write_file(tmp_path):
    out = tmp_path / "out.vcf"
    write_file(make_virtual_vcf(), "vcf", out)
    assert out.exists()
    assert out.read_bytes().startswith(b"##fileformat=VCFv4.2")


@pytest.mark.generate_vcf
def test_synthetic_vcf_data_to_file(tmp_path):
    out = tmp_path / "out.vcf"
    synthetic_vcf_data(
        synthetic_vcf_path=out,
        output_type="vcf",
        num_rows=2,
        num_samples=2,
        chromosomes=["chr1"],
        seed=42,
        sample_prefix="sample_",
        id_type="padded_count",
        phased=True,
        large_format=True,
        reference_dir_path=None,
    )
    assert out.exists()
    assert out.read_bytes().startswith(b"##fileformat=VCFv4.2")


@pytest.mark.generate_vcf
def test_synthetic_vcf_data_to_stdout(capfdbinary):
    synthetic_vcf_data(
        synthetic_vcf_path=None,
        output_type="vcf",
        num_rows=2,
        num_samples=2,
        chromosomes=["chr1"],
        seed=42,
        sample_prefix="sample_",
        id_type="padded_count",
        phased=True,
        large_format=True,
        reference_dir_path=None,
    )
    captured = capfdbinary.readouterr()
    assert captured.out.startswith(b"##fileformat=VCFv4.2")


def _batch_kwargs(**overrides):
    kwargs = dict(
        output_type="vcf",
        num_vcfs=1,
        vcf_prefix="",
        num_rows=2,
        num_samples=2,
        chromosomes=["chr1"],
        seed=42,
        sample_prefix="sample_",
        id_type="padded_count",
        phased=True,
        large_format=True,
        reference_dir_path=None,
        num_threads=1,
    )
    kwargs.update(overrides)
    return kwargs


@pytest.mark.generate_batch
def test_batch_synthetic_vcf_data_creates_dir(tmp_path):
    out_dir = tmp_path / "nested" / "batch"
    assert not out_dir.exists()
    batch_synthetic_vcf_data(synthetic_vcf_dir=out_dir, **_batch_kwargs())
    assert out_dir.is_dir()
    assert len(list(out_dir.glob("*.vcf"))) == 1


@pytest.mark.generate_batch
def test_batch_synthetic_vcf_data_file_count(tmp_path):
    out_dir = tmp_path / "batch"
    batch_synthetic_vcf_data(
        synthetic_vcf_dir=out_dir,
        **_batch_kwargs(num_vcfs=3, vcf_prefix="pfx_"),
    )
    files = sorted(out_dir.glob("*.vcf"))
    assert len(files) == 3
    for f in files:
        assert f.name.startswith("pfx_")


@pytest.mark.generate_batch
def test_batch_synthetic_vcf_data_gzip_ext(tmp_path):
    out_dir = tmp_path / "batch"
    batch_synthetic_vcf_data(
        synthetic_vcf_dir=out_dir,
        **_batch_kwargs(output_type="gzip", num_vcfs=2),
    )
    files = sorted(out_dir.glob("*.vcf.gz"))
    assert len(files) == 2
    for f in files:
        assert f.read_bytes()[:2] == b"\x1f\x8b"
