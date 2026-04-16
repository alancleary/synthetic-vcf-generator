import random

import fastrand
import pytest

from synthetic_vcf_generator import variant_types


# --- parse_weights --------------------------------------------------------


@pytest.mark.generate_vcf
def test_parse_weights_default_when_none():
    result = variant_types.parse_weights(
        None,
        variant_types.VALID_TYPE_KEYS,
        variant_types.DEFAULT_TYPE_WEIGHTS,
    )
    assert result == variant_types.DEFAULT_TYPE_WEIGHTS
    # The returned dict is a copy — mutating it must not affect the default
    result["snp"] = 0
    assert variant_types.DEFAULT_TYPE_WEIGHTS["snp"] == 80


@pytest.mark.generate_vcf
def test_parse_weights_single_key_sums_to_100():
    result = variant_types.parse_weights(
        "snp=100",
        variant_types.VALID_TYPE_KEYS,
        variant_types.DEFAULT_TYPE_WEIGHTS,
    )
    assert result == {"snp": 100, "mnp": 0, "indel": 0, "sv": 0}


@pytest.mark.generate_vcf
def test_parse_weights_full_distribution():
    result = variant_types.parse_weights(
        "snp=70,mnp=10,indel=10,sv=10",
        variant_types.VALID_TYPE_KEYS,
        variant_types.DEFAULT_TYPE_WEIGHTS,
    )
    assert result == {"snp": 70, "mnp": 10, "indel": 10, "sv": 10}


@pytest.mark.generate_vcf
def test_parse_weights_zero_allowed():
    result = variant_types.parse_weights(
        "snp=50,mnp=50,indel=0,sv=0",
        variant_types.VALID_TYPE_KEYS,
        variant_types.DEFAULT_TYPE_WEIGHTS,
    )
    assert result == {"snp": 50, "mnp": 50, "indel": 0, "sv": 0}


@pytest.mark.generate_vcf
def test_parse_weights_tolerates_whitespace():
    result = variant_types.parse_weights(
        " snp = 60 , mnp = 40 ",
        variant_types.VALID_TYPE_KEYS,
        variant_types.DEFAULT_TYPE_WEIGHTS,
    )
    assert result["snp"] == 60
    assert result["mnp"] == 40


@pytest.mark.generate_vcf
def test_parse_weights_sum_less_than_100_errors():
    with pytest.raises(ValueError, match="sum to 100"):
        variant_types.parse_weights(
            "snp=80,mnp=5,indel=10",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_sum_greater_than_100_errors():
    with pytest.raises(ValueError, match="sum to 100"):
        variant_types.parse_weights(
            "snp=90,indel=15",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_rejects_float():
    with pytest.raises(ValueError, match="integer"):
        variant_types.parse_weights(
            "snp=50.5,mnp=49.5",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_rejects_non_numeric():
    with pytest.raises(ValueError, match="integer"):
        variant_types.parse_weights(
            "snp=abc",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_rejects_negative():
    with pytest.raises(ValueError, match="non-negative"):
        variant_types.parse_weights(
            "snp=-5,mnp=105",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_unknown_key_errors():
    with pytest.raises(ValueError, match="Unknown weight key"):
        variant_types.parse_weights(
            "foo=100",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_missing_equals_errors():
    with pytest.raises(ValueError, match="expected 'key=value'"):
        variant_types.parse_weights(
            "snp100",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_empty_pair_errors():
    with pytest.raises(ValueError, match="Empty"):
        variant_types.parse_weights(
            "snp=100,",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


@pytest.mark.generate_vcf
def test_parse_weights_duplicate_key_errors():
    with pytest.raises(ValueError, match="Duplicate"):
        variant_types.parse_weights(
            "snp=50,snp=50",
            variant_types.VALID_TYPE_KEYS,
            variant_types.DEFAULT_TYPE_WEIGHTS,
        )


# --- pick_variant_kind ----------------------------------------------------


@pytest.mark.generate_vcf
def test_pick_variant_kind_all_snp():
    rng = random.Random(42)
    for _ in range(100):
        assert (
            variant_types.pick_variant_kind(
                {"snp": 100, "mnp": 0, "indel": 0, "sv": 0},
                {"ins": 50, "del": 50},
                {"del": 40, "ins": 20, "dup": 20, "inv": 20},
                rng,
            )
            == "snp"
        )


@pytest.mark.generate_vcf
def test_pick_variant_kind_all_mnp():
    rng = random.Random(42)
    for _ in range(100):
        assert (
            variant_types.pick_variant_kind(
                {"snp": 0, "mnp": 100, "indel": 0, "sv": 0},
                {"ins": 50, "del": 50},
                {"del": 40, "ins": 20, "dup": 20, "inv": 20},
                rng,
            )
            == "mnp"
        )


@pytest.mark.generate_vcf
def test_pick_variant_kind_indel_routing():
    rng = random.Random(42)
    kinds = {
        variant_types.pick_variant_kind(
            {"snp": 0, "mnp": 0, "indel": 100, "sv": 0},
            {"ins": 50, "del": 50},
            {"del": 40, "ins": 20, "dup": 20, "inv": 20},
            rng,
        )
        for _ in range(200)
    }
    assert kinds == {"ins_small", "del_small"}


@pytest.mark.generate_vcf
def test_pick_variant_kind_sv_routing():
    rng = random.Random(42)
    kinds = {
        variant_types.pick_variant_kind(
            {"snp": 0, "mnp": 0, "indel": 0, "sv": 100},
            {"ins": 50, "del": 50},
            {"del": 25, "ins": 25, "dup": 25, "inv": 25},
            rng,
        )
        for _ in range(400)
    }
    assert kinds == {"sv_del", "sv_ins", "sv_dup", "sv_inv"}


# --- max_variant_length ---------------------------------------------------


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("kind", "expected"),
    [
        ("snp", 1),
        ("mnp", variant_types.SIZE_RANGES["mnp"][1]),
        ("ins_small", variant_types.SIZE_RANGES["indel"][1] + 1),
        ("del_small", variant_types.SIZE_RANGES["indel"][1] + 1),
        ("sv_del", variant_types.SIZE_RANGES["sv"][1] + 1),
        ("sv_ins", 1),
        ("sv_dup", variant_types.SIZE_RANGES["sv"][1] + 1),
        ("sv_inv", variant_types.SIZE_RANGES["sv"][1] + 1),
    ],
)
def test_max_variant_length(kind, expected):
    assert variant_types.max_variant_length(kind) == expected


# --- variant generators ---------------------------------------------------


@pytest.fixture(autouse=True)
def _seed_fastrand():
    fastrand.pcg32_seed(42)


@pytest.mark.generate_vcf
def test_gen_snp_shape():
    ref, alt, info = variant_types.gen_snp(100, reference_data=None)
    assert len(ref) == 1
    assert len(alt) == 1
    assert ref != alt
    assert ref in variant_types.ALLELES and alt in variant_types.ALLELES
    assert info == {}


@pytest.mark.generate_vcf
def test_gen_mnp_shape():
    ref, alt, info = variant_types.gen_mnp(100, reference_data=None)
    low, high = variant_types.SIZE_RANGES["mnp"]
    assert low <= len(ref) <= high
    assert len(ref) == len(alt)
    for r, a in zip(ref, alt):
        assert r != a
    assert info == {}


@pytest.mark.generate_vcf
def test_gen_ins_small_shape():
    ref, alt, info = variant_types.gen_ins_small(100, reference_data=None)
    assert len(ref) == 1
    low, high = variant_types.SIZE_RANGES["indel"]
    assert 1 + low <= len(alt) <= 1 + high
    assert alt[0] == ref
    assert info == {}


@pytest.mark.generate_vcf
def test_gen_del_small_shape():
    ref, alt, info = variant_types.gen_del_small(100, reference_data=None)
    low, high = variant_types.SIZE_RANGES["indel"]
    assert 1 + low <= len(ref) <= 1 + high
    assert len(alt) == 1
    assert ref[0] == alt
    assert info == {}


@pytest.mark.generate_vcf
@pytest.mark.parametrize(
    ("gen", "svtype", "expect_end"),
    [
        (variant_types.gen_sv_del, "DEL", True),
        (variant_types.gen_sv_ins, "INS", False),
        (variant_types.gen_sv_dup, "DUP", True),
        (variant_types.gen_sv_inv, "INV", True),
    ],
)
def test_gen_sv_shape(gen, svtype, expect_end):
    position = 1000
    ref, alt, info = gen(position, reference_data=None)
    assert len(ref) == 1 and ref in variant_types.ALLELES
    assert alt == f"<{svtype}>"
    assert info["SVTYPE"] == svtype
    assert "SVLEN" in info
    svlen = int(info["SVLEN"])
    low, high = variant_types.SIZE_RANGES["sv"]
    assert low <= abs(svlen) <= high
    if expect_end:
        assert "END" in info
        assert int(info["END"]) == position + abs(svlen)
    else:
        assert "END" not in info
