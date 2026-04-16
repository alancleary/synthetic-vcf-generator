from __future__ import annotations

import random

import fastrand


ALLELES = ("A", "C", "G", "T")

SIZE_RANGES = {
    "mnp": (2, 5),
    "indel": (1, 50),
    "sv": (50, 10_000),
}

VALID_TYPE_KEYS = ("snp", "mnp", "indel", "sv")
VALID_INDEL_KEYS = ("ins", "del")
VALID_SV_KEYS = ("del", "ins", "dup", "inv")

DEFAULT_TYPE_WEIGHTS = {"snp": 80, "mnp": 5, "indel": 10, "sv": 5}
DEFAULT_INDEL_WEIGHTS = {"ins": 50, "del": 50}
DEFAULT_SV_WEIGHTS = {"del": 40, "ins": 20, "dup": 20, "inv": 20}

SV_HEADER_LINES = (
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
)


def parse_weights(
    value: str | None,
    valid_keys: tuple[str, ...],
    default: dict[str, int],
) -> dict[str, int]:
    """
    Parse a CSV of `key=int` pairs into a dict of weights.

    Rules:
      - Each value must be a non-negative integer.
      - Sum must equal exactly 100.
      - Unspecified keys default to 0.
      - Unknown keys raise ValueError.
      - If `value` is None, return a copy of `default`.
    """
    if value is None:
        return dict(default)

    weights = {k: 0 for k in valid_keys}
    seen: set[str] = set()
    for pair in value.split(","):
        pair = pair.strip()
        if not pair:
            raise ValueError(f"Empty weight entry in '{value}'")
        if "=" not in pair:
            raise ValueError(f"Invalid weight syntax: '{pair}' (expected 'key=value')")
        key, _, raw = pair.partition("=")
        key = key.strip()
        raw = raw.strip()
        if key not in valid_keys:
            raise ValueError(
                f"Unknown weight key: '{key}'. Valid keys: {', '.join(valid_keys)}"
            )
        if key in seen:
            raise ValueError(f"Duplicate weight key: '{key}'")
        seen.add(key)
        try:
            parsed = int(raw)
        except ValueError:
            raise ValueError(f"Weight value must be an integer: '{key}={raw}'")
        if parsed < 0:
            raise ValueError(f"Weight value must be non-negative: '{key}={raw}'")
        weights[key] = parsed

    total = sum(weights.values())
    if total != 100:
        raise ValueError(f"Weights must sum to 100, got {total}: '{value}'")
    return weights


def pick_variant_kind(
    type_weights: dict[str, int],
    indel_weights: dict[str, int],
    sv_weights: dict[str, int],
    rng: random.Random,
) -> str:
    """Pick a leaf variant kind based on hierarchical weights."""
    top_keys = list(type_weights.keys())
    top_vals = list(type_weights.values())
    top = rng.choices(top_keys, weights=top_vals, k=1)[0]
    if top == "snp":
        return "snp"
    if top == "mnp":
        return "mnp"
    if top == "indel":
        sub_keys = list(indel_weights.keys())
        sub_vals = list(indel_weights.values())
        sub = rng.choices(sub_keys, weights=sub_vals, k=1)[0]
        return f"{sub}_small"
    sub_keys = list(sv_weights.keys())
    sub_vals = list(sv_weights.values())
    sub = rng.choices(sub_keys, weights=sub_vals, k=1)[0]
    return f"sv_{sub}"


def max_variant_length(kind: str) -> int:
    """
    Conservative upper bound on the span (in bases) the variant can occupy on
    the reference, used to avoid running past the chromosome end.
    """
    if kind == "snp":
        return 1
    if kind == "mnp":
        return SIZE_RANGES["mnp"][1]
    if kind in ("ins_small", "del_small"):
        return SIZE_RANGES["indel"][1] + 1
    if kind == "sv_ins":
        # Pure insertion — only the anchor base is on the reference.
        return 1
    # sv_del / sv_dup / sv_inv each span anchor + SVLEN bases.
    return SIZE_RANGES["sv"][1] + 1


def _random_base() -> str:
    return ALLELES[fastrand.pcg32randint(0, 3)]


def _random_different_base(ref: str) -> str:
    candidates = [a for a in ALLELES if a != ref]
    return candidates[fastrand.pcg32randint(0, len(candidates) - 1)]


def _random_bases(length: int) -> str:
    return "".join(_random_base() for _ in range(length))


def _get_ref_bases(position: int, length: int, reference_data) -> str:
    if reference_data is not None:
        bases = reference_data.get_ref_at_pos(position - 1, length=length).upper()
        return bases
    return _random_bases(length)


def _safe_anchor(ref: str) -> str:
    return ref if ref in ALLELES else _random_base()


def gen_snp(position, reference_data):
    ref = _get_ref_bases(position, 1, reference_data)
    if ref in ALLELES:
        alt = _random_different_base(ref)
    else:
        alt = _random_base()
    return ref, alt, {}


def gen_mnp(position, reference_data):
    length = fastrand.pcg32randint(*SIZE_RANGES["mnp"])
    ref = _get_ref_bases(position, length, reference_data)
    alt_chars = []
    for r in ref:
        if r in ALLELES:
            alt_chars.append(_random_different_base(r))
        else:
            alt_chars.append(_random_base())
    return ref, "".join(alt_chars), {}


def gen_ins_small(position, reference_data):
    length = fastrand.pcg32randint(*SIZE_RANGES["indel"])
    anchor = _safe_anchor(_get_ref_bases(position, 1, reference_data))
    alt = anchor + _random_bases(length)
    return anchor, alt, {}


def gen_del_small(position, reference_data):
    length = fastrand.pcg32randint(*SIZE_RANGES["indel"])
    ref = _get_ref_bases(position, length + 1, reference_data)
    anchor = _safe_anchor(ref[0])
    ref = anchor + ref[1:]
    return ref, anchor, {}


def _gen_sv_common(position, reference_data):
    length = fastrand.pcg32randint(*SIZE_RANGES["sv"])
    anchor = _safe_anchor(_get_ref_bases(position, 1, reference_data))
    return length, anchor


def gen_sv_del(position, reference_data):
    length, anchor = _gen_sv_common(position, reference_data)
    return (
        anchor,
        "<DEL>",
        {"SVTYPE": "DEL", "END": str(position + length), "SVLEN": f"-{length}"},
    )


def gen_sv_ins(position, reference_data):
    length, anchor = _gen_sv_common(position, reference_data)
    return anchor, "<INS>", {"SVTYPE": "INS", "SVLEN": str(length)}


def gen_sv_dup(position, reference_data):
    length, anchor = _gen_sv_common(position, reference_data)
    return (
        anchor,
        "<DUP>",
        {"SVTYPE": "DUP", "END": str(position + length), "SVLEN": str(length)},
    )


def gen_sv_inv(position, reference_data):
    length, anchor = _gen_sv_common(position, reference_data)
    return (
        anchor,
        "<INV>",
        {"SVTYPE": "INV", "END": str(position + length), "SVLEN": str(length)},
    )


VARIANT_GENERATORS = {
    "snp": gen_snp,
    "mnp": gen_mnp,
    "ins_small": gen_ins_small,
    "del_small": gen_del_small,
    "sv_del": gen_sv_del,
    "sv_ins": gen_sv_ins,
    "sv_dup": gen_sv_dup,
    "sv_inv": gen_sv_inv,
}
