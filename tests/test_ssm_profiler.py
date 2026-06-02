import math

import pytest

from uht_tooling.workflows.ssm_profiler import (
    parse_scheme_map,
    parse_target_sites,
    theoretical_aa_distribution,
)


def test_parse_target_sites_sorts_and_deduplicates():
    assert parse_target_sites([7, 3, 7, 1]) == [1, 3, 7]


def test_parse_scheme_map_parses_multiple_entries():
    assert parse_scheme_map(["45:NNK", "46:NNW"]) == {45: "NNK", 46: "NNW"}


def test_parse_scheme_map_rejects_invalid_codes():
    with pytest.raises(ValueError):
        parse_scheme_map(["45:NNZ"])


def test_theoretical_aa_distribution_nnk_sums_to_one():
    dist = theoretical_aa_distribution("NNK")
    assert math.isclose(sum(dist.values()), 1.0)
    assert "*" in dist


def test_theoretical_aa_distribution_atg_is_single_methionine():
    assert theoretical_aa_distribution("ATG") == {"M": 1.0}
