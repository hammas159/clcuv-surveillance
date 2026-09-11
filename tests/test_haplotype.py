"""Tests for clonal collapse and the stratified emergence test.

Both exist because of false positives found on real GenBank data, so each is tested by
first reproducing the false positive and then showing the control removes it. A fix
demonstrated only on data where nothing was wrong has not been demonstrated.
"""

from __future__ import annotations

from clcuv.atlas import (
    Isolate,
    build_atlas,
    emerging_variants,
    rises_within_locations,
)
from clcuv.haplotype import collapse_clonal, effective_sample_sizes, identity

REFERENCE = "ACGT" * 25  # 100 bp


def _isolate(name, sequence, period="2020", location="Punjab"):
    return Isolate(name=name, sequence=sequence, period=period, location=location)


def _with_variant(position: int, base: str) -> str:
    chars = list(REFERENCE)
    chars[position] = base
    return "".join(chars)


def _tagged(index: int, *, variant: bool) -> str:
    """A genome unique to one isolate, optionally carrying the variant at column 61."""
    chars = list(REFERENCE)
    position = 80 + index
    chars[position] = "C" if chars[position] != "C" else "A"
    if variant:
        chars[60] = "T"
    return "".join(chars)


# --- identity -------------------------------------------------------------


def test_identical_sequences_score_one():
    assert identity(REFERENCE, REFERENCE) == 1.0


def test_one_difference_in_a_hundred():
    assert identity(REFERENCE, _with_variant(10, "T")) == 0.99


def test_gaps_and_ambiguity_are_skipped_not_counted_as_differences():
    """A partially sequenced genome must not look divergent because of what is
    missing, or every short submission becomes a new haplotype."""
    assert identity(REFERENCE, "N" * 50 + REFERENCE[50:]) == 1.0


def test_nothing_comparable_is_zero_rather_than_an_error():
    assert identity(REFERENCE, "-" * 100) == 0.0
    assert identity("", "") == 0.0


# --- collapse -------------------------------------------------------------


def test_a_clonal_batch_collapses_to_one_observation():
    """The real case: eight genomes, one submission, one haplotype."""
    isolates = [_isolate(f"ON31278{i}", REFERENCE, period="2021") for i in range(8)]
    kept, report = collapse_clonal(isolates)
    assert len(kept) == 1
    assert (report.isolates_in, report.isolates_out) == (8, 1)
    assert report.largest_clone == 8
    assert report.inflation == 8.0


def test_genuinely_distinct_sequences_survive():
    isolates = [
        _isolate("a", REFERENCE),
        _isolate("b", _with_variant(10, "T")),
        _isolate("c", _with_variant(20, "T")),
    ]
    kept, report = collapse_clonal(isolates, threshold=0.995)
    assert len(kept) == 3
    assert report.inflation == 1.0


def test_identical_sequences_in_different_places_are_kept():
    """The line between duplication and signal: one haplotype in Punjab and in Sindh
    is two observations of spread, not one sample counted twice."""
    isolates = [
        _isolate("a", REFERENCE, location="Punjab"),
        _isolate("b", REFERENCE, location="Sindh"),
        _isolate("c", REFERENCE, location="KPK"),
    ]
    assert len(collapse_clonal(isolates)[0]) == 3


def test_identical_sequences_in_different_years_are_kept():
    isolates = [
        _isolate("a", REFERENCE, period="2019"),
        _isolate("b", REFERENCE, period="2021"),
    ]
    assert len(collapse_clonal(isolates)[0]) == 2


def test_the_threshold_decides_what_counts_as_the_same_infection():
    near = _with_variant(10, "T")  # 99% identical
    isolates = [_isolate("a", REFERENCE), _isolate("b", near)]
    assert len(collapse_clonal(isolates, threshold=0.98)[0]) == 1
    assert len(collapse_clonal(isolates, threshold=0.999)[0]) == 2


def test_the_report_names_the_worst_strata():
    isolates = [_isolate(f"x{i}", REFERENCE, period="2021") for i in range(11)]
    isolates.append(_isolate("y", _with_variant(4, "T"), period="2019"))
    _, report = collapse_clonal(isolates)
    summary = report.summary()
    assert summary["isolates"] == 12
    assert summary["haplotypes"] == 2
    assert summary["largest_clonal_group"] == 11
    assert summary["worst_strata"] == {"2021/Punjab": "11 -> 1"}


def test_effective_sample_size_flags_strata_too_clonal_to_test():
    isolates = [_isolate(f"x{i}", REFERENCE, period="2021") for i in range(8)]
    isolates += [_isolate(f"y{i}", _with_variant(i, "T"), period="2019") for i in range(6)]
    sizes = effective_sample_sizes(isolates)
    assert sizes[("2021", "Punjab")] == {
        "sequences": 8,
        "haplotypes": 1,
        "inflation": 8.0,
        "usable_for_statistics": False,
    }
    assert sizes[("2019", "Punjab")]["usable_for_statistics"] is True


def test_collapse_of_nothing():
    kept, report = collapse_clonal([])
    assert kept == []
    assert report.inflation == 0.0


# --- the false positive, end to end ---------------------------------------


def _clonal_survey():
    """A reconstruction of the real dataset's shape.

    2019: eight genomes, one submission, all identical, none carrying the variant.
    2021: eight genomes, one submission, all identical, all carrying it.

    Pooled that reads 0/8 -> 8/8, a hundred-point rise at z = 4. It is also two
    isolates, and two isolates cannot support a claim about a province.
    """
    return [_isolate(f"MW1834{i:02d}", REFERENCE, period="2019") for i in range(8)] + [
        _isolate(f"ON3127{i:02d}", _with_variant(40, "T"), period="2021") for i in range(8)
    ]


def test_raw_counts_produce_the_false_positive():
    """Shows the bug is real before testing the fix. A control demonstrated only to be
    harmless has not been demonstrated to be necessary."""
    found = emerging_variants(build_atlas(_clonal_survey(), REFERENCE), min_samples=8)
    assert [e.variant.label for e in found] == ["A41T"]
    assert found[0].change == 1.0
    assert found[0].z > 3


def test_collapsing_first_removes_it():
    kept, report = collapse_clonal(_clonal_survey())
    assert report.inflation == 8.0
    assert emerging_variants(build_atlas(kept, REFERENCE), min_samples=8) == []


def test_collapsing_keeps_a_rise_that_has_independent_support():
    """The control must not simply suppress everything. Twelve distinct genomes per
    year are twelve observations, and the signal in them survives."""
    isolates = [_isolate(f"a{i}", _tagged(i, variant=False), period="2019") for i in range(12)]
    isolates += [_isolate(f"b{i}", _tagged(i, variant=True), period="2021") for i in range(12)]

    kept, report = collapse_clonal(isolates)
    assert report.inflation == 1.0
    found = emerging_variants(build_atlas(kept, REFERENCE), min_samples=8)
    assert "A61T" in [e.variant.label for e in found]


# --- stratification -------------------------------------------------------


def _geographic_confounder():
    """The other confounder the real data contained: the survey moved province.

    A variant fixed in Sindh and absent from Punjab appears to emerge when 2021 adds
    Sindh to a survey that was Punjab-only in 2019. Nothing changed but the map.
    """
    isolates = [
        _isolate(f"p19_{i}", _with_variant(i % 8, "A"), period="2019", location="Punjab")
        for i in range(10)
    ]
    isolates += [
        _isolate(f"p21_{i}", _with_variant(i % 8, "A"), period="2021", location="Punjab")
        for i in range(10)
    ]
    isolates += [
        _isolate(f"s21_{i}", _with_variant(50, "A"), period="2021", location="Sindh")
        for i in range(10)
    ]
    return isolates


def test_the_pooled_test_is_fooled_by_geography():
    found = emerging_variants(build_atlas(_geographic_confounder(), REFERENCE))
    assert any(e.variant.label == "G51A" for e in found)


def test_stratifying_discards_it():
    found = emerging_variants(build_atlas(_geographic_confounder(), REFERENCE), stratify=True)
    assert not any(e.variant.label == "G51A" for e in found)


def test_a_rise_inside_one_location_is_confirmed():
    isolates = [_isolate(f"a{i}", REFERENCE, period="2019") for i in range(10)]
    isolates += [_isolate(f"b{i}", _with_variant(50, "A"), period="2021") for i in range(10)]
    found = emerging_variants(build_atlas(isolates, REFERENCE), stratify=True)
    assert [e.variant.label for e in found] == ["G51A"]
    assert found[0].confirmed_in == ["Punjab"]
    assert found[0].stratified is True


def test_confirmation_requires_both_periods_present_in_the_location():
    """A location sampled in one year alone can say nothing about a trend and must not
    be reported as confirming one."""
    atlas = build_atlas(_geographic_confounder(), REFERENCE)
    sindh_only = next(v for v in atlas if v.label == "G51A")
    assert rises_within_locations(sindh_only) == []


def test_summary_exposes_whether_it_was_confirmed():
    isolates = [_isolate(f"a{i}", REFERENCE, period="2019") for i in range(10)]
    isolates += [_isolate(f"b{i}", _with_variant(50, "A"), period="2021") for i in range(10)]
    summary = emerging_variants(build_atlas(isolates, REFERENCE), stratify=True)[0].summary()
    assert summary["variant"] == "G51A"
    assert summary["confirmed_in"] == ["Punjab"]
