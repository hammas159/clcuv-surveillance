"""Tests for the banded centre-star aligner."""

import random

import pytest

from clcuv.align import (
    GAP,
    AlignmentImpossible,
    Scoring,
    align,
    align_pair,
    alignment_report,
    check_comparable,
    choose_centre,
)

BASES = "ACGT"


def _random_sequence(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(n))


def _mutate(sequence: str, *, sites: int, seed: int) -> str:
    rng = random.Random(seed)
    chars = list(sequence)
    for position in rng.sample(range(len(chars)), sites):
        chars[position] = rng.choice([b for b in BASES if b != chars[position]])
    return "".join(chars)


# --- pairwise -------------------------------------------------------------


def test_identical_sequences_align_without_gaps():
    sequence = _random_sequence(300, seed=1)
    top, bottom = align_pair(sequence, sequence)
    assert top == bottom == sequence
    assert GAP not in top


def test_substitutions_do_not_introduce_gaps():
    a = _random_sequence(400, seed=2)
    b = _mutate(a, sites=12, seed=3)
    top, bottom = align_pair(a, b)
    assert len(top) == len(bottom) == 400
    assert GAP not in top and GAP not in bottom
    assert sum(1 for x, y in zip(top, bottom, strict=True) if x != y) == 12


def test_a_deletion_is_placed_as_a_gap():
    a = _random_sequence(300, seed=4)
    b = a[:150] + a[156:]  # six bases removed from the middle
    top, bottom = align_pair(a, b)
    assert len(top) == len(bottom)
    assert bottom.count(GAP) == 6
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_an_insertion_gaps_the_other_row():
    a = _random_sequence(300, seed=5)
    b = a[:100] + "ACGTACGT" + a[100:]
    top, bottom = align_pair(a, b)
    assert top.count(GAP) == 8
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_original_sequences_are_always_recoverable():
    """The invariant that makes an alignment an alignment: no base is lost."""
    a = _random_sequence(250, seed=6)
    b = _mutate(a[:120] + a[125:], sites=8, seed=7)
    top, bottom = align_pair(a, b)
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_ambiguity_codes_are_neither_rewarded_nor_punished():
    scoring = Scoring()
    assert scoring.pair("A", "A") == scoring.match
    assert scoring.pair("A", "C") == scoring.mismatch
    assert scoring.pair("N", "A") == 0
    assert scoring.pair("A", "N") == 0


def test_empty_input_degrades_rather_than_crashing():
    top, bottom = align_pair("", "ACGT")
    assert top == GAP * 4
    assert bottom == "ACGT"


# --- the guard ------------------------------------------------------------


def test_rotated_circular_genomes_are_refused():
    """The failure this guard exists for: two rotations are biologically identical
    and share no aligned column, and nothing else in the pipeline would notice."""
    sequence = _random_sequence(600, seed=8)
    rotated = sequence[300:] + sequence[:300]
    with pytest.raises(AlignmentImpossible):
        check_comparable([sequence, sequence, rotated, rotated, rotated])


def test_co_oriented_genomes_pass_the_guard():
    a = _random_sequence(600, seed=9)
    check_comparable([a, _mutate(a, sites=30, seed=10), _mutate(a, sites=45, seed=11)])


def test_wildly_different_lengths_are_refused():
    a = _random_sequence(600, seed=12)
    with pytest.raises(AlignmentImpossible):
        check_comparable([a, a, a[:400]])


def test_a_single_sequence_is_not_an_alignment():
    with pytest.raises(AlignmentImpossible):
        check_comparable(["ACGT"])


def test_one_odd_sequence_among_many_does_not_trip_the_guard():
    """80% agreement is the bar: one bad submission should not block the analysis."""
    a = _random_sequence(600, seed=13)
    others = [_mutate(a, sites=20, seed=s) for s in range(14, 24)]
    rotated = a[300:] + a[:300]
    check_comparable([a, *others, rotated])


# --- multiple alignment ---------------------------------------------------


def test_centre_is_the_sequence_nearest_the_median_length():
    assert choose_centre(["A" * 100, "A" * 300, "A" * 305]) == 1


def test_multiple_alignment_returns_equal_lengths():
    a = _random_sequence(400, seed=30)
    sequences = [a, _mutate(a, sites=10, seed=31), a[:200] + a[205:], a + "ACGTA"]
    aligned = align(sequences)
    assert len({len(row) for row in aligned}) == 1


def test_multiple_alignment_preserves_every_base():
    a = _random_sequence(400, seed=32)
    sequences = [a, _mutate(a, sites=10, seed=33), a[:200] + a[205:], a + "ACGTA"]
    for original, row in zip(sequences, align(sequences), strict=True):
        assert row.replace(GAP, "") == original


def test_insertions_in_different_sequences_get_separate_columns():
    """The merge step's actual job: two sequences with insertions at the same place
    must not be forced to share one column, or a deletion appears from nowhere."""
    a = _random_sequence(300, seed=34)
    sequences = [a, a[:150] + "TTTT" + a[150:], a[:150] + "GGGGGG" + a[150:]]
    aligned = align(sequences)
    assert all(row.replace(GAP, "") == s for row, s in zip(aligned, sequences, strict=True))
    assert len(aligned[0]) >= 306


def test_identical_sequences_align_to_themselves():
    a = _random_sequence(200, seed=35)
    aligned = align([a, a, a])
    assert aligned == [a, a, a]


def test_report_counts_invariant_columns():
    a = _random_sequence(200, seed=36)
    report = alignment_report(align([a, a, a]))
    assert report["sequences"] == 3
    assert report["width"] == 200
    assert report["invariant_columns"] == 200
    assert report["invariant_fraction"] == 1.0
    assert report["gap_fraction"] == 0.0


def test_report_notices_variation():
    a = _random_sequence(300, seed=37)
    report = alignment_report(align([a, _mutate(a, sites=30, seed=38)]))
    assert report["invariant_columns"] == 270
    assert report["invariant_fraction"] < 1.0


def test_report_on_nothing():
    assert alignment_report([]) == {"sequences": 0}
