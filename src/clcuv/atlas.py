"""Mutation atlas: what is changing, where, and whether it is spreading.

A mutation atlas is not a list of mutations. Every genome differs from the reference at
dozens of positions and almost none of them matter. The question surveillance exists to
answer is narrower:

    **which variant is rising in frequency, and how fast?**

A mutation at 40% frequency that has sat at 40% for three seasons is background. A
mutation that went 2% → 8% → 25% in three seasons is a lineage winning, and it is worth
knowing about while it is still at 25%.

So the atlas is organised around **trajectories**, not snapshots, and the alerting is on
the derivative rather than the level.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Isolate:
    """One sequenced sample, with the metadata that makes it surveillance."""

    name: str
    sequence: str
    # Collection period - a season, month or year. Ordered lexically, so use a sortable
    # form such as "2026-Q1".
    period: str = ""
    location: str = ""
    host: str = ""


@dataclass
class Variant:
    position: int          # 0-based alignment column
    reference: str
    alternate: str
    counts_by_period: dict[str, int] = field(default_factory=dict)
    totals_by_period: dict[str, int] = field(default_factory=dict)
    locations: Counter = field(default_factory=Counter)

    @property
    def label(self) -> str:
        return f"{self.reference}{self.position + 1}{self.alternate}"

    def frequency(self, period: str) -> float:
        total = self.totals_by_period.get(period, 0)
        return round(self.counts_by_period.get(period, 0) / total, 6) if total else 0.0

    def trajectory(self) -> list[tuple[str, float, int]]:
        """(period, frequency, sample size), in period order."""
        return [
            (p, self.frequency(p), self.totals_by_period.get(p, 0))
            for p in sorted(self.totals_by_period)
        ]

    @property
    def overall_frequency(self) -> float:
        total = sum(self.totals_by_period.values())
        return round(sum(self.counts_by_period.values()) / total, 6) if total else 0.0


def build_atlas(
    isolates: Sequence[Isolate], reference: str, *, min_frequency: float = 0.01
) -> list[Variant]:
    """Every variant against the reference, with its trajectory.

    `min_frequency` drops singletons. In a set of a few hundred genomes, a variant seen
    once is more likely a sequencing error than a lineage, and an atlas full of
    sequencing errors is one nobody reads.
    """
    if not isolates:
        return []

    length = len(reference)
    for isolate in isolates:
        if len(isolate.sequence) != length:
            raise ValueError(
                f"{isolate.name}: length {len(isolate.sequence)} does not match the "
                f"reference ({length}) - sequences must be aligned"
            )

    variants: dict[tuple[int, str], Variant] = {}
    period_totals_by_position: dict[int, Counter] = defaultdict(Counter)

    # Denominators are counted per position, not per isolate: a genome with an N at a
    # position contributes no information there and must not inflate the denominator.
    for isolate in isolates:
        for i, base in enumerate(isolate.sequence.upper()):
            if base in "ACGT":
                period_totals_by_position[i][isolate.period] += 1

    for isolate in isolates:
        for i, base in enumerate(isolate.sequence.upper()):
            ref = reference[i].upper()
            if base not in "ACGT" or ref not in "ACGT" or base == ref:
                continue
            key = (i, base)
            variant = variants.setdefault(key, Variant(position=i, reference=ref, alternate=base))
            variant.counts_by_period[isolate.period] = (
                variant.counts_by_period.get(isolate.period, 0) + 1
            )
            variant.locations[isolate.location] += 1

    for (position, _), variant in variants.items():
        variant.totals_by_period = dict(period_totals_by_position[position])

    return sorted(
        (v for v in variants.values() if v.overall_frequency >= min_frequency),
        key=lambda v: (-v.overall_frequency, v.position),
    )


def two_proportion_z(
    successes_a: int, total_a: int, successes_b: int, total_b: int
) -> float:
    """Z-statistic for a difference between two observed proportions.

    Necessary because frequency estimates are noisy and the noise scales with sample
    size. Two seasons of 100 genomes each can easily differ by ten points with nothing
    happening at all, and a detector that fires on that is one nobody trusts by the
    third false alarm.
    """
    if total_a <= 0 or total_b <= 0:
        return 0.0
    pooled = (successes_a + successes_b) / (total_a + total_b)
    if pooled in (0.0, 1.0):
        return 0.0
    standard_error = math.sqrt(pooled * (1 - pooled) * (1 / total_a + 1 / total_b))
    if standard_error == 0:
        return 0.0
    return (successes_b / total_b - successes_a / total_a) / standard_error


@dataclass
class Emerging:
    variant: Variant
    first_period: str
    last_period: str
    first_frequency: float
    last_frequency: float
    change: float
    fold: float | None
    z: float = 0.0

    def summary(self) -> dict:
        return {
            "variant": self.variant.label,
            "from": f"{self.first_frequency:.1%} ({self.first_period})",
            "to": f"{self.last_frequency:.1%} ({self.last_period})",
            "change": round(self.change, 4),
            "fold": self.fold,
            "z": round(self.z, 3),
            "locations": dict(self.variant.locations.most_common(5)),
        }


def emerging_variants(
    variants: Sequence[Variant], *, min_change: float = 0.10, min_samples: int = 10,
    min_z: float = 1.96,
) -> list[Emerging]:
    """Variants whose frequency is rising, and rising by more than chance.

    Two filters, both necessary and neither sufficient.

    `min_samples` excludes periods with too few sequences. A variant at "100%" in a
    period where three genomes were sequenced is not at 100% — it is unmeasured.

    `min_z` requires the rise to be statistically distinguishable from sampling noise.
    Without it the detector fires on the ordinary wobble of two finite samples: this
    was a real false positive here, where a variant sitting at a constant 40% was
    reported as rising because two seasons of 100 genomes happened to land at 33% and
    45%. Effect size alone is not evidence.
    """
    out: list[Emerging] = []

    for variant in variants:
        usable = [
            (period, freq, n) for period, freq, n in variant.trajectory()
            if n >= min_samples
        ]
        if len(usable) < 2:
            continue

        (first_period, first_freq, first_n) = usable[0]
        (last_period, last_freq, last_n) = usable[-1]

        change = last_freq - first_freq
        if change < min_change:
            continue

        z = two_proportion_z(
            variant.counts_by_period.get(first_period, 0), first_n,
            variant.counts_by_period.get(last_period, 0), last_n,
        )
        if z < min_z:
            continue

        out.append(Emerging(
            variant=variant, first_period=first_period, last_period=last_period,
            first_frequency=first_freq, last_frequency=last_freq, change=change,
            # Fold change is undefined from zero. Reporting it as infinite, or as a
            # large number, turns "newly detected" into "exploding" — a different claim.
            fold=round(last_freq / first_freq, 3) if first_freq > 0 else None,
            z=round(z, 4),
        ))

    return sorted(out, key=lambda e: -e.change)


def geographic_spread(variant: Variant) -> dict:
    """Where a variant has been seen, and how concentrated it is.

    A variant confined to one district is a local lineage; the same variant across five
    districts is spreading, and the difference changes what anyone should do about it.
    """
    total = sum(variant.locations.values())
    if not total:
        return {"locations": 0, "concentration": 0.0}

    shares = [n / total for n in variant.locations.values()]
    return {
        "variant": variant.label,
        "locations": len(variant.locations),
        "top": dict(variant.locations.most_common(5)),
        # Herfindahl index: 1.0 means a single location, low means widely spread.
        "concentration": round(sum(s * s for s in shares), 6),
        "spreading": len(variant.locations) >= 3 and sum(s * s for s in shares) < 0.5,
    }
