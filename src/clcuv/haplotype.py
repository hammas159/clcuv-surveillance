"""Collapsing clonal sequences, because eight identical genomes are one observation.

The failure this exists for was found on real GenBank data, not imagined.

A survey of Cotton leaf curl Multan virus in Punjab looked like this:

    2019   8 genomes   all from submission MW183409-MW183416
    2020  11 genomes   all from submission MW654014-MW654024
    2021   8 genomes   all from submission ON312781-ON312788

and **every one of the 28 pairs within the 2021 set was 100% identical**. Eight
sequences, one haplotype, one field, one submission. The effective sample size is one.

Fed to a two-proportion z-test as n=8 against n=8, that produced nine variants
"significantly emerging" between 2019 and 2021 at z > 3. Every number was computed
correctly. The conclusion was worthless, because the test was told there were sixteen
independent observations when there were two.

This is exactly the failure a research agent has when it counts a wire story
republished by twelve outlets as twelve corroborating sources. Same bug, different
field: **the unit of replication is not the row**.

Collapsing happens within a (period, location) stratum, which is where clonal expansion
and repeated sampling of one field actually inflate the count. Identical haplotypes in
*different* places or times are left alone — that is signal, not duplication.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from .atlas import Isolate


def identity(a: str, b: str) -> float:
    """Fraction of comparable sites that agree.

    Gaps and ambiguity codes are skipped rather than counted as differences — a
    partially sequenced genome should not look divergent because of what is missing.
    """
    compared = same = 0
    # strict=False: a partial genome is still comparable over the part it covers.
    for x, y in zip(a.upper(), b.upper(), strict=False):
        if x not in "ACGT" or y not in "ACGT":
            continue
        compared += 1
        if x == y:
            same += 1
    return same / compared if compared else 0.0


@dataclass
class Haplotype:
    """One distinct sequence within a stratum, and every isolate carrying it."""

    representative: Isolate
    members: list[Isolate] = field(default_factory=list)

    @property
    def count(self) -> int:
        return 1 + len(self.members)

    @property
    def accessions(self) -> list[str]:
        return [self.representative.name] + [m.name for m in self.members]


@dataclass
class CollapseReport:
    isolates_in: int
    isolates_out: int
    strata: int
    largest_clone: int
    by_stratum: dict[tuple[str, str], tuple[int, int]] = field(default_factory=dict)

    @property
    def inflation(self) -> float:
        """How much the raw count overstated the independent evidence."""
        return round(self.isolates_in / self.isolates_out, 3) if self.isolates_out else 0.0

    def summary(self) -> dict:
        return {
            "isolates": self.isolates_in,
            "haplotypes": self.isolates_out,
            "inflation": self.inflation,
            "largest_clonal_group": self.largest_clone,
            "strata": self.strata,
            "worst_strata": {
                f"{period}/{location}": f"{raw} -> {kept}"
                for (period, location), (raw, kept) in sorted(
                    self.by_stratum.items(), key=lambda kv: kv[1][1] - kv[1][0]
                )[:5]
                if raw != kept
            },
        }


def collapse_clonal(
    isolates: Sequence[Isolate], *, threshold: float = 0.999
) -> tuple[list[Isolate], CollapseReport]:
    """Reduce near-identical sequences within each (period, location) to one.

    `threshold` defaults to 0.999 — on a 2,800 bp genome that is roughly "differing at
    at most two sites". Sequences that close within one field and one season are the
    same infection sampled repeatedly, not independent evidence about the population.

    Returns the collapsed isolates and a report, because silently discarding five sixths
    of a dataset without saying so would be its own kind of dishonesty.
    """
    by_stratum: dict[tuple[str, str], list[Isolate]] = {}
    for isolate in isolates:
        by_stratum.setdefault((isolate.period, isolate.location), []).append(isolate)

    kept: list[Isolate] = []
    report = CollapseReport(
        isolates_in=len(isolates), isolates_out=0, strata=len(by_stratum), largest_clone=0
    )

    for stratum, members in by_stratum.items():
        haplotypes: list[Haplotype] = []
        for isolate in members:
            for haplotype in haplotypes:
                if identity(isolate.sequence, haplotype.representative.sequence) >= threshold:
                    haplotype.members.append(isolate)
                    break
            else:
                haplotypes.append(Haplotype(representative=isolate))

        report.by_stratum[stratum] = (len(members), len(haplotypes))
        report.largest_clone = max(
            report.largest_clone, max((h.count for h in haplotypes), default=0)
        )
        kept.extend(h.representative for h in haplotypes)

    report.isolates_out = len(kept)
    return kept, report


def effective_sample_sizes(
    isolates: Sequence[Isolate], *, threshold: float = 0.999
) -> dict[tuple[str, str], dict]:
    """Raw against effective n for every stratum.

    Worth printing before any statistical test. A stratum where eleven sequences reduce
    to one cannot support a claim about a population, and the number that makes that
    obvious is the ratio, not either count alone.
    """
    collapsed, report = collapse_clonal(isolates, threshold=threshold)
    return {
        stratum: {
            "sequences": raw,
            "haplotypes": kept,
            "inflation": round(raw / kept, 2) if kept else 0.0,
            "usable_for_statistics": kept >= 5,
        }
        for stratum, (raw, kept) in sorted(report.by_stratum.items())
    }
