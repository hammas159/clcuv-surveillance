"""Strain classification, with an explicit "I do not recognise this" answer.

A classifier that must choose one of its known labels will assign a novel recombinant
to whichever strain it resembles least-badly, and report it with the same confidence as
a genuine match. For surveillance that is the worst possible failure: **the novel
strain is the one you built the system to find**, and forced-choice classification is
precisely the mechanism that hides it.

So classification returns `None` above a distance threshold. An unassigned genome is a
finding, not a gap.

k-mer profiles rather than alignment, because they are alignment-free: a recombinant
whose segments come from two parents has no single good alignment, and that is the
exact case where an alignment-based classifier degrades silently.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field


def kmer_profile(sequence: str, k: int = 6) -> Counter:
    """Counts of each k-mer. Ambiguous windows are skipped, not imputed."""
    seq = sequence.upper()
    counts: Counter = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if set(kmer) <= set("ACGT"):
            counts[kmer] += 1
    return counts


def cosine_distance(a: Counter, b: Counter) -> float:
    """1 − cosine similarity between two k-mer profiles.

    Cosine rather than Euclidean because it is length-invariant: a partial genome
    should be classified by its composition, not pushed away from every centroid
    because it is short.
    """
    if not a or not b:
        return 1.0
    shared = set(a) & set(b)
    dot = sum(a[k] * b[k] for k in shared)
    norm_a = math.sqrt(sum(v * v for v in a.values()))
    norm_b = math.sqrt(sum(v * v for v in b.values()))
    if not norm_a or not norm_b:
        return 1.0
    return round(1 - dot / (norm_a * norm_b), 6)


@dataclass
class Strain:
    name: str
    centroid: Counter = field(default_factory=Counter)
    members: int = 0
    # Distance spread within the strain, used to calibrate the novelty threshold.
    mean_internal_distance: float = 0.0
    max_internal_distance: float = 0.0


@dataclass
class Classification:
    strain: str | None
    distance: float
    runner_up: str | None = None
    runner_up_distance: float = 1.0
    novel: bool = False
    reason: str = ""

    @property
    def margin(self) -> float:
        """How much better the best match is than the second best.

        A small margin means the genome sits between two strains, which for a
        recombinant is not a tie to be broken but the answer itself.
        """
        return round(self.runner_up_distance - self.distance, 6)


@dataclass
class StrainClassifier:
    k: int = 6
    # Multiple of a strain's own internal spread beyond which a genome is novel.
    novelty_factor: float = 2.0
    # Floor, for strains whose members are nearly identical and whose spread is ~0.
    min_novelty_threshold: float = 0.05
    # Below this margin the genome sits between strains.
    min_margin: float = 0.01
    strains: dict[str, Strain] = field(default_factory=dict)

    def fit(self, labelled: Sequence[tuple[str, str]]) -> StrainClassifier:
        """`labelled` is (strain name, sequence)."""
        by_strain: dict[str, list[Counter]] = {}
        for name, sequence in labelled:
            by_strain.setdefault(name, []).append(kmer_profile(sequence, self.k))

        self.strains = {}
        for name, profiles in by_strain.items():
            centroid: Counter = Counter()
            for profile in profiles:
                centroid.update(profile)

            distances = [cosine_distance(p, centroid) for p in profiles]
            self.strains[name] = Strain(
                name=name, centroid=centroid, members=len(profiles),
                mean_internal_distance=round(sum(distances) / len(distances), 6),
                max_internal_distance=round(max(distances), 6),
            )
        return self

    def threshold_for(self, strain: Strain) -> float:
        return max(
            self.min_novelty_threshold,
            strain.max_internal_distance * self.novelty_factor,
        )

    def classify(self, sequence: str) -> Classification:
        if not self.strains:
            raise ValueError("classifier has not been fitted")

        profile = kmer_profile(sequence, self.k)
        ranked = sorted(
            ((name, cosine_distance(profile, s.centroid)) for name, s in self.strains.items()),
            key=lambda pair: pair[1],
        )

        best_name, best_distance = ranked[0]
        runner_up, runner_up_distance = ranked[1] if len(ranked) > 1 else (None, 1.0)

        result = Classification(
            strain=best_name, distance=best_distance,
            runner_up=runner_up, runner_up_distance=runner_up_distance,
        )

        if best_distance > self.threshold_for(self.strains[best_name]):
            # Unlike anything known. This is the finding, not a failure to classify.
            return Classification(
                strain=None, distance=best_distance, runner_up=best_name,
                runner_up_distance=runner_up_distance, novel=True,
                reason=(
                    f"distance {best_distance:.4f} exceeds the threshold "
                    f"{self.threshold_for(self.strains[best_name]):.4f} for its nearest "
                    f"strain {best_name}"
                ),
            )

        if runner_up is not None and result.margin < self.min_margin:
            return Classification(
                strain=None, distance=best_distance, runner_up=runner_up,
                runner_up_distance=runner_up_distance, novel=True,
                reason=(
                    f"sits between {best_name} and {runner_up} "
                    f"(margin {result.margin:.4f}) - possible recombinant"
                ),
            )

        return result


def recombination_signal(
    query: str, parent_a: str, parent_b: str, *, window: int = 200, step: int = 50
) -> dict:
    """Does similarity switch parent along the genome?

    A point mutation shifts similarity to both parents a little. **Recombination moves
    a block**, so similarity to one parent rises while similarity to the other falls,
    at a specific coordinate. Detecting the switch — rather than the overall distance —
    is what distinguishes a recombinant from a merely divergent isolate, and
    recombination is how begomoviruses acquire resistance-breaking traits wholesale
    rather than one mutation at a time.
    """
    if not (len(query) == len(parent_a) == len(parent_b)):
        raise ValueError("sequences must be aligned and the same length")

    windows = []
    for start in range(0, max(1, len(query) - window + 1), step):
        q = query[start : start + window]
        a = parent_a[start : start + window]
        b = parent_b[start : start + window]

        matches_a = sum(1 for x, y in zip(q, a) if x == y and x in "ACGT")
        matches_b = sum(1 for x, y in zip(q, b) if x == y and x in "ACGT")
        windows.append({
            "start": start,
            "similarity_a": round(matches_a / len(q), 4),
            "similarity_b": round(matches_b / len(q), 4),
            "closer_to": "a" if matches_a > matches_b else "b" if matches_b > matches_a else "tie",
        })

    parents = [w["closer_to"] for w in windows if w["closer_to"] != "tie"]
    switches = sum(1 for x, y in zip(parents, parents[1:]) if x != y)

    breakpoints = [
        windows[i + 1]["start"]
        for i in range(len(windows) - 1)
        if windows[i]["closer_to"] != windows[i + 1]["closer_to"]
        and "tie" not in (windows[i]["closer_to"], windows[i + 1]["closer_to"])
    ]

    return {
        "windows": windows,
        "switches": switches,
        "breakpoints": breakpoints,
        # One switch is a breakpoint. Many switches is noise, not ten recombinations.
        "recombinant": 1 <= switches <= 3,
    }
