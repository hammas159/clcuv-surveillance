"""Codons, translation, and selection pressure.

The question surveillance actually has to answer is not *"is this virus changing?"* —
every virus is always changing — but **"is something selecting for the change?"**

The genetic code is redundant: most amino acids have several codons, so many
substitutions change the DNA and not the protein. Those are *synonymous*, and they
accumulate at roughly the neutral mutation rate. Substitutions that do change the
protein are *non-synonymous*, and selection acts on them.

    dN/dS < 1    purifying selection — the protein is constrained, changes get removed
    dN/dS ≈ 1    neutral drift — nothing is selecting
    dN/dS > 1    positive selection — something is rewarding change

For a plant virus under host-resistance pressure, dN/dS > 1 in the coat protein is the
signature of a strain learning to break resistance. That is the signal worth alerting
on — and it is invisible if you only count mutations, because a region can be mutating
fast and be under strong purifying selection at the same time.

Nei-Gojobori counting, implemented directly. It is a loop over codons and a table.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

BASES = "TCAG"

# Standard genetic code, ordered to match the TCAG index expansion below.
_AMINO = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

CODON_TABLE: dict[str, str] = {
    a + b + c: _AMINO[i * 16 + j * 4 + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}

STOP = "*"


class CodonError(ValueError):
    pass


def translate(sequence: str, *, to_stop: bool = False) -> str:
    """Translate a nucleotide sequence in frame 1.

    A trailing partial codon is dropped rather than padded: padding invents a residue
    that is not in the data, and one fabricated amino acid at the end of every protein
    is exactly the kind of small lie that survives into a published figure.
    """
    seq = sequence.upper().replace("U", "T")
    protein = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        residue = CODON_TABLE.get(codon, "X")  # ambiguous or gapped codon
        if to_stop and residue == STOP:
            break
        protein.append(residue)
    return "".join(protein)


def is_synonymous(codon_a: str, codon_b: str) -> bool:
    a = CODON_TABLE.get(codon_a.upper())
    b = CODON_TABLE.get(codon_b.upper())
    if a is None or b is None:
        return False
    return a == b


def synonymous_sites(codon: str) -> float:
    """Expected synonymous sites in a codon (Nei-Gojobori).

    Each of the three positions contributes the fraction of its three possible
    substitutions that leave the amino acid unchanged. A codon therefore has a
    fractional number of synonymous sites — which is the point: the denominator has to
    account for the fact that third positions are mostly free and second positions
    almost never are.
    """
    codon = codon.upper()
    if codon not in CODON_TABLE:
        return 0.0

    total = 0.0
    for position in range(3):
        synonymous = 0
        for base in "ACGT":
            if base == codon[position]:
                continue
            mutant = codon[:position] + base + codon[position + 1 :]
            if is_synonymous(codon, mutant):
                synonymous += 1
        total += synonymous / 3
    return total


@dataclass(frozen=True)
class Selection:
    synonymous_differences: float
    nonsynonymous_differences: float
    synonymous_sites: float
    nonsynonymous_sites: float
    codons_compared: int

    @property
    def ps(self) -> float:
        return (
            self.synonymous_differences / self.synonymous_sites
            if self.synonymous_sites else 0.0
        )

    @property
    def pn(self) -> float:
        return (
            self.nonsynonymous_differences / self.nonsynonymous_sites
            if self.nonsynonymous_sites else 0.0
        )

    @property
    def dn_ds(self) -> float | None:
        """None when there are no synonymous differences.

        Reporting infinity, or silently substituting a large number, turns "we cannot
        tell" into "strong positive selection" — which is the wrong direction to be
        wrong in for an alerting system.
        """
        if self.ps == 0:
            return None
        return round(self.pn / self.ps, 6)

    def interpretation(self) -> str:
        ratio = self.dn_ds
        if ratio is None:
            return "undetermined - no synonymous differences observed"
        if ratio < 0.5:
            return "purifying selection - the protein is constrained"
        if ratio < 1.0:
            return "weak purifying selection"
        if ratio <= 1.2:
            return "approximately neutral"
        return "positive selection - something is rewarding change"


def selection_pressure(seq_a: str, seq_b: str) -> Selection:
    """Nei-Gojobori dN/dS between two aligned coding sequences."""
    a = seq_a.upper().replace("U", "T")
    b = seq_b.upper().replace("U", "T")
    if len(a) != len(b):
        raise CodonError("sequences must be aligned and the same length")

    syn_diff = non_syn_diff = 0.0
    syn_sites = non_syn_sites = 0.0
    compared = 0

    for i in range(0, min(len(a), len(b)) - min(len(a), len(b)) % 3, 3):
        codon_a, codon_b = a[i : i + 3], b[i : i + 3]
        if codon_a not in CODON_TABLE or codon_b not in CODON_TABLE:
            continue  # gapped or ambiguous codons contribute nothing

        compared += 1
        sites_a = synonymous_sites(codon_a)
        sites_b = synonymous_sites(codon_b)
        syn_sites += (sites_a + sites_b) / 2
        non_syn_sites += (3 - sites_a + 3 - sites_b) / 2

        if codon_a == codon_b:
            continue

        differences = sum(1 for x, y in zip(codon_a, codon_b) if x != y)
        if differences == 1:
            if is_synonymous(codon_a, codon_b):
                syn_diff += 1
            else:
                non_syn_diff += 1
        else:
            # Multiple substitutions in one codon have several possible mutational
            # paths. Nei-Gojobori averages over them; with no phylogeny to pick a path
            # this is the honest treatment rather than assuming an order.
            paths = _average_pathways(codon_a, codon_b)
            syn_diff += paths[0]
            non_syn_diff += paths[1]

    return Selection(
        synonymous_differences=round(syn_diff, 6),
        nonsynonymous_differences=round(non_syn_diff, 6),
        synonymous_sites=round(syn_sites, 6),
        nonsynonymous_sites=round(non_syn_sites, 6),
        codons_compared=compared,
    )


def _average_pathways(codon_a: str, codon_b: str) -> tuple[float, float]:
    """Average synonymous/non-synonymous counts over every mutational order."""
    positions = [i for i in range(3) if codon_a[i] != codon_b[i]]
    if not positions:
        return 0.0, 0.0

    import itertools

    syn_total = non_syn_total = 0.0
    orders = list(itertools.permutations(positions))

    for order in orders:
        current = codon_a
        syn = non_syn = 0
        for position in order:
            nxt = current[:position] + codon_b[position] + current[position + 1 :]
            # A path through a stop codon is not a path a lineage took.
            if CODON_TABLE.get(nxt) == STOP and CODON_TABLE.get(codon_b) != STOP:
                syn = non_syn = 0
                break
            if is_synonymous(current, nxt):
                syn += 1
            else:
                non_syn += 1
            current = nxt
        syn_total += syn
        non_syn_total += non_syn

    return syn_total / len(orders), non_syn_total / len(orders)


def amino_acid_changes(seq_a: str, seq_b: str) -> list[tuple[int, str, str]]:
    """Protein-level differences as (1-based residue, from, to).

    Reported in residue coordinates rather than nucleotide coordinates because that is
    how the literature names resistance-breaking mutations, and a surveillance report
    nobody can cross-reference is a surveillance report nobody uses.
    """
    protein_a = translate(seq_a)
    protein_b = translate(seq_b)
    return [
        (i + 1, x, y)
        for i, (x, y) in enumerate(zip(protein_a, protein_b))
        if x != y
    ]


def codon_positions(sequence: Sequence[str] | str) -> list[int]:
    """Codon position (1, 2 or 3) for each nucleotide index."""
    return [(i % 3) + 1 for i in range(len(sequence))]
