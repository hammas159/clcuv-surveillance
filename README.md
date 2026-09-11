# clcuv-surveillance

[![ci](https://github.com/hammas159/clcuv-surveillance/actions/workflows/ci.yml/badge.svg)](https://github.com/hammas159/clcuv-surveillance/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.10%2B-blue)
![dependencies](https://img.shields.io/badge/dependencies-none-success)
![license](https://img.shields.io/badge/license-MIT-green)

**Genomic surveillance for Cotton Leaf Curl Virus — which variant is winning, why, and
whether our diagnostics can still see it.**

Mutation atlas · selection pressure · phylogeny · novel-strain detection ·
recombination. Zero dependencies.

---

## The problem

Cotton Leaf Curl Virus is the most serious threat to Punjab's cotton crop.
Resistance-breaking strains have already defeated resistant cultivars **and** the PCR
primers used to detect them — a field can test clean while the crop fails.

Surveillance has to answer four questions, and each one has a way of being answered
wrongly that is worse than not answering it at all.

## 1. Which variant is rising — and is that real?

A mutation sitting at 40% for three seasons is background. One that went 2% → 8% → 25%
is a lineage winning, and worth knowing about **while it is still at 25%**. So the atlas
is built around trajectories, and alerting is on the derivative rather than the level.

**But effect size alone is not evidence.** Two seasons of 100 genomes can differ by ten
points with nothing happening at all. That was a real false positive here: a variant
held at a constant 40% was reported as rising because two samples happened to land at
33% and 45%. A two-proportion z-test now runs alongside the effect-size threshold, and
the noise case is a test:

```python
assert emerging_variants(variants) == []        # z below 1.96
assert emerging_variants(variants, min_z=0.0)   # effect size alone would have fired
```

Periods with too few sequences are excluded outright — a variant at "100%" in a period
where three genomes were sequenced is not at 100%, it is unmeasured. And fold change
from zero is reported as `None`, because *"newly detected"* and *"exploding"* are
different claims.

## 2. Is something selecting for it?

Every virus is always changing. The question is whether anything is **rewarding** the
change.

```
dN/dS < 1    purifying — the protein is constrained
dN/dS ≈ 1    neutral drift
dN/dS > 1    positive selection — something is rewarding change
```

For a plant virus under host-resistance pressure, dN/dS > 1 in the coat protein is the
signature of a strain learning to break resistance. It is **invisible if you only count
mutations**: a region can mutate fast and be under strong purifying selection at the
same time.

Nei-Gojobori counting, implemented directly. When there are no synonymous differences
the ratio is `None` — reporting infinity would turn *"we cannot tell"* into *"strong
positive selection"*, which is the wrong direction to be wrong in for an alert.

## 3. How are the strains related?

**Neighbour-joining, not UPGMA — and the difference is the point.**

UPGMA assumes a molecular clock: every lineage evolving at the same rate. Real viral
lineages do not, and the one under resistance pressure evolves *fastest* — which is
precisely the lineage surveillance cares about. UPGMA places it wrongly, and does so
confidently.

Both are implemented, and the comparison is a test on an exact additive matrix where
the true tree is known:

```python
def test_neighbour_joining_recovers_the_true_branch_lengths():
    assert lengths == pytest.approx({"A": 5.0, "B": 2.0, "C": 2.0, "D": 3.0})

def test_upgma_gets_it_wrong_when_rates_differ():
    assert not TRUTH <= clades(upgma(NAMES, MATRIX))

def test_upgma_is_correct_when_the_clock_holds():
    ...   # it is not a broken method, it is a method with an assumption
```

Distances use **Jukes-Cantor**, because counting differences undercounts divergence —
the same site can mutate twice and the second change hides the first. At 5% divergence
that hardly matters; begomoviruses routinely reach 30%. At p ≥ 0.75 the distance is
`inf`, not a large number: the sequences are statistically indistinguishable from
random and the honest statement is that it cannot be estimated.

## 4. Is this something we have never seen?

**A classifier that must choose will assign a novel recombinant to whichever strain it
resembles least-badly — and report it with the same confidence as a real match.** For
surveillance that is the worst possible failure, because the novel strain is the one the
system exists to find.

So classification returns `None` above a distance threshold, calibrated from each
strain's own internal spread. An unassigned genome is a **finding**, not a gap:

```python
{"strain": None, "novel": True,
 "reason": "distance 0.8602 exceeds the threshold 0.0500 for its nearest strain Burewala"}
```

A small margin between the top two strains is also flagged — a genome sitting between
two lineages is not a tie to be broken, it is the answer.

k-mer profiles rather than alignment, because a recombinant whose segments come from two
parents has **no single good alignment** — exactly the case where an alignment-based
classifier degrades silently.

### Recombination is detected as a switch, not a distance

A point mutation shifts similarity to both parents a little. Recombination **moves a
block**: similarity to one parent rises while similarity to the other falls, at a
specific coordinate. Detecting the switch is what separates a recombinant from a merely
divergent isolate — and recombination is how begomoviruses acquire resistance-breaking
traits wholesale rather than one mutation at a time.

## It closes the loop with primer-designer

The mutation atlas feeds
[`primer-designer`](https://github.com/hammas159/primer-designer) directly:

> **This one finds the virus escaping. That one tells you whether your test can still
> see it, and designs a replacement.**

An emerging variant under a primer binding site is a primer-health alert waiting to
happen, and a 3′-terminal mismatch there means the assay goes blind while still
reporting cleanly.

## Tests

**61 tests. No dependencies, no sequence downloads, no BLAST.**

Phylogenetics and selection are exact — an additive matrix has one correct tree, and
dN/dS on synonymous-only changes is zero — so those are asserted rather than
approximated.

| Covered | |
|---|---|
| Genetic code | all 64 codons, translation, stop handling, partial codons dropped not padded |
| Selection | Met/Trp have no synonymous sites, fourfold degeneracy, synonymous-only = 0, undetermined ≠ infinite, multi-substitution averaging |
| Distance | gaps as missing data, Jukes-Cantor magnitude, saturation to `inf`, matrix symmetry |
| Trees | **NJ recovers exact topology and branch lengths**, **UPGMA fails off-clock**, UPGMA correct on-clock, Newick output |
| Atlas | variant naming, singleton removal, rising detected, **sampling noise rejected**, small periods excluded, fold from zero |
| Geography | spread vs concentration, Herfindahl index |
| Classification | known assigned, **novel refused**, unfitted refuses, length-invariant distance |
| Recombination | block swap and breakpoint located, non-recombinant clean |

## Limits

- **No aligner.** Sequences must arrive aligned; MAFFT or MUSCLE does that job and this
  repo does not duplicate it.
- Jukes-Cantor assumes equal base frequencies and equal substitution rates. Kimura
  two-parameter or GTR fit real data better; JC is the one you can read in four lines.
- Neighbour-joining is distance-based. Maximum likelihood and Bayesian inference are
  more accurate and need a substitution model, an optimiser, and far more compute.
- Recombination detection compares against **two named parents**. Screening all pairs in
  a population is the real task, and it is combinatorially larger.
- dN/dS is pairwise Nei-Gojobori. Site-specific and branch-specific selection need a
  phylogeny and a codon model.
- **No real data ships with this repo.** CLCuV genomes are public in NCBI Virus; the
  analysis is here, the sequences are not.

## License

MIT
