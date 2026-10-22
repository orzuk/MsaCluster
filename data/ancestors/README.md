# Reconstructed XCL1 (lymphotactin) ancestors — Dishman et al. 2021

## Source

> Dishman, A. F., Tyler, R. C., Fox, J. C., Kleist, A. B., Prehoda, K. E.,
> Babu, M. M., Peterson, F. C., Volkman, B. F.
> **Evolution of fold switching in a metamorphic protein.**
> *Science* **371**, 86–90 (2021). doi:10.1126/science.abd8700

Sequences transcribed from the supplementary materials, **fig. S1C**
("Ancestral reconstruction of the XC chemokine family", panel C — the
72-column alignment of the reconstructed ancestors against human XCL1).

## Files

| File | Contents |
|---|---|
| `dishman_ancestors.fasta` | the five sequences, ungapped |
| `dishman_ancestors_aligned.fasta` | the same five rows in the 72-column fig. S1C alignment (`-` = gap) |

Record order is evolutionary (root → tip): `Anc.0`, `Anc.2`, `Anc.3`,
`Anc.4`, `XCL1`. `Anc.4` is the direct ancestor of the modern human tip
`XCL1`. Anc.1 was not reconstructed by the authors (it carries two
disulfides and was judged unlikely to be metamorphic).

## Experimentally measured fold populations (the calibration ground truth)

NMR, 37 °C / 150 mM NaCl unless noted. "Alternative fold" = the all-beta
XCL1 dimer fold (PDB 2N54); "chemokine fold" = the canonical monomeric
chemokine fold (PDB 2HDM / 1J8I).

| Sequence | Chemokine fold | Alternative fold | Note |
|---|---|---|---|
| Anc.0 | 100 % | 0 % | monomorphic; NMR structure PDB **7JH1** |
| Anc.2 | 100 % | 0 % | single fold at all tested conditions |
| Anc.3 | 92 ± 3.7 % | ~8 % | metamorphic |
| Anc.4 | 9.3 ± 2.0 % | ~91 % | metamorphic |
| XCL1 | ~50 % | ~50 % | metamorphic (modern human tip) |

So the expected ordering of preference for the **alternative** fold is

    Anc.0 ~ Anc.2  <  Anc.3  <  XCL1  <  Anc.4

This table is also hard-coded in `scripts/ancestor_calibration.py`
(`GROUND_TRUTH`), which is the single place the numbers are consumed.

## Verification of the transcription

- **Anc.0** is an exact, residue-for-residue match (64/64) to the RCSB
  sequence of **PDB 7JH1** (`>7JH1_1|Chain A|XCL1 ancestor`), the NMR
  structure the paper deposited for Anc.0.
- The **XCL1** row (69 aa) is an exact prefix of the UniProt **P47992**
  mature chain (residues 22–114), i.e. mature residues 1–69, with no
  substitutions — wild-type sequence, not the CC3 / W55D disulfide-locked
  variant.
- The percent identities to XCL1 implied by the alignment (41.0 / 48.4 /
  60.6 / 68.2 % for Anc.0/2/3/4, identities over columns occupied by both
  rows) reproduce the values printed in the figure (40 / 49 / 60 / 68 %)
  to within ~1 %.
- Cysteine counts: Anc.0 has 4 (two disulfides); Anc.3, Anc.4 and XCL1
  have 2 (one disulfide). **Anc.2 has 3** — the 1–3 disulfide is lost by
  a single C→V substitution at alignment column 10, leaving the former
  Cys3 (column 35) free and unpaired. That is what fig. S1C shows.

## Reference folds used by the calibration

The fold-switch pair `2n54B_2hdmA` in this repository:

- **F1 = 2N54 chain B** — the alternative (all-beta dimer) XCL1 fold.
- **F2 = 2HDM chain A** — the canonical chemokine fold.

Per the repo-wide convention, `Delta = score(F1) - score(F2)`, so
**positive Delta = prefers the alternative fold**.
