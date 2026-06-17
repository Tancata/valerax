# DTL + LORe: lineage-specific rediploidization with transfers

This documents the extension of the **LORe** (lineage-specific ordered
rediploidization) model from the DL model to the **DTL** (duplication–transfer–
loss) model in `kalerax`. Before this, `--lore` was DL-only; it now works with
`--rec-model UndatedDTL` as well.

## What LORe models (recap)

After a WGD the two ohnolog copies need not diverge immediately. A gene lineage
carries a hidden state — **R** (resolved/disomic) or **U** (unresolved/
tetrasomic: a single lineage still holding a not-yet-diverged ohnolog pair). At
the WGD branch the retained duplicate enters **U** with probability `q`; a U
lineage passes through speciation **as a unit** and, on each branch, **resolves
with probability `r`** (resolution = ohnolog divergence, recorded on the branch
where it happens).

- `q = 0` → no WGD (the null);
- `r = 1` → resolution certain at the WGD branch → **AORe** (immediate, WHALE-like);
- `r < 1` → resolution can be **delayed** into descendant lineages → **LORe**.

`r` is a single global parameter, fitted jointly with `q` by maximum likelihood.

## How the DTL port works

The DL LORe machinery is mirrored into `UndatedDTLMultiModel`:

- a per-clade unresolved CLV `_uu` (parallel to the resolved `_uq`), an
  unresolved extinction `_uEU`, and a per-branch resolution probability
  `_resolutionProbs`;
- the **WGD transform** the parent consumes becomes the R/U mixture
  `P_top = (1−q)·P + q·U` (and `E_top = (1−q)·E + q·EU`), where the doubling
  bracket `U`/`EU` is the resolved value when `r = 1`;
- a **U-aware backtrace** (`sampleResolutionCommits`) that reads out, as a
  posterior marginal under the fitted `r`, the species branch on which each
  ohnolog pair resolved. An always-on consistency guard (`check=true`) asserts
  that every sampled-term weight equals the exact inside-CLV (`_uu`, `_uqTop`).

### The one modelling choice: U × transfers

The U-state propagates **vertically only** — through speciation (both daughters
inherit U as a unit), speciation-with-unresolved-loss, and the duplication
bracket. **Transfers and highways act on resolved lineages only**: a tetrasomic
(U) lineage does not transfer, and a transferred copy is born **resolved**,
landing within the destination branch *below* any WGD there. This is the natural
extension of the existing AORe-DTL convention ("a transferred lineage lands below
the WGD, so transfer reads stay raw") and keeps the U recursion transfer-free.
Consequently `_uEU` has no transfer-extinction term, and the backtrace continues
a transferred copy with `btR` (resolved) at the destination.

This is a deliberate first-cut assumption, not a derived necessity; a model in
which a tetrasomic locus can itself transfer (carrying U) would be a strictly
richer — and much harder to identify — alternative.

## Reduction to AORe (the regression gate)

At `r = 1` every U term collapses to the resolved doubling bracket, so
`P_top` and `E_top` become exactly the AORe/WHALE transforms. The model is
therefore **bit-for-bit identical to the existing DTL+WGD (AORe) likelihood**
when `r = 1`, which is the nested null an LRT on `r` rests on.

**Validation performed:**

1. **Bit-for-bit r=1 == AORe.** On a 30-family salmonid set with fixed rates, the
   post-port DTL null and DTL+WGD joint log-likelihoods reproduce the pre-port
   (baseline) values exactly: `−1863.892400` (no WGD) and `−1854.694800`
   (+WGD, q̂=0.36303). Adding the U-state machinery does not perturb the AORe /
   no-WGD likelihood.
2. **Backtrace consistency.** With the `check=true` guard enabled, the U-cell and
   WGD R/U-coin sampling weights equal the inside CLVs at every cell across many
   families × samples (no assertion fired).
3. **Resolution profile localisation.** All inferred ohnolog-resolution events
   fall within the WGD clade; every branch outside it reads exactly 0.
4. **Bundled regression gate.** `tests/dtl_lore_regression.sh` runs DTL null,
   DTL+WGD (AORe), and DTL+WGD+`--lore` on `validation/example-1` and asserts the
   nesting `null ≤ AORe ≤ LORe` with `r̂ ∈ [0,1]`. On the moss (PPAT) terminal
   branch the WGD is detected and LORe correctly reverts to `r = 1` (a terminal
   branch has no descendant lineages, so there is nothing to delay).

## Using it

```sh
kalerax -f families.txt -s species_tree.nw --rec-model UndatedDTL \
  --gene-tree-rooting UNIFORM --species-tree-search SKIP \
  --wgd TAXONA,TAXONB --lore -p out
#   -> out/reconciliations/{totalSpeciesResolutionCounts,wgdSummary}.txt
```

## Caveats

- Inherits all the LORe caveats (identifiability lives in gene-tree topology;
  the resolution *branch* is a posterior readout under one global `r`, not a
  dated time; undated model infers *where*, not *when*).
- The U × transfer interaction is the assumption above; results should be read
  with that in mind. As a robustness check, compare DL+LORe and DTL+LORe: if the
  resolution-timing profile is stable across the two, transfers are not driving
  the LORe signal.
- Like the rest of the WGD/LORe code, this is a research prototype.
