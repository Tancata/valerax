# Whole-genome duplication (WGD) inference in AleRax

This is a working report on an extension to [AleRax](https://github.com/BenoitMorel/AleRax)
that adds **whole-genome duplication (WGD) modelling** to its undated
gene-tree/species-tree reconciliation likelihood, with a per-branch retention
parameter `q`. The goal is to let AleRax test hypothetical WGDs (and estimate
how many duplicated copies were retained) the way
[WHALE](https://github.com/arzwa/Whale.jl) (Zwaenepoel & Van de Peer 2019) does,
while reusing AleRax's fast amalgamated-likelihood machinery and extending it to
the duplication–transfer–loss (DTL) model.

> **Status:** research prototype. The model is implemented and validated against
> WHALE on a small published dataset. It is a *scaffold* — see *Limitations*.

---

## 1. What was implemented

### The model
At a WGD, a surviving lineage instantaneously **doubles**; the duplicate is
**retained with probability `q`** and lost with probability `1 − q` (the
asymmetric "1-or-2" retention model used by WHALE / Rabier et al. 2014).
Crucially, **`q = 0` recovers the ordinary no-WGD model**, which is the null
hypothesis a likelihood-ratio test rests on.

This gives two closed-form transforms applied at the top of a species branch
(`E` = extinction probability of the branch, `P` = clade likelihood at the
branch, `Σ_splits` = sum over conditional-clade splits):

```
E_top = (1 − q)·E + q·E²
P_top = (1 − q)·P + q·[ 2·E·P + Σ_splits ]
```

### No tree surgery
Rather than inserting a degree-2 node into the species tree (which fights
corax's node arrays), the WGD is attached to the **top of an existing branch**
and the transform is applied only to the values a branch's **parent** consumes.
The rule that keeps this correct:

- reads of a **child** branch (during speciation / speciation-loss) go through
  the transformed `…Top` arrays;
- reads of a branch's **own** events (duplication, loss, transfer, origination)
  stay raw, because the WGD sits *above* the branch while those events happen
  *within* it.

With no WGD declared, every `…Top` array is an exact copy of its raw
counterpart, so the likelihood is **bit-for-bit identical to stock AleRax**
(this is the primary regression test).

### Where it lives
- `src/ale/UndatedDLMultiModel.hpp` — DL model (duplication + loss).
- `src/ale/UndatedDTLMultiModel.hpp` — DTL model (adds transfers). Here a
  transferred lineage is taken to land *within* the destination branch, below
  any WGD there, so transfer/highway reads stay raw; only the vertical
  parent→child reads use the transformed arrays.
- `src/ale/MultiModel.hpp`, `AleEvaluator.*`, `AleArguments.*`, `ale.cpp` —
  plumbing: a `setWGD` broadcast, the CLI option, and a retention optimizer.

### Using it
```sh
alerax -f families.txt -s species_tree.nw --rec-model UndatedDL \
       --gene-tree-rooting UNIFORM --species-tree-search SKIP \
       --wgd PPAT --wgd ATHA,ATRI -p output_dir
```
`--wgd LABEL` places a WGD on a terminal branch; `--wgd LABEL1,LABEL2` places it
on the branch leading to the LCA of two taxa; an optional `:q0` sets the
starting retention. The retention `q` of each declared WGD is then estimated by
maximum likelihood (reparameterised as `q = r/(1+r)`, `r ≥ 0`), in a dedicated
optimisation pass that runs even when D/L rates are fixed — so `q` can be
estimated against a fixed-rate background for a clean likelihood-ratio test.

### Tests
`tests/wgd_regression.cpp` (CMake target `wgd_regression`) checks:
1. **no-WGD == stock AleRax** likelihood, to the bit
   (`−58.689981412938` on the bundled `simulated_2` data), and
2. declaring a WGD and sweeping `q → 0` converges monotonically back to that
   value (confirming `q = 0` is the no-WGD null).

---

## 2. Benchmarking against WHALE

We validated against WHALE on its own **`example-1`** dataset (12 land-plant
gene families over 9 species; the inputs are staged here under
`validation/example-1/`). Both tools use the **same gene→species mapping** (the
species code is the label prefix before the first `_`) and the **same two
hypothetical WGDs** WHALE's documentation uses:

- **PPAT** — the *Physcomitrella* (moss) terminal branch;
- **ATHA–ATRI LCA** — the angiosperm crown branch.

As a sanity check, the WHALE install reproduces its documented reference
likelihoods exactly (e.g. `−592.0186` summed over the 12 families), so the
comparison is against a faithful WHALE.

**Maximum-likelihood retention estimates** (WHALE: constant-rate DLWGD, η fixed;
AleRax: UndatedDL). p-values use the boundary-correct ½χ²₀+½χ²₁ mixture, because
`q = 0` sits on the edge of the parameter space:

| tool | WGD branch | q̂ | 2·ΔLL | p | detect @0.05 |
|---|---|---|---|---|---|
| **AleRax** | PPAT | 0.31 | 3.31 | 0.034 | **yes** |
| **WHALE** | PPAT | 0.19 | 0.51 | 0.24 | no |
| **AleRax** | angiosperm | 0.00 | 0.00 | 1.0 | no |
| **WHALE** | angiosperm | 0.00 | 0.00 | 1.0 | no |

**What agrees:** both tools put the retention signal on the **same branch**
(PPAT), and both drive the **angiosperm `q̂` to exactly 0** — i.e. they agree
there is no evidence for the deeper WGD in this tiny dataset. The agreement on
the null is exact.

**What differs:** AleRax estimates a *stronger* moss WGD (higher `q̂`, ~6× larger
LRT) and would call it significant, where WHALE would not. This is expected and
not a bug — WHALE is a **dated** model (it uses the branch lengths and slices
time), whereas AleRax `UndatedDL` ignores branch lengths; the two also use
different root/origination priors and place the WGD at slightly different points
on the branch. WHALE's `q̂` is stable across its time discretisation, so the gap
is a genuine model difference, not a numerical artifact. WHALE's own docs treat
`example-1` as illustrative and caution against over-reading an LRT on so little
data — good context for both tools' verdicts.

**Bottom line:** the implementation is validated on the null (no false WGD) and
agrees with WHALE on *where* the signal is; it is simply more liberal about
*how strong* it is, consistent with being an undated model.

---

## 3. Effect of allowing transfers (DL vs DTL)

A natural worry: in the DTL model, could **horizontal transfers explain away** a
WGD — i.e. account for the extra gene copies that would otherwise be evidence
for genome doubling? We ported the WGD transform into the DTL model and re-ran
the same analysis (transfers with the realistic `PARENTS` constraint).

| model | null LL | +PPAT LL | q̂[PPAT] | 2·ΔLL | p | detect |
|---|---|---|---|---|---|---|
| **DL** (no transfers) | −362.4 | −360.7 | 0.31 | 3.31 | 0.034 | yes |
| **DTL** (transfers) | −330.2 | −327.9 | **0.34** | **4.71** | **0.015** | yes |

Two observations:

1. **DTL fits the data much better overall** (null −330 vs −362): transfers
   absorb ~32 log-units of gene-tree structure that the DL model otherwise
   charges to duplication and loss.
2. **The WGD signal survives — and slightly strengthens** (`q̂` 0.31 → 0.34,
   LRT 3.3 → 4.7). Allowing transfers did **not** explain away the moss WGD.

**Why** transfers can't absorb it: the duplicated *Physcomitrella* genes form
**moss–moss sister pairs** in the gene trees. A transfer would make a moss gene
sister to *another species'* gene, not to a second moss gene, so the ohnolog
topology is only explainable by duplication/WGD. Transfers soak up the *other*
incongruence in the dataset (hence the much better baseline fit) but leave the
ohnolog pairs to the WGD. The same four multi-copy families drive the signal in
both models, and the angiosperm WGD stays undetected (`q̂ = 0`) either way. The
finding is robust to the transfer constraint (`NONE` gives essentially the same
`q̂`).

---

## 4. Limitations

- **Tiny dataset.** `example-1` is 12 families; significance values here are
  illustrative, not conclusions about plant evolution.
- **Undated model.** Unlike WHALE, AleRax `UndatedDL/DTL` ignore branch lengths,
  which is why absolute `q̂`/LRT differ from WHALE. Only relative comparisons
  (DL vs DTL, branch vs branch) should be trusted across the two tools.
- **No synteny.** From gene trees alone, a WGD is **not identifiable** against a
  burst of small-scale duplications on the same branch. This DL/DTL+WGD model is
  a scaffold; distinguishing the two would need a synteny-block coupling layer
  (a planned future direction).
- The DTL port currently has an empirical regression check (no-WGD likelihood
  unchanged); a dedicated C++ unit test like the DL model's is still to add.

---

## 5. Reproducing

```sh
# build (no MPI):
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DDISABLE_MPI=ON
cmake --build . --target alerax wgd_regression -j4

# regression test (no-WGD == stock, and q->0 convergence):
./bin/wgd_regression <species_tree> <gene_trees> <mapping>

# WHALE benchmark inputs and the exact run recipe:
#   validation/example-1/README.md
```

---

*This extension is a prototype for research use. The retention estimates should
be interpreted alongside the limitations above, and cross-checked against a
dated and/or synteny-aware method (e.g. WHALE) before drawing biological
conclusions.*
