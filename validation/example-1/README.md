# WGD cross-validation — Whale.jl example-1 (STEP 4, Tier A)

Inputs for cross-validating the AleRax WGD extension against WHALE
(Zwaenepoel & Van de Peer 2019) on the `Whale.jl` `example/example-1` dataset.

> The `*.ale` files and `tree.nw` are redistributed from `arzwa/Whale.jl`
> (LGPL-2.1). See [`SOURCE.md`](SOURCE.md) for attribution and citation.

## Contents

- `*.ale` — 12 ALE-format CCP files (gene-tree distributions), copied from
  `arzwa/Whale.jl` `example/example-1/ale/`.
- `tree.nw` — the 9-taxon species tree (copied from the same example).
- `make_inputs.py` — regenerates `mappings/` and `families.txt` (below).
- `mappings/<family>.link` — per-family gene→species maps, "treerecs" format
  (`<gene> <species>`), one line per gene leaf.
- `families.txt` — AleRax families file (each family: `gene_tree = <.ale>`,
  `mapping = <.link>`).

Species tree taxa: MPOL, PPAT, SMOE, OSAT, ATHA, CPAP, ATRI, GBIL, PABI.

## Mapping rule

Species code = the gene-label prefix before the first `_`
(e.g. `ATHA_AT5G48120.1_AT5G48120` → `ATHA`, `GBIL_Gb_13638` → `GBIL`). This is
the same rule AleRax uses for automatic mapping, written out explicitly here so
the mapping is auditable and identical to what WHALE uses.

## CCP rooting

The `.ale` files are unrooted CCPs. AleRax detects the leading `#` and parses
them as ALE format, which requires `--gene-tree-rooting UNIFORM` (all root
positions equally weighted = unrooted). Do not use ROOTED/MAD here.

## WGD placement

Two WGDs (matching the WHALE example-1 setup):
- `--wgd PPAT` — top of the PPAT (moss) terminal branch.
- `--wgd ATHA,ATRI` — top of the branch leading to the LCA of ATHA and ATRI,
  i.e. the angiosperm crown clade {OSAT, ATHA, CPAP, ATRI}
  (resolves to the internal node auto-labeled `Node_OSAT_ATRI_0`).

## Running

NOTE: the repository path contains a space ("My Drive"), which AleRax's input
reader does not handle. Stage to a space-free directory first and regenerate
`families.txt` with absolute paths into that directory:

```sh
cp -R validation/example-1 /tmp/wgd_example1
cd /tmp/wgd_example1
python3 make_inputs.py --root /tmp/wgd_example1   # rewrites families.txt + mappings
```

Then run kalerax (species tree fixed; WGD node indices are tied to the topology):

```sh
<repo>/build/bin/kalerax \
  -f /tmp/wgd_example1/families.txt \
  -s /tmp/wgd_example1/tree.nw \
  --rec-model UndatedDL \
  --gene-tree-rooting UNIFORM \
  --species-tree-search SKIP \
  --wgd PPAT \
  --wgd ATHA,ATRI \
  -p /tmp/wgd_example1/out
```

To estimate q with D/L fixed (for the WHALE q-grid / LRT comparison) add
`--fix-rates --d <D> --l <L>`; the retention optimizer runs regardless of
`--fix-rates`. For the q=0 nested null (LRT), the retention optimizer drives
r→0 ⇒ q→0 on a branch with no signal; compute `2*(LL(q-hat) - LL(0))` and use
the 50:50 χ²₀/χ²₁ boundary mixture (halve the naive χ²₁ p-value).

## Reference (smoke check, default DL rates)

A verification run of the command above produced (per-branch retentions):
`q[PPAT]≈0.31`, `q[angiosperm]≈1e-10`, joint `ll≈-360.72` over the 12 families.
These are just a wiring check, not the validated estimates — STEP 4 drives the
WHALE comparison and the LRT.
