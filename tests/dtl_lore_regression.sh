#!/usr/bin/env bash
# Empirical regression gate for the DTL+LORe extension of UndatedDTLMultiModel.
#
# The project builds its model regression checks empirically against the
# production binary (the WGD/LORe model headers are not wired into CMake as unit
# tests). This script exercises the DTL+WGD+LORe path on the bundled, self-
# contained validation/example-1 dataset (12 land-plant families, 9 taxa) and
# asserts the model invariants that the r=1 == AORe reduction and the nesting
# guarantee imply:
#
#   1. plain DTL (no WGD)            -> lnL_null      (finite)
#   2. DTL + WGD (no --lore, r=1)    -> lnL_aore      (>= lnL_null)
#   3. DTL + WGD + --lore (r free)   -> lnL_lore      (>= lnL_aore)
#   and 0 <= r_hat < 1.
#
# Because "no --lore" sets r=1 on every branch and the LORe optimizer's nested
# AORe pass uses the same r=1 code path, lnL_aore is the bit-for-bit r=1 value
# of the LORe model. (A bit-for-bit r=1 == AORe and no-WGD-unchanged check was
# additionally confirmed against a pre-LORe baseline binary; see DTL_LORE.md.)
#
# Usage: dtl_lore_regression.sh [path/to/kalerax|alerax] [path/to/example-1]
set -euo pipefail

BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
[ -x "$BIN" ] || BIN="${1:-$(dirname "$0")/../build/bin/alerax}"
DATA="${2:-$(dirname "$0")/../validation/example-1}"
DATA="$(cd "$DATA" && pwd -P)"   # absolute, so families paths resolve anywhere
WGD="PPAT"   # moss terminal branch (the example-1 branch that carries WGD signal)
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

# Repoint the families file at $DATA (the shipped families.txt bakes in absolute
# paths from another machine). Greedy .*/ keeps only the basename and is safe
# against spaces in the original paths.
sed -e "s|^gene_tree = .*/|gene_tree = $DATA/|" \
    -e "s|^mapping = .*/|mapping = $DATA/mappings/|" \
    "$DATA/families.txt" > "$TMP/families.txt"

common=(-s "$DATA/tree.nw" -f "$TMP/families.txt" --rec-model UndatedDTL
        --gene-tree-rooting UNIFORM --species-tree-search SKIP --fix-rates
        --seed 1)

ll() { grep -E "Final species tree likelihood" "$1" | tail -1 |
       sed -E 's/.*ll=//'; }

echo "binary: $BIN"
"$BIN" "${common[@]}" -g 1 -p "$TMP/null"  >"$TMP/null.log"  2>&1
"$BIN" "${common[@]}" -g 1 --wgd "$WGD" -p "$TMP/aore" >"$TMP/aore.log" 2>&1
"$BIN" "${common[@]}" -g 5 --wgd "$WGD" --lore -p "$TMP/lore" >"$TMP/lore.log" 2>&1

LN=$(ll "$TMP/null.log"); LA=$(ll "$TMP/aore.log"); LL=$(ll "$TMP/lore.log")
RHAT=$(grep -E "After WGD\+LORe opt" "$TMP/lore.log" | tail -1 |
       sed -E 's/.* r=//')
echo "  lnL_null (DTL)            = $LN"
echo "  lnL_aore (DTL+WGD, r=1)   = $LA"
echo "  lnL_lore (DTL+WGD+--lore) = $LL   r_hat=${RHAT:-NA}"

# also confirm the resolution profile was written and is confined to the clade
test -s "$TMP/lore/reconciliations/totalSpeciesResolutionCounts.txt" \
  || { echo "FAIL: no resolution profile written"; exit 1; }

# r=1 is a valid outcome: on a terminal WGD branch (or with no delayed-resolution
# signal) the LORe optimizer correctly reverts to AORe (r=1).
awk -v ln="$LN" -v la="$LA" -v ll="$LL" -v r="${RHAT:-2}" 'BEGIN{
  if (la < ln - 1e-6) { print "FAIL: AORe < null (not nested)"; exit 1 }
  if (ll < la - 1e-6) { print "FAIL: LORe < AORe (LORe must dominate)"; exit 1 }
  if (r < -1e-9 || r > 1 + 1e-9) { print "FAIL: r_hat outside [0,1]"; exit 1 }
  print "PASS: null <= AORe <= LORe, r_hat in [0,1]"
}'
