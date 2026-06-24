#!/usr/bin/env bash
# Empirical regression gate for per-event LORe r (multiple --wgd, disjoint case),
# for BOTH UndatedDL and UndatedDTL. Models the structure of
# tests/dtl_lore_regression.sh but uses a small SYNTHETIC two-clade input built
# here (no external data needed), so the simulated truth is known exactly.
#
# Synthetic setup (species tree: (OG,((A1,A2),(B1,B2)))):
#   * Clade A (stem = LCA(A1,A2)): gene trees are SPECIES-GROUPED ohnologs
#     ((A1_1,A1_2),(A2_1,A2_2)) -> the ohnolog divergence happened independently
#     in each lineage AFTER speciation => true LORe => r_A < 1.
#   * Clade B (stem = LCA(B1,B2)): gene trees are OHNOLOG-GROUPED
#     ((B1_1,B2_1),(B1_2,B2_2)) -> divergence at the WGD, before speciation
#     => true AORe => r_B ~ 1.
# The two WGD subtrees are disjoint, so per-event r is identifiable and the
# optimiser must recover two DISTINCT r-hat (A noticeably < 1, B ~ 1) and two
# independent q-hat.
#
# Invariants checked (see SPEC_multi_wgd_per_event_r.md section 5, items 5-8):
#   5. two disjoint internal WGDs -> two independent q-hat and two DISTINCT r-hat
#      matching truth, in wgdSummary.txt, for UndatedDL AND UndatedDTL.
#   6. a WGD on a terminal branch has its r pinned: log shows
#      "r[node X]=1(pinned;terminal)" and the info line counts it as pinned.
#   7. two nested *free-r* WGDs + --lore aborts with the DISJOINT message; nested
#      WGDs WITHOUT --lore run and report two independent q-hat.
#   8. the run log prints per-WGD "q[node ..]=" AND "r[node ..]=" (not one shared
#      r), and wgdSummary.txt carries the matching per-row q/r.
#   9. --lore-wgd fits r for the named WGD only and pins every other declared WGD
#      to AORe r=1, even an internal one (9a) and even when the target is nested
#      inside an AORe WGD (9b, the headline relaxation -- no abort).
#  10. --summary-only keeps the summaries/totals/wgdSummary but drops the
#      per-sample reconciliations/all/ files and the CCP binaries.
#
# Usage: multi_wgd_regression.sh [path/to/kalerax]
set -euo pipefail

BIN="${1:-$(dirname "$0")/../build/bin/kalerax}"
[ -x "$BIN" ] || BIN="${1:-$(dirname "$0")/../build/bin/alerax}"
[ -x "$BIN" ] || { echo "FAIL: binary not found/executable: $BIN"; exit 1; }
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

echo "binary: $BIN"

# --- build the synthetic two-clade input --------------------------------------
mkdir -p "$TMP/gtrees" "$TMP/maps"
SP="$TMP/species.nw"
echo '(OG:2,((A1:1,A2:1):1,(B1:1,B2:1):1):1);' > "$SP"
NA=24
NB=24
FAM="$TMP/families.txt"
echo '[FAMILIES]' > "$FAM"
for i in $(seq 1 "$NA"); do
  f="famA$i"
  # clade A doubled, SPECIES-grouped (true LORe); clade B single-copy
  echo "(OG_1,(((A1_1,A1_2),(A2_1,A2_2)),(B1_1,B2_1)));" > "$TMP/gtrees/$f.nw"
  printf "A1_1\tA1\nA1_2\tA1\nA2_1\tA2\nA2_2\tA2\nB1_1\tB1\nB2_1\tB2\nOG_1\tOG\n" \
    > "$TMP/maps/$f.link"
  printf -- "- %s\ngene_tree = %s/gtrees/%s.nw\nmapping = %s/maps/%s.link\n" \
    "$f" "$TMP" "$f" "$TMP" "$f" >> "$FAM"
done
for i in $(seq 1 "$NB"); do
  f="famB$i"
  # clade B doubled, OHNOLOG-grouped (true AORe); clade A single-copy
  echo "(OG_1,((A1_1,A2_1),((B1_1,B2_1),(B1_2,B2_2))));" > "$TMP/gtrees/$f.nw"
  printf "A1_1\tA1\nA2_1\tA2\nB1_1\tB1\nB1_2\tB1\nB2_1\tB2\nB2_2\tB2\nOG_1\tOG\n" \
    > "$TMP/maps/$f.link"
  printf -- "- %s\ngene_tree = %s/gtrees/%s.nw\nmapping = %s/maps/%s.link\n" \
    "$f" "$TMP" "$f" "$TMP" "$f" >> "$FAM"
done

common=(-f "$FAM" -s "$SP" --gene-tree-rooting UNIFORM
        --species-tree-search SKIP --fix-rates --seed 1)

# Last "After WGD+LORe opt" log line (the final, refined values).
loreLine() { grep -E "After WGD\+LORe opt" "$1" | tail -1; }

fail() { echo "FAIL: $*"; exit 1; }

# --- invariants 5 & 8: two disjoint internal WGDs, per model ------------------
for MODEL in UndatedDL UndatedDTL; do
  echo "=== $MODEL: two disjoint internal WGDs (A=LORe, B=AORe) ==="
  out="$TMP/out_${MODEL}"
  "$BIN" "${common[@]}" --rec-model "$MODEL" -g 1 \
    --wgd A1,A2 --wgd B1,B2 --lore -p "$out" > "$out.log" 2>&1 \
    || fail "$MODEL two-WGD --lore run exited non-zero"

  sum="$out/reconciliations/wgdSummary.txt"
  [ -s "$sum" ] || fail "$MODEL: no wgdSummary.txt written"

  # Two data rows (skip the header), columns: label node q r commits tetra.
  nrows=$(grep -cvE '^#' "$sum")
  [ "$nrows" -eq 2 ] || fail "$MODEL: expected 2 WGD rows, got $nrows"

  # Row 1 == clade A (declared first, LORe); row 2 == clade B (AORe).
  qA=$(awk '!/^#/{print $3; exit}' "$sum")
  rA=$(awk '!/^#/{print $4; exit}' "$sum")
  qB=$(awk '!/^#/{c++; if(c==2){print $3; exit}}' "$sum")
  rB=$(awk '!/^#/{c++; if(c==2){print $4; exit}}' "$sum")
  echo "  qA=$qA rA=$rA   qB=$qB rB=$rB"

  awk -v qa="$qA" -v ra="$rA" -v qb="$qB" -v rb="$rB" 'BEGIN{
    if (qa < 0.3) { print "FAIL: clade A q-hat too small (no WGD detected)"; exit 1 }
    if (qb < 0.3) { print "FAIL: clade B q-hat too small (no WGD detected)"; exit 1 }
    if (ra > 0.85) { print "FAIL: clade A r-hat not < 1 (LORe not recovered)"; exit 1 }
    if (rb < 0.95) { print "FAIL: clade B r-hat not ~ 1 (AORe not recovered)"; exit 1 }
    if (rb - ra < 0.2) { print "FAIL: r-hat A and B not DISTINCT"; exit 1 }
  }' || exit 1

  # Invariant 8: the log prints per-WGD q AND r (two of each, not one shared r).
  line=$(loreLine "$out.log")
  [ -n "$line" ] || fail "$MODEL: no 'After WGD+LORe opt' log line"
  nq=$(grep -oE 'q\[node [0-9]+\]=' <<<"$line" | wc -l | tr -d ' ')
  nr=$(grep -oE 'r\[node [0-9]+\]=' <<<"$line" | wc -l | tr -d ' ')
  [ "$nq" -eq 2 ] || fail "$MODEL: log has $nq per-WGD q[node ..]= (expected 2)"
  [ "$nr" -eq 2 ] || fail "$MODEL: log has $nr per-WGD r[node ..]= (expected 2)"
  echo "  PASS $MODEL: two independent q-hat, two DISTINCT r-hat, per-WGD q & r logged"
done

# --- invariant 6: terminal-WGD r is pinned ------------------------------------
echo "=== terminal-WGD r pinning (--wgd A1 terminal + --wgd B1,B2 internal) ==="
out="$TMP/out_term"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 1 \
  --wgd A1 --wgd B1,B2 --lore -p "$out" > "$out.log" 2>&1 \
  || fail "terminal-WGD --lore run exited non-zero"
grep -qE 'pinned:terminal/non-target\)' "$out.log" \
  || fail "info line did not report a pinned WGD"
grep -qE 'r\[node [0-9]+\]=1\(pinned;terminal\)' "$out.log" \
  || fail "log did not show r[node X]=1(pinned;terminal) for the terminal WGD"
# The free-parameter count must exclude the terminal WGD (1 free, 1 pinned).
grep -qE 'per-event LORe r \(1 free, 1 pinned:terminal/non-target\)' "$out.log" \
  || fail "info line free/pinned r count is wrong for the terminal WGD"
echo "  PASS: terminal-WGD r pinned to 1 and excluded from the free parameters"

# --- invariant 7a: two nested *free-r* WGDs + --lore aborts -------------------
# A1,A2 (clade-A stem) is nested inside A1,B1 (the A+B ancestor); bare --lore
# makes BOTH free-r, so their subtrees overlap and r is unidentifiable.
echo "=== nested internal WGDs + --lore (expect abort with DISJOINT message) ==="
set +e
"$BIN" "${common[@]}" --rec-model UndatedDL -g 1 \
  --wgd A1,A2 --wgd A1,B1 --lore -p "$TMP/out_nest" > "$TMP/nest.log" 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || fail "nested + --lore did NOT abort (exit 0)"
grep -qE 'DISJOINT subtrees' "$TMP/nest.log" \
  || fail "nested + --lore aborted without the DISJOINT message"
echo "  PASS: nested free-r WGDs + --lore abort with the DISJOINT subtrees error"

# --- invariant 7b: nested WGDs WITHOUT --lore run fine ------------------------
echo "=== nested WGDs without --lore (expect a clean run, two q-hat) ==="
out="$TMP/out_nest2"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 1 \
  --wgd A1,A2 --wgd A1 -p "$out" > "$out.log" 2>&1 \
  || fail "nested WGDs without --lore exited non-zero"
sum="$out/reconciliations/wgdSummary.txt"
nrows=$(grep -cvE '^#' "$sum")
[ "$nrows" -eq 2 ] || fail "nested-no-lore: expected 2 WGD rows, got $nrows"
echo "  PASS: nested WGDs without --lore run and report two independent q-hat"

# --- invariant 9a: --lore-wgd pins a non-target INTERNAL WGD ------------------
# Two disjoint internal WGDs, but only clade A is the LORe target; clade B is an
# internal WGD yet must be pinned to AORe r=1 (not fitted).
echo "=== --lore-wgd selects clade A only (clade B internal but pinned) ==="
out="$TMP/out_sel"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 1 \
  --wgd A1,A2 --wgd B1,B2 --lore-wgd A1,A2 -p "$out" > "$out.log" 2>&1 \
  || fail "--lore-wgd run exited non-zero"
sum="$out/reconciliations/wgdSummary.txt"
nrows=$(grep -cvE '^#' "$sum")
[ "$nrows" -eq 2 ] || fail "--lore-wgd: expected 2 WGD rows, got $nrows"
rA=$(awk '!/^#/{print $4; exit}' "$sum")             # clade A (target)
rB=$(awk '!/^#/{c++; if(c==2){print $4; exit}}' "$sum") # clade B (pinned)
echo "  rA(target)=$rA  rB(pinned)=$rB"
awk -v ra="$rA" -v rb="$rB" 'BEGIN{
  if (rb != 1) { print "FAIL: non-target clade B r not pinned to exactly 1"; exit 1 }
  if (ra > 0.85) { print "FAIL: target clade A r not fitted < 1"; exit 1 }
}' || exit 1
grep -qE 'per-event LORe r \(1 free, 1 pinned:terminal/non-target\)' "$out.log" \
  || fail "--lore-wgd: info line should report 1 free, 1 pinned"
grep -qE 'r\[node [0-9]+\]=1\(pinned;non-target\)' "$out.log" \
  || fail "--lore-wgd: non-target WGD not logged as (pinned;non-target)"
echo "  PASS: --lore-wgd fits only the target; non-target internal WGD pinned r=1"

# --- invariant 9b: --lore-wgd target NESTED inside an (AORe) WGD --------------
# A1,A2 (target, free-r) is nested inside A1,B1 (pinned AORe). This is the
# headline relaxation: it must NOT abort, and the outer WGD must pin to r=1.
echo "=== --lore-wgd target nested inside an AORe WGD (expect clean run) ==="
out="$TMP/out_nest_sel"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 1 \
  --wgd A1,B1 --wgd A1,A2 --lore-wgd A1,A2 -p "$out" > "$out.log" 2>&1 \
  || fail "--lore-wgd nested target exited non-zero (should be allowed)"
sum="$out/reconciliations/wgdSummary.txt"
nrows=$(grep -cvE '^#' "$sum")
[ "$nrows" -eq 2 ] || fail "nested --lore-wgd: expected 2 WGD rows, got $nrows"
# Row 1 == A1,B1 (outer, pinned); row 2 == A1,A2 (inner, target).
rOuter=$(awk '!/^#/{print $4; exit}' "$sum")
awk -v r="$rOuter" 'BEGIN{
  if (r != 1) { print "FAIL: outer (non-target) WGD r not pinned to 1"; exit 1 }
}' || exit 1
grep -qE 'per-event LORe r \(1 free, 1 pinned:terminal/non-target\)' "$out.log" \
  || fail "nested --lore-wgd: info line should report 1 free, 1 pinned"
echo "  PASS: --lore-wgd target nested in an AORe WGD runs and pins the outer WGD"

# --- invariant 10: --summary-only suppresses the bulky per-family output ------
echo "=== --summary-only keeps summaries, drops per-sample bulk + CCPs ==="
out="$TMP/out_sum"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 3 \
  --wgd A1,A2 --wgd B1,B2 --lore --summary-only --cleanup-ccp \
  -p "$out" > "$out.log" 2>&1 \
  || fail "--summary-only run exited non-zero"
# kept: the global summary, per-family means, and totals. (totalTransfers.txt
# is legitimately empty under UndatedDL -- no transfer events -- so we assert on
# the per-species event totals, which are always populated.)
[ -s "$out/reconciliations/wgdSummary.txt" ] \
  || fail "--summary-only dropped wgdSummary.txt"
[ -s "$out/reconciliations/totalSpeciesEventCounts.txt" ] \
  || fail "--summary-only dropped totalSpeciesEventCounts.txt"
ls "$out"/reconciliations/summaries/*_meanSpeciesEventCounts.txt >/dev/null 2>&1 \
  || fail "--summary-only dropped the per-family mean summaries"
# dropped: every per-sample artefact under reconciliations/all/ and the CCPs.
nbulk=$(find "$out/reconciliations/all" -type f 2>/dev/null | wc -l | tr -d ' ')
[ "$nbulk" -eq 0 ] \
  || fail "--summary-only left $nbulk per-sample files in reconciliations/all/"
nccp=$(find "$out/ccps" -name '*.ccp' 2>/dev/null | wc -l | tr -d ' ')
[ "$nccp" -eq 0 ] || fail "--summary-only left $nccp .ccp files"
# contrast: a normal run DOES write per-sample files.
out2="$TMP/out_full"
"$BIN" "${common[@]}" --rec-model UndatedDL -g 3 \
  --wgd A1,A2 --wgd B1,B2 --lore -p "$out2" > "$out2.log" 2>&1 \
  || fail "full (non-summary) run exited non-zero"
nbulk2=$(find "$out2/reconciliations/all" -type f 2>/dev/null | wc -l | tr -d ' ')
[ "$nbulk2" -gt 0 ] \
  || fail "sanity: full run wrote no per-sample files (test is not meaningful)"
echo "  PASS: --summary-only keeps summaries ($nbulk all/ files vs $nbulk2 full); CCPs cleaned"

echo "PASS: all per-event-r multi-WGD invariants hold (UndatedDL + UndatedDTL)."
