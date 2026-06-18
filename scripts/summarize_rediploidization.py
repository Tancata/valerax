#!/usr/bin/env python3
"""
Summarize per-family rediploidization (and tetrasomy) from a kalerax --lore run.

kalerax writes, per gene family, a
  reconciliations/summaries/<family>_meanResolutionCounts.txt
with columns:
  species_branch   expected_WGD_resolution_events   expected_tetrasomic_at_tip
(the second is the posterior-mean number of this family's ohnolog pairs that
rediploidized -- U->R committed -- on that branch; the third is the expected
number still tetrasomic at that extant tip). This script consolidates those
16k+ files into two tidy tables:

  <prefix>_distribution.tsv : long format, one row per (family, branch)
        family  species_branch  expected_resolutions  probability
     probability = expected_resolutions / (family total); sums to 1 per family.

  <prefix>_summary.tsv : one row per family
        family
        expected_n_rediploidizations      total U->R events (>1 when delayed:
                                           a delayed pair resolves once per
                                           descendant lineage it commits in)
        expected_on_WGD_branch            events at the WGD (immediate / AORe)
        expected_in_descendant_lineages   events delayed into descendants (LORe)
        frac_immediate                    on_WGD / total
        n_resolution_branches             # branches the family resolved on
        effective_n_branches              exp(Shannon entropy) of the per-branch share
        expected_tetrasomic_at_tips       expected pairs still tetrasomic today
        frac_still_tetrasomic             tetrasomic / (rediploidizations + tetrasomic)
        class                             immediate (>=0.8) / delayed (<=0.2) / mixed
        top3_branches(share;need_not_sum_to_1)

Usage
-----
  python3 summarize_rediploidization.py \
      --summaries-dir runs/dtl_lore/reconciliations/summaries \
      --wgd-branch Node_TaxonA_TaxonB_0 \
      --out-prefix results/dtl_family_rediploidization

The --wgd-branch is the kalerax internal label of the WGD branch (see
wgdSummary.txt); resolution there counts as "immediate".
"""
import argparse
import csv
import glob
import math
import os

SUFFIX = "_meanResolutionCounts.txt"


def read_family(path):
    """Return (resolutions {branch: val}, tetrasomy {branch: val})."""
    res, tet = {}, {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            br = p[0]
            r = float(p[1])
            t = float(p[2]) if len(p) > 2 else 0.0
            if r > 0:
                res[br] = r
            if t > 0:
                tet[br] = t
    return res, tet


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summaries-dir", required=True,
                    help="reconciliations/summaries dir with *_meanResolutionCounts.txt")
    ap.add_argument("--wgd-branch", required=True,
                    help="WGD branch label (resolution here = immediate); see wgdSummary.txt")
    ap.add_argument("--out-prefix", required=True,
                    help="output prefix; writes <prefix>_distribution.tsv and _summary.tsv")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.summaries_dir, "*" + SUFFIX)))
    if not files:
        raise SystemExit(f"no *{SUFFIX} files under {args.summaries_dir}")

    dist_path = args.out_prefix + "_distribution.tsv"
    summ_path = args.out_prefix + "_summary.tsv"
    n_imm = n_del = n_mix = 0

    with open(dist_path, "w", newline="") as dh, open(summ_path, "w", newline="") as sh:
        dw = csv.writer(dh, delimiter="\t")
        sw = csv.writer(sh, delimiter="\t")
        dw.writerow(["family", "species_branch", "expected_resolutions", "probability"])
        sw.writerow(["family", "expected_n_rediploidizations", "expected_on_WGD_branch",
                     "expected_in_descendant_lineages", "frac_immediate",
                     "n_resolution_branches", "effective_n_branches",
                     "expected_tetrasomic_at_tips", "frac_still_tetrasomic",
                     "class", "top3_branches(share;need_not_sum_to_1)"])
        for path in files:
            fam = os.path.basename(path)[:-len(SUFFIX)]
            res, tet = read_family(path)
            tot = sum(res.values())
            tetra = sum(tet.values())
            if tot <= 0 and tetra <= 0:
                continue
            if tot > 0:
                for b in sorted(res, key=lambda x: -res[x]):
                    dw.writerow([fam, b, f"{res[b]:.5g}", f"{res[b] / tot:.5g}"])
            imm = res.get(args.wgd_branch, 0.0)
            deld = tot - imm
            fimm = (imm / tot) if tot > 0 else 0.0
            probs = [v / tot for v in res.values()] if tot > 0 else []
            H = -sum(p * math.log(p) for p in probs if p > 0)
            eff = math.exp(H) if probs else 0.0
            denom = tot + tetra
            frac_tet = (tetra / denom) if denom > 0 else 0.0
            cls = "immediate" if fimm >= 0.8 else "delayed" if fimm <= 0.2 else "mixed"
            n_imm += cls == "immediate"
            n_del += cls == "delayed"
            n_mix += cls == "mixed"
            top3 = "; ".join(f"{b}({res[b] / tot:.2f})"
                             for b in sorted(res, key=lambda x: -res[x])[:3]) if tot > 0 else ""
            sw.writerow([fam, f"{tot:.4g}", f"{imm:.4g}", f"{deld:.4g}",
                         f"{fimm:.3g}", len(res), f"{eff:.3g}",
                         f"{tetra:.4g}", f"{frac_tet:.3g}", cls, top3])

    print(f"wrote {dist_path}")
    print(f"wrote {summ_path}  ({n_imm} immediate, {n_del} delayed, {n_mix} mixed)")


if __name__ == "__main__":
    main()
