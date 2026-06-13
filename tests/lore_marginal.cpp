// LORe resolution-branch MARGINAL readout (WGD_LORE_marginal.md, STEP 5).
//
// Recovers, per lineage, the posterior over the branch on which the U->R
// resolution (ohnolog divergence) committed, under the FITTED global r, via the
// U-aware backtrace (UndatedDLMultiModel::sampleResolutionCommits). This is the
// proper marginal that replaces the r in {0,1} profile proxy.
//
// Checks:
//   5.1  r=1 => every commit lands on the WGD branch (regardless of topology).
//   5.2  weight==inside consistency (always-on guard inside sampleResolution-
//        Commits; we pass check=true everywhere).
//   5.3/5.4  under the FITTED global r, the highest-posterior commit branch
//        matches the simulated truth (species-grouped -> terminals; ohnolog-
//        grouped -> WGD branch).
//
// usage: lore_marginal <speciesTree> <geneTree> <taxonA> <taxonB> <AORe|LORe>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <IO/GeneSpeciesMapping.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/Random.hpp>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>
#include <util/enums.hpp>

#include "../src/ale/UndatedDLMultiModel.hpp"

static RecModelInfo makeInfo() {
  return RecModelInfo(RecModel::UndatedDL, RecOpt::LBFGSB,
                      /*perFamilyRates=*/false, /*gammaCategories=*/1,
                      OriginationStrategy::UNIFORM, /*pruneSpeciesTree=*/false,
                      /*rootedGeneTree=*/false, /*forceGeneTreeRoot=*/false,
                      /*madRooting=*/false, /*branchLengthThreshold=*/-1.0,
                      TransferConstaint::NONE, /*noDup=*/false, /*noDL=*/false,
                      /*noTL=*/false, /*fractionMissingFile=*/"",
                      /*memorySavings=*/false);
}

int main(int argc, char **argv) {
  if (argc < 6) {
    std::fprintf(stderr, "usage: %s <speciesTree> <geneTree> <taxonA> <taxonB> "
                         "<AORe|LORe>\n",
                 argv[0]);
    return 2;
  }
  Random::setSeed(42);
  const std::string speciesTreeFile = argv[1];
  const std::string geneTreeFile = argv[2];
  const std::string taxonA = argv[3];
  const std::string taxonB = argv[4];
  const std::string expected = argv[5];
  const double q = 0.5; // fixed retention; we fit r only

  PLLRootedTree speciesTree(speciesTreeFile, /*isFile=*/true);
  GeneSpeciesMapping mapping;
  mapping.fill("", geneTreeFile);
  ConditionalClades ccp(geneTreeFile, "", CCPRooting::UNIFORM);
  const std::string ccpFile = geneTreeFile + ".loremarg.ccp";
  ccp.serialize(ccpFile);
  auto info = makeInfo();
  UndatedDLMultiModel<ScaledValue> model(speciesTree, mapping, info, ccpFile);
  std::remove(ccpFile.c_str());

  auto labelToNode = speciesTree.getLabelToNode(false);
  auto *wNode = speciesTree.getLCA(labelToNode[taxonA], labelToNode[taxonB]);
  unsigned int w = wNode->node_index;
  unsigned int bA = labelToNode[taxonA]->node_index;
  unsigned int bB = labelToNode[taxonB]->node_index;
  model.setWGD(w, q);
  std::printf("WGD branch = node %u; lineage branches %s=%u, %s=%u\n", w,
              taxonA.c_str(), bA, taxonB.c_str(), bB);

  const unsigned int N = 40000;
  std::vector<double> commits;

  // --- STEP 5.1: r = 1 => all commits on the WGD branch ---
  model.setResolutionProb(1.0);
  model.computeLogLikelihood();
  model.sampleResolutionCommits(N, commits, /*check=*/true);
  double tot1 = 0.0, onW = 0.0;
  for (unsigned int e = 0; e < commits.size(); ++e) {
    tot1 += commits[e];
    if (e == w) onW += commits[e];
  }
  std::printf("STEP 5.1 (r=1): commits total=%.0f, on WGD branch=%.0f (%.4f)\n",
              tot1, onW, tot1 > 0 ? onW / tot1 : 0.0);
  if (tot1 == 0 || onW != tot1) {
    std::fprintf(stderr, "FAIL 5.1: not all commits on the WGD branch at r=1\n");
    return 1;
  }
  std::printf("PASS 5.1: at r=1 every commit lands on the WGD branch.\n");

  // --- fit the global r by a likelihood grid (q fixed) ---
  double bestR = 1.0, bestLL = -1e300;
  for (int i = 1; i <= 50; ++i) {
    double r = i / 50.0;
    model.setResolutionProb(r);
    double ll = model.computeLogLikelihood();
    if (ll > bestLL) {
      bestLL = ll;
      bestR = r;
    }
  }
  std::printf("fitted global r-hat = %.3f (ll=%.5f)\n", bestR, bestLL);

  // --- STEP 5.3/5.4: marginal under the fitted r ---
  model.setResolutionProb(bestR);
  model.computeLogLikelihood();
  model.sampleResolutionCommits(N, commits, /*check=*/true);
  double tot = 0.0;
  for (double v : commits) {
    tot += v;
  }
  auto frac = [&](unsigned int e) { return tot > 0 ? commits[e] / tot : 0.0; };
  std::printf("commit fractions: WGD(node %u)=%.3f  %s=%.3f  %s=%.3f  (n=%.0f "
              "commits)\n",
              w, frac(w), taxonA.c_str(), frac(bA), taxonB.c_str(), frac(bB),
              tot);
  // simple Monte Carlo SE on a fraction p: sqrt(p(1-p)/n)
  auto se = [&](unsigned int e) {
    double p = frac(e);
    return tot > 0 ? std::sqrt(p * (1 - p) / tot) : 0.0;
  };
  std::printf("  (SE: WGD=%.4f %s=%.4f %s=%.4f)\n", se(w), taxonA.c_str(),
              se(bA), taxonB.c_str(), se(bB));

  bool ok = false;
  if (expected == "AORe") {
    // resolution should sit on the WGD branch, fitted r -> ~1
    ok = (bestR > 0.9) && (frac(w) > 0.9);
    std::printf("expected AORe: r-hat~1 and commits on WGD branch\n");
  } else { // LORe
    // resolution should sit on the two lineage (terminal) branches, not WGD
    ok = (bestR < 1.0) && (frac(bA) > 0.3) && (frac(bB) > 0.3) &&
         (frac(w) < 0.1);
    std::printf("expected LORe: r-hat<1 and commits on the lineage branches\n");
  }
  if (!ok) {
    std::fprintf(stderr,
                 "FAIL: marginal does not match the simulated placement (%s)\n",
                 expected.c_str());
    return 1;
  }
  std::printf("PASS: highest-posterior commit branch matches the simulated "
              "truth (%s) under fitted r.\n",
              expected.c_str());
  return 0;
}
