// Regression test for the WGD extension of UndatedDLMultiModel.
//
// Item #1 (WGD_IMPLEMENTATION.md): "_hasWGD all-false must reproduce the
//   original DL likelihood exactly." With no WGD declared, the post-WGD
//   transform arrays (_uEtop, _dlclvsTop) are exact pass-throughs of the raw
//   arrays (_uE, _dlclvs), so the log-likelihood equals the original pre-WGD
//   value on identical inputs. GOLDEN_LL was captured by building this test
//   against the original header.
//
// Item #2 (WGD_PATCH4.md STEP 0): under the WHALE DLWGD parameterization
//   (q=0 recovers no-WGD), declaring a WGD and sweeping q -> 0 must make the
//   log-likelihood converge monotonically back to GOLDEN_LL.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <IO/GeneSpeciesMapping.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>
#include <util/enums.hpp>

#include "../src/ale/UndatedDLMultiModel.hpp"

// Captured from the original (pre-WGD) code; see header comment.
// (simulated_2 species tree, family 0_pruned, default 0.2 dup/loss rates.)
static const double GOLDEN_LL = -58.689981412938;
static const double TOL = 1e-9;

// Build a RecModelInfo matching the defaults the test relies on:
// UndatedDL, single gamma category, uniform origination, no species-tree
// pruning (keeps indexing simple), DL enabled.
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
  if (argc < 4) {
    std::fprintf(stderr,
                 "usage: %s <speciesTree.newick> <geneTrees.newick> "
                 "<mapping.link>\n",
                 argv[0]);
    return 2;
  }
  const std::string speciesTreeFile = argv[1];
  const std::string geneTreeFile = argv[2];
  const std::string mappingFile = argv[3];

  PLLRootedTree speciesTree(speciesTreeFile, /*isFile=*/true);

  GeneSpeciesMapping mapping;
  mapping.fill(mappingFile, geneTreeFile);

  // Build the conditional clade probabilities from the gene tree distribution
  // and serialize them so MultiModel can unserialize.
  ConditionalClades ccp(geneTreeFile, /*likelihoods=*/"", CCPRooting::UNIFORM);
  const std::string ccpFile = std::string(geneTreeFile) + ".regtest.ccp";
  ccp.serialize(ccpFile);

  auto info = makeInfo();

  UndatedDLMultiModel<ScaledValue> model(speciesTree, mapping, info, ccpFile);
  // Use the model's default rates (0.2 dup / 0.2 loss); do not optimize, so the
  // result is fully deterministic and comparable to the captured golden value.
  double ll = model.computeLogLikelihood();

  std::remove(ccpFile.c_str());

  std::printf("WGD regression: all-false log-likelihood = %.12f\n", ll);
  std::printf("                golden (original)        = %.12f\n", GOLDEN_LL);

  if (std::fabs(ll - GOLDEN_LL) > TOL) {
    std::fprintf(stderr,
                 "FAIL: _hasWGD all-false does NOT reproduce the original "
                 "likelihood (|diff| = %.3e > %.3e)\n",
                 std::fabs(ll - GOLDEN_LL), TOL);
    return 1;
  }
  std::printf("PASS: _hasWGD all-false reproduces the original likelihood.\n");

  // Validation item #2: q->0 convergence to the no-WGD likelihood.
  // Under the WHALE DLWGD form (q=0 recovers no-WGD), declaring a WGD and
  // sweeping q -> 0 must make the log-likelihood converge monotonically to the
  // all-false value GOLDEN_LL. (This check is only meaningful because STEP 0
  // switched the closed forms from the Binomial(2,q) model to the WHALE form.)
  corax_rnode_t *wgdNode = nullptr;
  for (auto node : speciesTree.getInnerNodes()) {
    if (node->parent != nullptr && node->left != nullptr) {
      wgdNode = node; // internal, non-root branch
      break;
    }
  }
  if (!wgdNode) {
    std::fprintf(stderr, "FAIL: no internal non-root branch for WGD sweep\n");
    return 1;
  }

  const double qs[] = {0.5, 0.1, 0.01, 1e-4};
  const unsigned int nq = sizeof(qs) / sizeof(qs[0]);
  double prevDist = -1.0; // |ll(q) - GOLDEN_LL| from the previous (larger) q
  double lastLL = 0.0;
  std::printf("WGD q-sweep on branch %u:\n", wgdNode->node_index);
  for (unsigned int i = 0; i < nq; ++i) {
    model.setWGD(wgdNode->node_index, qs[i]);
    double llq = model.computeLogLikelihood();
    double dist = std::fabs(llq - GOLDEN_LL);
    std::printf("  q=%.4g -> ll=%.12f  |ll-noWGD|=%.3e\n", qs[i], llq, dist);
    if (!std::isfinite(llq)) {
      std::fprintf(stderr, "FAIL: WGD likelihood is not finite at q=%.4g\n",
                   qs[i]);
      return 1;
    }
    // Monotone convergence: each smaller q must move strictly closer to the
    // no-WGD likelihood than the previous (larger) q did.
    if (prevDist >= 0.0 && dist >= prevDist) {
      std::fprintf(stderr,
                   "FAIL: not converging to no-WGD as q->0 (|diff| %.3e at "
                   "q=%.4g >= %.3e at the previous larger q)\n",
                   dist, qs[i], prevDist);
      return 1;
    }
    prevDist = dist;
    lastLL = llq;
  }
  // At q=1e-4 the WGD is a tiny perturbation; the likelihood must be within
  // CONV_TOL of the no-WGD value.
  const double CONV_TOL = 1e-2;
  if (std::fabs(lastLL - GOLDEN_LL) > CONV_TOL) {
    std::fprintf(stderr,
                 "FAIL: q=1e-4 did not converge to no-WGD (|diff| = %.3e > "
                 "%.3e)\n",
                 std::fabs(lastLL - GOLDEN_LL), CONV_TOL);
    return 1;
  }
  std::printf("PASS: likelihood converges monotonically to the no-WGD value as "
              "q->0.\n");

  // ---- LORe (delayed rediploidization) checks (WGD_LORE.md) ----
  // STEP 5 HARD GATE: at r = 1 (default), the LORe build must reproduce the
  // Patch-4 WHALE log-likelihood. References captured from the pre-LORe binary
  // for a WGD on this branch (simulated_2, family 0_pruned).
  struct QRef {
    double q, ll;
  };
  const QRef whaleRef[] = {
      {0.5, -59.140189160263},
      {0.1, -58.767016303227},
      {0.01, -58.697442712799},
      {1e-4, -58.690055767353},
  };
  const double GATE_TOL = 1e-9;
  model.setResolutionProb(1.0);
  std::printf("LORe r=1 gate (must reproduce WHALE):\n");
  for (const auto &ref : whaleRef) {
    model.setWGD(wgdNode->node_index, ref.q);
    double llq = model.computeLogLikelihood();
    double d = std::fabs(llq - ref.ll);
    std::printf("  q=%.4g r=1 -> ll=%.12f  (WHALE %.12f, |diff|=%.2e)\n", ref.q,
                llq, ref.ll, d);
    if (d > GATE_TOL) {
      std::fprintf(stderr,
                   "FAIL: LORe r=1 does NOT reproduce WHALE at q=%.4g "
                   "(|diff|=%.3e > %.3e)\n",
                   ref.q, d, GATE_TOL);
      return 1;
    }
  }
  std::printf("PASS: LORe r=1 reproduces the WHALE log-likelihood (hard gate).\n");

  // Sanity: with a WGD (q=0.5) fixed, varying the resolution probability r must
  // change the likelihood (the LORe machinery is live) and stay finite.
  model.setWGD(wgdNode->node_index, 0.5);
  const double rs[] = {1.0, 0.5, 0.1, 1e-3};
  double prevLLr = 0.0;
  bool movedWithR = false;
  std::printf("LORe r-sweep at q=0.5 on branch %u:\n", wgdNode->node_index);
  for (unsigned int i = 0; i < sizeof(rs) / sizeof(rs[0]); ++i) {
    model.setResolutionProb(rs[i]);
    double llr = model.computeLogLikelihood();
    std::printf("  r=%.4g -> ll=%.12f\n", rs[i], llr);
    if (!std::isfinite(llr)) {
      std::fprintf(stderr, "FAIL: LORe likelihood not finite at r=%.4g\n", rs[i]);
      return 1;
    }
    if (i > 0 && std::fabs(llr - prevLLr) > 1e-9) {
      movedWithR = true;
    }
    prevLLr = llr;
  }
  if (!movedWithR) {
    std::fprintf(stderr,
                 "FAIL: likelihood does not respond to r (LORe machinery dead)\n");
    return 1;
  }
  std::printf("PASS: LORe machinery is live (likelihood responds to r).\n");

  model.setResolutionProb(1.0); // restore default
  return 0;
}
