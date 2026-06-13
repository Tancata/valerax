// LORe AORe-vs-LORe MODEL-SELECTION sanity check (WGD_LORE_marginal.md STEP 0).
//
// NOTE: this is NOT the resolution-branch marginal. It forces r in {0,1} on
// branches (a profile-likelihood comparison between two hypotheses) and is kept
// only as a coarse AORe-vs-LORe sanity check. The actual per-lineage commit-
// branch POSTERIOR under the fitted global r is computed by the U-aware
// backtrace (sampleResolutionCommits) and tested in lore_marginal.
//
// For a *species-grouped* gene tree ((Aa,Ab),(Ba,Bb)) the resolution belongs on
// the terminal branches (LORe); for an *ohnolog-grouped* ((Aa,Ba),(Ab,Bb)) it
// belongs at the WGD branch (AORe). This declares the WGD on the
// LCA(taxonA,taxonB) branch, fixes q, and checks the higher-inside-likelihood
// hypothesis matches the gene-tree topology.
//
// usage: lore_modelselection_sanity <speciesTree> <geneTree> <taxonA> <taxonB>
//                                   <AORe|LORe>     (the expected winner)

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
  const std::string speciesTreeFile = argv[1];
  const std::string geneTreeFile = argv[2];
  const std::string taxonA = argv[3];
  const std::string taxonB = argv[4];
  const std::string expected = argv[5];

  PLLRootedTree speciesTree(speciesTreeFile, /*isFile=*/true);
  GeneSpeciesMapping mapping;
  mapping.fill("", geneTreeFile); // auto: species = label prefix before first _
  ConditionalClades ccp(geneTreeFile, "", CCPRooting::UNIFORM);
  const std::string ccpFile = geneTreeFile + ".lorerec.ccp";
  ccp.serialize(ccpFile);

  auto info = makeInfo();
  UndatedDLMultiModel<ScaledValue> model(speciesTree, mapping, info, ccpFile);
  std::remove(ccpFile.c_str());

  // WGD branch = the branch leading to the LCA of taxonA and taxonB.
  auto labelToNode = speciesTree.getLabelToNode(false);
  if (!labelToNode.count(taxonA) || !labelToNode.count(taxonB)) {
    std::fprintf(stderr, "FAIL: taxa not found in species tree\n");
    return 1;
  }
  auto *wNode = speciesTree.getLCA(labelToNode[taxonA], labelToNode[taxonB]);
  if (!wNode || !wNode->left) {
    std::fprintf(stderr, "FAIL: WGD LCA is not an internal node\n");
    return 1;
  }
  unsigned int w = wNode->node_index;
  unsigned int cl = wNode->left->node_index;
  unsigned int cr = wNode->right->node_index;
  std::printf("WGD on branch to node %u (LCA %s,%s); daughter branches %u,%u\n",
              w, taxonA.c_str(), taxonB.c_str(), cl, cr);

  const double q = 0.5; // fixed retention: isolate the resolution-placement signal
  model.setWGD(w, q);

  // Config AORe: resolve immediately at the WGD branch (r=1 everywhere).
  model.setResolutionProb(1.0);
  double llAORe = model.computeLogLikelihood();

  // Config LORe: do not resolve at the WGD branch (r[w]=0); resolve on the two
  // descendant lineage branches (r=1 there, inherited from the global setting).
  model.setResolutionProb(1.0);
  model.setResolutionProbBranch(w, 0.0);
  double llLORe = model.computeLogLikelihood();

  std::printf("  AORe (resolve at WGD branch)        ll = %.10f\n", llAORe);
  std::printf("  LORe (resolve on lineage branches)  ll = %.10f\n", llLORe);
  const char *winner = (llLORe > llAORe) ? "LORe" : "AORe";
  std::printf("  -> inferred resolution placement: %s (expected %s)\n", winner,
              expected.c_str());

  if (!std::isfinite(llAORe) || !std::isfinite(llLORe)) {
    std::fprintf(stderr, "FAIL: non-finite likelihood\n");
    return 1;
  }
  if (expected != winner) {
    std::fprintf(stderr,
                 "FAIL: resolution placement not recovered (got %s, want %s)\n",
                 winner, expected.c_str());
    return 1;
  }
  std::printf("PASS: recovered the simulated resolution placement (%s).\n",
              winner);
  return 0;
}
