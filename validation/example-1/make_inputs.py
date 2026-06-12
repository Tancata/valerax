#!/usr/bin/env python3
"""
Generate the AleRax inputs for the Whale.jl example-1 WGD cross-validation
(WGD_PATCH4.md STEP 4, Tier A).

For each *.ale file in this directory:
  - parse the gene-tree leaf labels from the ALE constructor string (line 2),
  - map each gene to its species (the prefix before the first '_'),
  - write a per-family mapping file mappings/<family>.link in the "treerecs"
    format (one "<gene> <species>" per line).
Then write families.txt referencing each .ale file (as gene_tree) and its
mapping. AleRax detects the leading '#' and parses the .ale as an ALE-format
CCP, which requires UNIFORM (unrooted) CCP rooting.

Paths in families.txt are written relative to --root (default: a copy of this
directory), so the file can be staged anywhere -- importantly to a path without
spaces, since AleRax's input reader does not handle spaces in paths.
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.realpath(__file__))


def parse_leaf_labels(ale_file):
    """Extract leaf labels from the newick constructor string of an .ale file."""
    with open(ale_file) as f:
        first = f.readline()
        if not first.startswith("#"):
            raise ValueError("%s does not look like an ALE file" % ale_file)
        constructor = f.readline().strip()
    # leaf labels are the tokens that immediately precede a ':' (branch length)
    # and follow a '(' or ',' (i.e. not internal-node labels, which follow ')')
    labels = re.findall(r"[(,]([^(),:;]+):", constructor)
    return labels


def species_of(gene):
    """Species code = prefix before the first underscore."""
    return gene.split("_", 1)[0]


def family_name(ale_basename):
    # OG0004512.fasta.nex.treesample.ale -> OG0004512
    return ale_basename.split(".", 1)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=HERE,
                    help="path prefix to write into families.txt (default: this dir)")
    ap.add_argument("--out", default=os.path.join(HERE, "families.txt"),
                    help="output families file")
    args = ap.parse_args()

    ale_files = sorted(f for f in os.listdir(HERE) if f.endswith(".ale"))
    if not ale_files:
        sys.exit("No .ale files found in %s" % HERE)

    mappings_dir = os.path.join(HERE, "mappings")
    os.makedirs(mappings_dir, exist_ok=True)

    all_species = set()
    families = []
    for ale in ale_files:
        fam = family_name(ale)
        labels = parse_leaf_labels(os.path.join(HERE, ale))
        if not labels:
            sys.exit("No leaf labels parsed from %s" % ale)
        link_path = os.path.join(mappings_dir, fam + ".link")
        with open(link_path, "w") as out:
            for gene in labels:
                sp = species_of(gene)
                all_species.add(sp)
                out.write("%s %s\n" % (gene, sp))
        families.append((fam, ale, "mappings/" + fam + ".link", len(labels)))

    root = os.path.abspath(args.root)
    with open(args.out, "w") as out:
        out.write("[FAMILIES]\n")
        for fam, ale, link, _ in families:
            out.write("- %s\n" % fam)
            out.write("gene_tree = %s\n" % os.path.join(root, ale))
            out.write("mapping = %s\n" % os.path.join(root, link))

    print("Wrote %d family mappings to %s" % (len(families), mappings_dir))
    print("Wrote %s (root=%s)" % (args.out, root))
    print("Species observed (%d): %s" %
          (len(all_species), ", ".join(sorted(all_species))))
    for fam, _, _, n in families:
        print("  %-12s %3d gene leaves" % (fam, n))


if __name__ == "__main__":
    main()
