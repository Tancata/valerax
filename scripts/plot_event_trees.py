#!/usr/bin/env python3
"""
Draw species-tree event maps from kalerax/AleRax reconciliation output.

For each model (e.g. DTL, DTL+WGD, DTL+LORe) this renders one PDF page showing
the species tree with, mapped to every branch:
    cp = total gene copies   D = duplications   T = transfers
    L  = losses              O = originations
(posterior means, summed over families), branches coloured by gene copies, and
the genome-wide D/T/L/O totals in the page header so models can be compared.

Input per model is a kalerax/AleRax `totalSpeciesEventCounts.txt`
(comma-separated, header:
 species_label, speciations, duplications, losses, transfers, presence,
 origination, copies, singletons, transfers_to),
found in each run's `reconciliations/` directory.

Internal species-tree branches are matched to the kalerax convention
`Node_<tipA>_<tipB>_0` by taking the MRCA of tipA and tipB, so the species tree
newick only needs the tip names (no internal labels required).

Usage
-----
  python3 plot_event_trees.py \
      --tree species_tree.nw --out event_maps.pdf [--png] \
      --model "DTL (null)::runs/dtl_null/reconciliations/totalSpeciesEventCounts.txt" \
      --model "DTL + WGD::runs/dtl_wgd/reconciliations/totalSpeciesEventCounts.txt" \
      --model "DTL + WGD + LORe::runs/dtl_lore/reconciliations/totalSpeciesEventCounts.txt" \
      --wgd-branch Taxon_A,Taxon_B          # highlight the WGD branch (MRCA of the two tips)

Requires: ete3, matplotlib.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors
from ete3 import Tree

EVENT_KEYS = ("duplications", "transfers", "losses", "origination", "copies")


def load_events(path):
    """Read a totalSpeciesEventCounts.txt into {species_label: {event: value}}."""
    d = {}
    with open(path) as fh:
        hdr = [h.strip() for h in fh.readline().split(",")]
        idx = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            p = [x.strip() for x in line.split(",")]
            if len(p) < len(hdr):
                continue
            d[p[idx["species_label"]]] = {k: float(p[idx[k]]) for k in EVENT_KEYS}
    return d


def fmt(v):
    if v >= 1000:
        return f"{v / 1000:.1f}k"
    if v >= 10:
        return f"{v:.0f}"
    return f"{v:.1f}"


def layout(t):
    """Cladogram layout: x = depth in edges, y = leaf order (ladderized)."""
    t.ladderize()
    leaves = t.get_leaves()
    ypos = {lf: i for i, lf in enumerate(leaves)}
    xpos = {}
    for n in t.traverse("preorder"):
        xpos[n] = 0 if n.is_root() else xpos[n.up] + 1
    for n in t.traverse("postorder"):
        if not n.is_leaf():
            ypos[n] = sum(ypos[c] for c in n.children) / len(n.children)
    return ypos, xpos, leaves


def resolve_wgd_node(t, spec):
    """Resolve --wgd-branch (an internal label, or 'tipA,tipB' MRCA) to a node."""
    if not spec:
        return None
    if "," in spec:
        a, b = spec.split(",", 1)
        return t.get_common_ancestor([a.strip(), b.strip()])
    # treat as a kalerax Node_<tipA>_<tipB>_0 label
    if spec.startswith("Node_"):
        parts = spec.split("_")
        return t.get_common_ancestor([parts[1], parts[2]])
    return None


def build_maps(t, ev):
    """Map species tree nodes <-> event-file labels (tips by name; internal by MRCA)."""
    node_to_label, label_to_node = {}, {}
    for n in t.get_leaves():
        if n.name in ev:
            node_to_label[n] = n.name
            label_to_node[n.name] = n
    for lbl in ev:
        if lbl.startswith("Node_"):
            parts = lbl.split("_")
            try:
                anc = t.get_common_ancestor([parts[1], parts[2]])
                node_to_label[anc] = lbl
                label_to_node[lbl] = anc
            except Exception:
                pass
    return node_to_label, label_to_node


def plot_model(ax, t, ypos, xpos, leaves, ev, node_to_label, label_to_node,
               title, wgd_node):
    allcopies = [ev[lbl]["copies"] for lbl in ev if lbl in label_to_node]
    norm = colors.Normalize(vmin=min(allcopies), vmax=max(allcopies))
    cmap = cm.viridis
    maxx = max(xpos.values())
    # branches (coloured by gene copies; WGD branch in crimson)
    for n in t.traverse():
        if n.is_root():
            continue
        lbl = node_to_label.get(n)
        col = cmap(norm(ev[lbl]["copies"])) if lbl in ev else "0.6"
        is_wgd = (n is wgd_node)
        ax.plot([xpos[n.up], xpos[n]], [ypos[n], ypos[n]],
                color=("crimson" if is_wgd else col),
                lw=(4.0 if is_wgd else 2.0), solid_capstyle="round", zorder=2)
    # vertical connectors
    for n in t.traverse():
        if n.is_leaf():
            continue
        ys = [ypos[c] for c in n.children]
        ax.plot([xpos[n], xpos[n]], [min(ys), max(ys)], color="0.4", lw=1.2, zorder=1)
    # per-branch annotation boxes
    for n in t.traverse():
        lbl = node_to_label.get(n)
        if lbl not in ev:
            continue
        e = ev[lbl]
        x, y = xpos[n], ypos[n]
        txt = (f"cp {fmt(e['copies'])}\nD {fmt(e['duplications'])}  "
               f"T {fmt(e['transfers'])}\nL {fmt(e['losses'])}  O {fmt(e['origination'])}")
        if n.is_leaf():
            ax.text(x + 0.12, y, "  " + n.name, va="center", ha="left",
                    fontsize=6.5, fontweight="bold")
            ax.text(x + 0.12, y - 0.42, txt, va="top", ha="left",
                    fontsize=4.7, color="#222")
        else:
            is_wgd = (n is wgd_node)
            box = dict(boxstyle="round,pad=0.15",
                       fc=("#ffecec" if is_wgd else "#eef3f7"),
                       ec=("crimson" if is_wgd else "#9bb"), lw=0.6, alpha=0.95)
            xpar = xpos[n.up] if n.up else x - 1
            ax.text((xpar + x) / 2, y + 0.10, txt, va="bottom", ha="center",
                    fontsize=4.6, color="#111", bbox=box, zorder=4)
    # genome-wide totals (for cross-model comparison)
    tot = {k: sum(ev[l][k] for l in ev if l in label_to_node)
           for k in ("duplications", "transfers", "losses", "origination")}
    sub = (f"genome-wide totals:  D {fmt(tot['duplications'])}   "
           f"T {fmt(tot['transfers'])}   L {fmt(tot['losses'])}   "
           f"O {fmt(tot['origination'])}")
    ax.set_title(title + "\n" + sub, fontsize=12, fontweight="bold")
    ax.set_xlim(-0.4, maxx + 3.2)
    ax.set_ylim(-1.2, len(leaves) + 0.5)
    ax.axis("off")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("total gene copies (branch)", fontsize=8)
    cb.ax.tick_params(labelsize=7)


def main():
    ap = argparse.ArgumentParser(description="Species-tree event maps from kalerax output.")
    ap.add_argument("--tree", required=True, help="species tree newick")
    ap.add_argument("--out", required=True, help="output PDF (one page per --model)")
    ap.add_argument("--model", action="append", required=True, metavar="LABEL::FILE",
                    help="'<title>::<totalSpeciesEventCounts.txt>' (repeatable)")
    ap.add_argument("--wgd-branch", default=None,
                    help="WGD branch to highlight: a 'Node_..' label or 'tipA,tipB' (MRCA)")
    ap.add_argument("--png", action="store_true", help="also write one PNG per page")
    ap.add_argument("--footer",
                    default="cp=gene copies  D=duplications  T=transfers  "
                            "L=losses  O=originations  (posterior means over families)",
                    help="footer caption")
    args = ap.parse_args()

    t = Tree(args.tree, format=1)
    ypos, xpos, leaves = layout(t)
    wgd_node = resolve_wgd_node(t, args.wgd_branch)

    with PdfPages(args.out) as pdf:
        for i, spec in enumerate(args.model):
            if "::" not in spec:
                raise SystemExit(f"--model must be 'LABEL::FILE', got: {spec}")
            title, path = spec.split("::", 1)
            ev = load_events(path)
            n2l, l2n = build_maps(t, ev)
            fig, ax = plt.subplots(figsize=(13, 16))
            plot_model(ax, t, ypos, xpos, leaves, ev, n2l, l2n, title, wgd_node)
            fig.text(0.5, 0.012, args.footer, ha="center", fontsize=7.5, color="#444")
            pdf.savefig(fig, bbox_inches="tight")
            if args.png and args.out.endswith(".pdf"):
                fig.savefig(f"{args.out[:-4]}_p{i + 1}.png", dpi=120, bbox_inches="tight")
            plt.close(fig)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
