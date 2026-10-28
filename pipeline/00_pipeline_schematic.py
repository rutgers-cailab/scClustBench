#!/usr/bin/env python3
# ============================================================
# pipeline/00_pipeline_schematic.py
#
# Draws the pipeline schematic: data processing and embeddings (panel A), the two stability
# tests and the Stage A / Stage B decision with its three verdicts (panel B), and the
# applications and synthetic validation (panel C). The figure is pure matplotlib drawing
# and reads no pipeline outputs, so the numbers shown in the boxes are typed into the labels.
#
# Paper:    Figure 1
# Outputs:  output/paper/Figure_1.png, output/paper/Figure_1.pdf
#           (or into the directory given as the first command-line argument)
# Run:      Rscript run_pipeline.R fig1   or   python3 pipeline/00_pipeline_schematic.py [outdir]
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os, sys

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
FIGDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "paper")
os.makedirs(FIGDIR, exist_ok=True)

# Box fill colours by role; verdict colours match the ACCEPT / EMBEDDING-DEPENDENT / REJECT scheme
C_DATA="#DCE6F1"; C_EMB="#CFE8DC"; C_BENCH="#FCE7C8"; C_TEST="#E7D9EF"
C_ACCEPT="#1B9E77"; C_REJECT="#D95F02"; C_ED="#E6AB02"; C_APP="#F2F2F2"; EDGE="#444444"

fig, ax = plt.subplots(figsize=(9.2, 12.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 14); ax.axis("off")

# Draw a rounded box at (x, y) with size (w, h) and centred text; returns nothing.
def box(x, y, w, h, text, fc, fs=10, bold=False, tc="black", ec=EDGE, r=0.12):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.02,rounding_size={r}",
                 linewidth=1.2, edgecolor=ec, facecolor=fc, zorder=2))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=tc, zorder=3, wrap=True)

# Draw an arrow from (x1, y1) to (x2, y2); returns nothing.
def arrow(x1, y1, x2, y2, color=EDGE, lw=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=16,
                 linewidth=lw, color=color, zorder=1, shrinkA=2, shrinkB=2))

# Draw a free text label (panel titles); returns nothing.
def label(x, y, text, fs=10, color="#333333", bold=True, ha="left"):
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs, color=color,
            fontweight="bold" if bold else "normal", zorder=3)

# ---------- Section A: data -> embeddings -> benchmark ----------
label(0.2, 13.7, "A   Build embeddings and benchmark clustering", fs=12)
box(1.4, 12.5, 7.2, 0.85, "Six mouse-retina scRNA-seq datasets (P15–P91)\nQC + scDblFinder doublet removal + SCTransform + Harmony integration", C_DATA, fs=9.5)
arrow(5, 12.5, 5, 12.1)
box(1.0, 11.25, 3.7, 0.8, "PCA embedding\n(linear)", C_EMB, fs=9.5, bold=True)
box(5.3, 11.25, 3.7, 0.8, "scDHA embedding\n(nonlinear)", C_EMB, fs=9.5, bold=True)
ax.text(5, 11.05, "Harmony re-applied to each embedding", ha="center", va="top", fontsize=8, style="italic", color="#555")
arrow(2.85, 11.25, 4.3, 10.55); arrow(7.15, 11.25, 5.7, 10.55)
box(1.4, 9.75, 7.2, 0.75, "Benchmark 7 configurations across k on 5 quality metrics\nROGUE  ·  balance  ·  per-cluster Jaccard  ·  marker specificity  ·  silhouette", C_BENCH, fs=9)
arrow(5, 9.75, 5, 9.35)

# ---------- Section B: two stability tests + decision ----------
label(0.2, 9.05, "B   Test candidate partitions for two kinds of reproducibility", fs=12)
box(0.8, 7.9, 3.9, 0.95, "Within-embedding stability\nbootstrap-resample cells,\nre-cluster (fixed embedding)", C_TEST, fs=9)
box(5.3, 7.9, 3.9, 0.95, "Across-embedding stability\nretrain or perturb the embedding\n(scDHA, scVI, two PCA perturbation families)", C_TEST, fs=8.6, bold=True)
arrow(2.75, 7.9, 4.2, 7.15); arrow(7.25, 7.9, 5.8, 7.15)
# Decision box (Stage A screen, Stage B tau test) and the three verdicts
box(2.0, 6.15, 6.0, 0.95, "Stage A: validity and stability screen\nStage B: max R(k) over valid k ≥ τ = 0.8?", "#FFFFFF", fs=9.5, bold=True, ec="#222222")
arrow(3.2, 6.15, 1.9, 5.5, color=C_ACCEPT); arrow(5, 6.15, 5, 5.5, color=C_ED); arrow(6.8, 6.15, 8.1, 5.5, color=C_REJECT)
ax.text(1.85, 5.82, "YES", fontsize=10, fontweight="bold", color=C_ACCEPT)
ax.text(7.75, 5.82, "NO", fontsize=10, fontweight="bold", color=C_REJECT)
box(0.4, 4.55, 2.95, 0.9, "ACCEPT\nreproducible sub-structure", C_ACCEPT, fs=9, bold=True, tc="white", ec=C_ACCEPT)
box(3.525, 4.55, 2.95, 0.9, "EMBEDDING-DEPENDENT\nembedding families disagree", C_ED, fs=9, bold=True, ec=C_ED)
box(6.65, 4.55, 2.95, 0.9, "REJECT\nmax R(k) < τ", C_REJECT, fs=9, bold=True, tc="white", ec=C_REJECT)

# ---------- Section C: applications / validation ----------
label(0.2, 4.05, "C   Applications and ground-truth validation", fs=12)
box(0.5, 2.4, 2.85, 1.05, "Phase 1\nWhole retina: k* = 6,\nPCA + Leiden, ARI 0.96\nvs reference annotation", "#E6EEF8", fs=8.6)
box(3.6, 2.4, 2.85, 1.05, "Phase 2\nBipolar cells: ACCEPT\nscDHA R = 0.95, k* = 16;\naccepted by four families", "#DDF0E8", fs=8.4, bold=False, ec=C_ACCEPT)
box(6.7, 2.4, 2.85, 1.05, "Phase 3\nRods: REJECT\nmax R 0.40–0.74 in\nfour embedding families", "#FBE3D2", fs=8.4, bold=False, ec=C_REJECT)
box(1.4, 0.85, 7.2, 1.05, "Synthetic ground-truth validation (Splatter)\nAccepts structured data at the true k, including a 5% rare subpopulation\n(ARI vs truth 0.97–1.00)  ·  rejects every null (max R < τ)", "#EFE7F5", fs=8.8, ec="#7E57C2")
arrow(5, 2.4, 5, 1.9, color="#7E57C2")

# ---------- Save ----------
plt.tight_layout(pad=0.5)
for ext in ("png","pdf"):
    fig.savefig(os.path.join(FIGDIR, f"Figure_1.{ext}"), dpi=300 if ext=="png" else None, bbox_inches="tight")
print("Saved Figure_1.png/pdf")
