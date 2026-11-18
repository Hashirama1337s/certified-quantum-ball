"""Draw the README figures of certified-quantum-ball from the data the repository ships.

Usage:  python figures/make_figures.py <repo root> [--png]
Needs Python 3.11 and matplotlib 3.11 (nothing else). Writes figures/*.svg next to this script.
(scorecard.svg at the repository root is drawn by tools/make_scorecard.py.)

Every plotted number is read from shipped JSON files (no pickle is loaded):
  results/result_qb1.json    QB1 (54caa4012c): 691 certified point values of b_nn (sign_map.rows; scale 'hat' =
                             b_nn e^{-4x} in the centre region, 'b' = b_nn elsewhere), the sonic values, c0
  results/result_qb1b.json   QB1b (4d3b12ef39): the 2261 certified cells (PQ2b.cells) and the zero (PQ2b.zero)
  results/result_qb2.json    QB2 (f2bcfadb3e): Psi2 enclosure of the first sealed run (PQ2-8.N8)
  results/result_qb2b.json   QB2b (373f306e9f): Psi2 route A (N8) and route B (C6), x_V (N7)
  results/result_qb2c.json   QB2c (7aab7efe5a): the blind replication rows (PQ2c-4.rows), converted
Plotting uses floats; the certified enclosures are far thinner than any line drawn here.
"""
import json
import math
import os
import sys
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle

BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
NEG, POS, GOLD, VIOLET = "#3987e5", "#d95926", "#c98500", "#9085e9"
W_IN = 12.5

plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 15, "axes.labelsize": 16, "xtick.labelsize": 14, "ytick.labelsize": 14,
    "axes.facecolor": BG, "figure.facecolor": "none", "savefig.facecolor": "none",
    "axes.edgecolor": AXC, "axes.labelcolor": TXT, "xtick.color": AXC, "ytick.color": AXC,
    "text.color": TXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 13, "lines.linewidth": 2.0, "axes.axisbelow": True,
})
PNG = "--png" in sys.argv


def new_fig(h_in):
    fig = plt.figure(figsize=(W_IN, h_in))
    fig.patches.append(FancyBboxPatch((0.004, 0.006), 0.992, 0.988, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, facecolor=BG, edgecolor=BORDER, linewidth=1.5,
                                      zorder=-10))
    return fig


def title(fig, main, sub):
    h = fig.get_figheight()
    fig.text(0.035, 1 - 0.40 / h, main, fontsize=20, fontweight="bold", color=TXT, va="top")
    fig.text(0.035, 1 - 0.84 / h, sub, fontsize=14, color=AXC, va="top")


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": None})
    with open(path, encoding="utf-8") as fh:
        svg = fh.read().replace("font-family: 'Segoe UI';", "font-family: 'Segoe UI', 'Helvetica', 'Arial', sans-serif;")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(svg)
    if PNG:
        prev = os.path.join(out, os.pardir, "_preview")
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")), format="png", dpi=110)
    plt.close(fig)


def mid(b):
    return float(b[0]) if isinstance(b, list) else float(b)


def load(root, rel):
    return json.load(open(os.path.join(root, rel), encoding="utf-8"))


def fig_bach(root, out):
    q1 = load(root, "results/result_qb1.json")
    q1b = load(root, "results/result_qb1b.json")
    rows = q1["sign_map"]["rows"]
    pts = sorted((mid(r["x"]), mid(r["value"]), r["scale"]) for r in rows)
    top = [(x, v) for x, v, s in pts if s == "b" and x >= -3]
    hat = [(x, v) for x, v, s in pts if s == "hat"] + [(x, v * math.exp(-4 * x)) for x, v, s in pts if s == "b" and x < -3]
    hat.sort()
    c0 = mid(q1["centre"]["c0"]) if "c0" in q1["centre"] else None
    z = q1b["PQ2b"]["zero"]["zero_enclosure_decimal"]
    x0 = 0.5 * (float(z[0]) + float(z[1]))
    bson = mid(q1["sonic"]["b_nn"])
    cells = q1b["PQ2b"]["cells"]
    nneg = sum(1 for c in cells if c[2] < 0)
    npos = sum(1 for c in cells if c[2] > 0)

    fig = new_fig(10.6)
    title(fig, "The quantum fingerprint: the Bach energy b_nn = r⁴B_(nn) of the critical ball",
          "Certified point values (QB1) and certified sign cells (QB1b). Negative in the core, positive in the shell, one zero.")
    ax1 = fig.add_axes([0.085, 0.47, 0.88, 0.38])
    ax1.plot([p[0] for p in top], [p[1] for p in top], color=ACC, lw=2.2)
    ax1.axhline(0, color=AXC, lw=1)
    ax1.axvline(0, color=GOLD, lw=1.4, ls="--")
    ax1.axvline(x0, color=VIOLET, lw=1.4)
    ax1.annotate("sonic point x = 0\nb_nn = %.4f" % bson, (0, bson), (1.0, bson + 0.12), color=GOLD, fontsize=13,
                 arrowprops=dict(arrowstyle="-", color=GOLD, lw=1))
    ax1.annotate("the one zero\nx0 = %.7f" % x0, (x0, 0), (-2.9, 0.33), color=VIOLET, fontsize=13,
                 arrowprops=dict(arrowstyle="-", color=VIOLET, lw=1))
    ax1.set_xlim(-3, 12)
    ax1.set_ylabel("b_nn")
    ax1.set_xlabel("x = ln(r / (−t))   (KHA time)")
    ax2 = fig.add_axes([0.085, 0.16, 0.88, 0.21])
    ax2.plot([p[0] for p in hat if p[0] <= -3], [p[1] for p in hat if p[0] <= -3], color=NEG, lw=2.2)
    if c0 is not None:
        ax2.axhline(c0, color=AXC, ls="--", lw=1.2)
        ax2.text(-11.9, c0 + 6, "c0 = %.4f" % c0, color=AXC, fontsize=13)
    ax2.set_xlim(-12, -3)
    ax2.set_ylabel("b_nn e^{−4x}")
    ax2.set_xlabel("x (approach to the centre)")
    ax3 = fig.add_axes([0.085, 0.035, 0.88, 0.04])
    for a, b, sg, _ in cells:
        lo, hi = float(Fraction(a)), float(Fraction(b))
        ax3.add_patch(Rectangle((lo, 0), hi - lo, 1, color=NEG if sg < 0 else POS, lw=0))
    ax3.set_xlim(-12, 12)
    ax3.set_ylim(0, 1)
    ax3.set_yticks([])
    ax3.grid(False)
    ax3.tick_params(labelsize=11)
    fig.text(0.085, 0.085, "%d certified cells on [−12, 12]: %d negative (blue), %d positive (orange); one zero, "
             "enclosure width %.1e" % (len(cells), nneg, npos, float(q1b["PQ2b"]["zero"]["zero_enclosure_width"])),
             fontsize=12.5, color=AXC)
    save(fig, out, "bach-profile.svg")


def fig_core(root, out):
    q1b = load(root, "results/result_qb1b.json")
    q2b = load(root, "results/result_qb2b.json")
    z = q1b["PQ2b"]["zero"]["zero_enclosure_decimal"]
    r0 = math.exp(0.5 * (float(z[0]) + float(z[1])))
    xv = mid(q2b["PQ2b-5"]["N8"]["x_V"])
    rv = math.exp(xv)
    fig = new_fig(7.6)
    title(fig, "Where the quantum source points, at one instant (t = −1, so r = e^x)",
          "Shading: sign of the Bach energy b_nn, the same for every free massless field of spin ≤ 1.")
    ax = fig.add_axes([0.03, 0.04, 0.55, 0.80])
    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1.18, color=POS, alpha=0.28, lw=0))
    ax.add_patch(Circle((0, 0), r0, color=NEG, alpha=0.55, lw=0))
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec=GOLD, lw=2, ls="--"))
    ax.add_patch(Circle((0, 0), rv, fill=False, ec=ACC, lw=2.4))
    ax.plot([0], [0], "o", color=TXT, ms=4)
    tx = 0.62
    fig.text(tx, 0.74, "core  r < %.5f:  b_nn < 0" % r0, color=NEG, fontsize=15)
    fig.text(tx, 0.66, "shell r > %.5f:  b_nn > 0" % r0, color=POS, fontsize=15)
    fig.text(tx, 0.56, "dashed: sonic sphere r = 1", color=GOLD, fontsize=15)
    fig.text(tx, 0.46, "solid: peak compactness, r = e^{x_V} = %.4f" % rv, color=ACC, fontsize=15)
    fig.text(tx, 0.405, "Ψ₂ is read here", color=ACC, fontsize=15)
    fig.text(tx, 0.26, "The picture shrinks self-similarly as t → 0⁻.\n"
             "Conditional on a Hadamard state respecting\nself-similarity (existence not proved).", color=AXC, fontsize=13)
    save(fig, out, "core-and-shell.svg")


def fig_psi(root, out):
    q2 = load(root, "results/result_qb2.json")["PQ2-8"]["N8"]["Psi2"]
    q2b = load(root, "results/result_qb2b.json")["PQ2b-5"]
    A, B = q2b["N8"]["Psi2"], q2b["C6"]["Psi2_routeB"]
    rows = load(root, "results/result_qb2c.json")["PQ2c-4"]["rows"]
    ref = float(A[0])
    fig = new_fig(6.4)
    title(fig, "Ψ₂: the certified sign and size of the log-growing push on the peak compactness",
          "Left: the enclosure against 0. Right: zoom, offsets from route A in units of 1e-12.")
    ax = fig.add_axes([0.06, 0.16, 0.36, 0.62])
    ax.axvline(0, color=AXC, lw=1.4)
    ax.plot([float(A[0])], [0], "o", color=ACC, ms=9)
    ax.set_xlim(-0.6, 4.2)
    ax.set_yticks([])
    ax.text(0.05, 0.25, "0", color=AXC, fontsize=14)
    ax.text(float(A[0]) - 1.3, 0.25, "Ψ₂ = %.10f" % float(A[0]), color=ACC, fontsize=13)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("Ψ₂ (per unit K·ε·s·e^{2s}, KHA clock)")
    ax2 = fig.add_axes([0.52, 0.16, 0.44, 0.62])
    items = [("route A (QB2b, N8)", A, ACC), ("route B (QB2b, C6)", B, GOLD), ("QB2 first run (±1.15e-7, drawn clipped)", q2, AXC)]
    for i, (lab, b, col) in enumerate(items):
        m, r = (float(b[0]) - ref) * 1e12, float(b[1]) * 1e12
        lo, hi = max(m - r, -30), min(m + r, 30)
        ax2.plot([lo, hi], [i, i], color=col, lw=5, solid_capstyle="butt")
        ax2.text(-29.5, i + 0.22, lab, color=col, fontsize=12.5)
    blind = [(k, (float(v["converted"][0]) - ref) * 1e12) for k, v in rows.items() if k != "headline Psi_per_unit_K_eps"]
    ax2.plot([b for _, b in blind], [3] * len(blind), "D", color=VIOLET, ms=7)
    ax2.text(-29.5, 3.22, "blind replication, 6 runs, converted (double precision, no enclosure)", color=VIOLET,
             fontsize=12.5)
    ax2.set_xlim(-30, 30)
    ax2.set_ylim(-0.6, 3.8)
    ax2.set_yticks([])
    ax2.set_xlabel("(value − route A) × 1e12")
    fig.text(0.06, 0.035, "Positive: toward trapping. The sign is the same for every free massless field of spin ≤ 1, "
             "in any Hadamard state respecting self-similarity.", fontsize=12, color=AXC)
    save(fig, out, "psi2-number-line.svg")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    root = os.path.abspath(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(os.getcwd(), "figures")
    os.makedirs(out, exist_ok=True)
    fig_bach(root, out)
    fig_core(root, out)
    fig_psi(root, out)
    print("wrote bach-profile.svg, core-and-shell.svg, psi2-number-line.svg in", out)


if __name__ == "__main__":
    main()
