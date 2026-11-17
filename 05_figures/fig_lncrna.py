#!/usr/bin/env python
"""fig_lncrna.py — Paper1 Fig 4：lncRNA 层与顺式关系

面板：
  (a) 长度分布（中位标注）
  (b) 类型构成（antisense / lincRNA / intronic）
  (c) 顺式关系证据类型
  (d) 到顺式靶基因的距离分布

用法：python fig_lncrna.py --paper-dir <PAPER_DIR>
输出：figs/fig4_lncrna.{png,pdf}
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pub_style import MM, DOUBLE, PALETTE  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-dir", required=True)
    a = ap.parse_args()
    outdir = os.path.join(a.paper_dir, "figs")
    tdir = os.path.join(a.paper_dir, "tables")
    os.makedirs(outdir, exist_ok=True)

    lf = pd.read_csv(f"{tdir}/lncrna_features.tsv", sep="\t")
    tg = pd.read_csv(f"{tdir}/lncrna_targets_curated.tsv", sep="\t")
    summ = json.load(open(f"{tdir}/lncrna_pas_summary.json"))

    fig, axes = plt.subplots(2, 2, figsize=MM(DOUBLE, 108),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.32})

    # ---- a: length ----
    ax = axes[0, 0]
    v = lf["length" if "length" in lf.columns else "exonic_len"].values
    bins = np.logspace(np.log10(200), np.log10(10000), 28)
    ax.hist(v, bins=bins, color=PALETTE["purple"], lw=0)
    ax.axvline(np.median(v), color="black", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("lncRNA transcript length (nt)")
    ax.set_ylabel("loci")
    ax.text(0.97, 0.86, f"range {v.min()}–{v.max()}", transform=ax.transAxes,
            fontsize=6, ha="right")
    ax.text(-0.16, 1.10, "a", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- b: class composition ----
    ax = axes[0, 1]
    cls = lf.lnc_class.value_counts()
    order = ["antisense", "sense-overlap", "lincRNA"]
    vals = [cls.get(k, 0) for k in order]
    cols = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"]]
    yy = np.arange(len(order))[::-1]
    ax.barh(yy, vals, color=cols, height=0.58, lw=0)
    for y0, n in zip(yy, vals):
        ax.text(n + 2, y0, f"{n}", va="center", fontsize=6.4)
    ax.set_yticks(yy)
    ax.set_yticklabels(order, fontsize=6.6)
    ax.set_xlim(0, max(vals) * 1.16)
    ax.set_xlabel("loci", fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.text(-0.24, 1.10, "b", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- c: cis-evidence types ----
    ax = axes[1, 0]
    ev = tg.evidence_type.value_counts()
    main_ev = [("antisense", ev.get("antisense", 0)),
               ("sense-overlap", ev.get("sense-overlap", 0)),
               ("divergent", ev.get("divergent", 0))]
    other = int(ev.sum() - sum(v for _, v in main_ev))
    labels = [m[0] for m in main_ev] + ["neighbour / adjacent"]
    vals = [m[1] for m in main_ev] + [other]
    yy = np.arange(len(labels))[::-1]
    ax.barh(yy, vals, color=[PALETTE["blue"], PALETTE["sky"], PALETTE["green"], PALETTE["grey"]],
            height=0.58, lw=0)
    for y0, n in zip(yy, vals):
        ax.text(n + 1.5, y0, str(n), va="center", fontsize=6.4)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels, fontsize=6.6)
    ax.set_xlim(0, max(vals) * 1.16)
    ax.set_xlabel("pairs", fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.text(-0.28, 1.10, "c", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- d: distance to target ----
    ax = axes[1, 1]
    d = tg.distance_bp.values
    bins = np.logspace(0, np.log10(max(d.max(), 1000)), 25)
    ax.hist(d, bins=bins, color=PALETTE["orange"], lw=0)
    ax.axvline(np.median(d), color="black", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("distance to cis target (bp)")
    ax.set_ylabel("pairs")
    ax.text(-0.16, 1.10, "d", transform=ax.transAxes, fontsize=10, fontweight="bold")

    fig.savefig(os.path.join(outdir, "fig4_lncrna.png"), dpi=400)
    fig.savefig(os.path.join(outdir, "fig4_lncrna.pdf"))
    plt.close(fig)
    print("wrote fig4_lncrna")


if __name__ == "__main__":
    main()
