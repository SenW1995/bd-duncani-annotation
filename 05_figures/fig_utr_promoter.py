#!/usr/bin/env python
"""fig_utr_promoter.py — Paper1 Fig 3：UTR、启动子与 PAS

面板：
  (a) 5'/3' UTR 长度分布直方图（对数分箱，标中位数与 B. bovis 参考）
  (b) TSS ±10 起始子 sequence logo（bits）
  (c) ATGGGG / CCCCAT 位置分布
  (d) ACACA / TGTGT / TATAT 周期性
  (e) 多聚腺苷化信号使用
  (f) 每基因 PAS 数

用法：python fig_utr_promoter.py --paper-dir <PAPER_DIR>
输出：figs/fig3_utr_promoter.{png,pdf}
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pub_style import MM, DOUBLE, PALETTE  # noqa: E402

BASE_COL = {"A": "#2E7D32", "C": "#1565C0", "G": "#EF6C00", "T": "#C62828"}


def draw_logo(ax, freq_df, pos_lo=-10, pos_hi=10):
    """标准 sequence logo：字母高度 = freq × information(bits)"""
    letters = ["A", "C", "G", "T"]
    sub = freq_df[freq_df.pos.between(pos_lo, pos_hi)]
    for _, r in sub.iterrows():
        total = sum(r[L] for L in letters)
        if total <= 0:
            continue
        p = {L: r[L] / total for L in letters}
        H = -sum(v * math.log2(v) for v in p.values() if v > 0)
        ic = 2.0 - H
        y = 0.0
        for L in sorted(letters, key=lambda x: p[x]):
            h = p[L] * ic
            if h < 0.03:
                y += h
                continue
            ax.add_patch(Rectangle((r.pos - 0.42, y), 0.84, h,
                                   facecolor=BASE_COL[L], edgecolor="none", alpha=0.35))
            ax.text(r.pos, y + h / 2, L, ha="center", va="center",
                    fontsize=5.0 + 9.0 * h, color=BASE_COL[L], fontweight="bold")
            y += h
    ax.set_xlim(pos_lo - 0.8, pos_hi + 0.8)
    ax.set_ylim(0, 2.05)
    ax.set_yticks([0, 1, 2])
    ax.axvline(0.5, ls=":", lw=0.7, c="#666666", zorder=0)
    ax.set_xticks(list(range(pos_lo, pos_hi + 1, 5)))
    ax.set_xlabel("position relative to TSS (nt)")
    ax.set_ylabel("bits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper-dir", required=True)
    a = ap.parse_args()
    outdir = os.path.join(a.paper_dir, "figs")
    tdir = os.path.join(a.paper_dir, "tables")
    os.makedirs(outdir, exist_ok=True)

    bf = pd.read_csv(f"{tdir}/tss_base_freq.tsv", sep="\t")
    prof = pd.read_csv(f"{tdir}/motif_position_profile.tsv", sep="\t")
    ul = pd.read_csv(f"{tdir}/utr_lengths.tsv", sep="\t")
    pas = pd.read_csv(f"{tdir}/pas_per_gene_curated.tsv", sep="\t")
    tss = pd.read_csv(os.path.join(a.paper_dir, "..", "Submission_Package", "02_data/tables/tss_per_gene.tsv"), sep="\t")
    tss = tss[tss.biotype == "protein_coding"]
    sig = pd.read_csv(f"{tdir}/pas_signal_usage.tsv", sep="\t")
    summ = json.load(open(f"{tdir}/utr_summary.json"))

    fig, axes = plt.subplots(3, 3, figsize=MM(DOUBLE, 168),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.52})

    # ---- a: UTR length histograms ----
    ax = axes[0, 0]
    pc = ul[ul.biotype == "protein_coding"]
    bins = np.logspace(1, 3.48, 30)
    ax.hist(pc.loc[pc.utr3_len > 0, "utr3_len"], bins=bins, color=PALETTE["orange"],
            alpha=0.65, lw=0, label="3′ UTR")
    ax.hist(pc.loc[pc.utr5_len > 0, "utr5_len"], bins=bins, color=PALETTE["blue"],
            alpha=0.65, lw=0, label="5′ UTR")
    ax.axvline(summ["utr5_median"], ls="--", lw=0.9, c=PALETTE["blue"])
    ax.axvline(summ["utr3_median"], ls="--", lw=0.9, c=PALETTE["orange"])
    ax.axvline(152, ls=":", lw=0.8, c="#999999")
    ax.axvline(116, ls=":", lw=0.8, c="#999999")
    ax.set_xscale("log")
    ax.set_xlim(10, 3000)
    ax.set_xlabel("UTR length (nt)")
    ax.set_ylabel("genes")
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([0], [0], ls=":", lw=0.8, c="#999999"))
    l.append("B. bovis 152 / 116")
    ax.legend(h, l, loc="upper left", fontsize=5.6, handlelength=1.3, frameon=False)
    ax.text(-0.18, 1.10, "a", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- b: TSS base composition ----
    ax = axes[0, 1]
    sub = bf[bf.pos.between(-20, 20)]
    for L in "ACGT":
        ax.plot(sub.pos, sub[L], color=BASE_COL[L], lw=1.2, label=L)
    ax.axvline(0, ls=":", lw=0.7, c="#666666")
    for x, y in [(-2, 0.55), (-1, 0.43), (0, 0.63)]:
        ax.annotate("", xy=(x, y), xytext=(x, y + 0.08),
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="#444444"))
    ax.set_xlim(-20, 21)
    ax.set_ylim(0, 0.72)
    ax.set_xticks([-20, -10, 0, 10, 20])
    ax.set_xlabel("position relative to TSS (nt)")
    ax.set_ylabel("base frequency")
    ax.legend(loc="upper right", ncol=4, fontsize=6, handlelength=0.9, columnspacing=0.7)
    ax.text(-0.18, 1.10, "b", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- c: ATGGGG / CCCCAT ----
    ax = axes[0, 2]
    win = prof[prof.offset.between(-200, 200)].copy()
    for m, colr in [("ATGGGG", PALETTE["red"]), ("CCCCAT", PALETTE["sky"])]:
        y = win[m].rolling(11, center=True, min_periods=1).sum() / 11 * 10
        ax.plot(win.offset, y, color=colr, lw=1.2, label=m)
    ax.axvline(0, ls=":", lw=0.7, c="#666666")
    ax.set_xlabel("position relative to TSS (nt)")
    ax.set_ylabel("occurrences per 10 nt")
    ax.legend(loc="upper right", fontsize=6, handlelength=1.2, frameon=False)
    ax.text(-0.18, 1.10, "c", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- d: periodicity ----
    ax = axes[1, 0]
    for m, colr in [("ACACA", PALETTE["blue"]), ("TGTGT", PALETTE["red"]), ("TATAT", PALETTE["green"])]:
        d = pd.read_csv(f"{tdir}/periodicity_{m}.tsv", sep="\t")
        y = d.ratio.rolling(3, center=True, min_periods=1).mean()
        ax.plot(d.bin_start, y, color=colr, lw=1.1, label=m)
    ax.axhline(1, ls=":", lw=0.7, c="#666666")
    ax.axvline(0, ls=":", lw=0.5, c="#BBBBBB")
    ax.set_xlim(-1000, 1000)
    ax.set_xlabel("position relative to TSS (nt)")
    ax.set_ylabel("observed / expected")
    ax.legend(loc="upper right", fontsize=6, ncol=3, columnspacing=0.8, handlelength=1.1)
    ax.text(-0.18, 1.10, "d", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- e: TSS sites per gene ----
    ax = axes[1, 1]
    cnt = tss.n_tss.value_counts().sort_index()
    total = cnt.sum()
    ax.bar(cnt.index, cnt.values / total * 100, color=PALETTE["blue"], width=0.62, lw=0)
    for x0, v in zip(cnt.index, cnt.values):
        ax.text(x0, v / total * 100 + 0.6, str(v), ha="center", fontsize=5.6)
    ax.set_xlabel("start sites per gene")
    ax.set_ylabel("% of coding genes")
    ax.set_xticks([x for x in cnt.index if x <= 9])
    ax.set_xlim(0.3, 11.6)
    ax.set_ylim(0, max(cnt.values) / total * 100 * 1.22)
    ax.text(-0.18, 1.10, "e", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- f: TSS usage class ----
    ax = axes[1, 2]
    cls_map = {
        "单一主起点(>=80%)": "single dominant",
        "一个主起点+次要起点": "dominant + minor",
        "上游起始_无同链上游基因": "upstream start",
        "上游基因通读": "read-through",
        "多个起点_无明显主峰": "no dominant site",
    }
    cls = tss["class"].map(cls_map).value_counts()
    order = ["single dominant", "dominant + minor", "upstream start",
             "no dominant site", "read-through"]
    vals = [cls.get(k, 0) for k in order]
    cols = [PALETTE["blue"], PALETTE["sky"], PALETTE["green"], PALETTE["orange"], PALETTE["grey"]]
    yy = np.arange(len(order))[::-1]
    ax.barh(yy, vals, color=cols, height=0.62, lw=0)
    for y0, n in zip(yy, vals):
        ax.text(n + 12, y0, str(n), va="center", fontsize=5.8)
    ax.set_yticks(yy)
    ax.set_yticklabels(order, fontsize=5.9)
    ax.set_xlim(0, max(vals) * 1.18)
    ax.set_xlabel("genes", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.tick_params(axis="y", length=0)
    ax.text(-0.42, 1.10, "f", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- g: dominant TSS read share ----
    ax = axes[2, 0]
    fr = tss.tss1_frac.values
    ax.hist(fr, bins=np.arange(0, 1.02, 0.05), color=PALETTE["blue"], lw=0)
    ax.axvline(np.median(fr), color="black", ls="--", lw=0.9)
    ax.set_xlabel("read share of the dominant TSS")
    ax.set_ylabel("genes")
    ax.text(-0.18, 1.10, "g", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- h: PAS signal usage ----
    ax = axes[2, 1]
    s = sig[sig.signal != "none"].sort_values("n")
    yy = np.arange(len(s))
    ax.barh(yy, s.n, color=PALETTE["blue"], height=0.7, lw=0)
    for y0, n in zip(yy, s.n):
        ax.text(n + 3, y0, str(n), va="center", fontsize=5.8)
    ax.set_yticks(yy)
    ax.set_yticklabels(s.signal, fontsize=6.4)
    ax.set_xlabel("genes", fontsize=7)
    ax.set_xlim(0, s.n.max() * 1.22)
    ax.tick_params(axis="x", labelsize=6.4)
    ax.tick_params(axis="y", length=0)
    ax.text(-0.26, 1.10, "h", transform=ax.transAxes, fontsize=10, fontweight="bold")

    # ---- i: PAS number per gene ----
    ax = axes[2, 2]
    cnt = pas.n_pas.value_counts().sort_index()
    total = cnt.sum()
    ax.bar(cnt.index, cnt.values / total * 100, color=PALETTE["purple"], width=0.62, lw=0)
    for x0, v in zip(cnt.index, cnt.values):
        ax.text(x0, v / total * 100 + 0.8, str(v), ha="center", fontsize=5.8)
    ax.set_xlabel("polyadenylation sites per gene")
    ax.set_ylabel("% of coding genes")
    ax.set_xticks(list(cnt.index))
    ax.set_ylim(0, max(cnt.values) / total * 100 * 1.2)
    ax.text(-0.18, 1.10, "i", transform=ax.transAxes, fontsize=10, fontweight="bold")

    fig.savefig(os.path.join(outdir, "fig3_utr_promoter.png"), dpi=400)
    fig.savefig(os.path.join(outdir, "fig3_utr_promoter.pdf"))
    plt.close(fig)
    print("wrote fig3_utr_promoter")


if __name__ == "__main__":
    main()
