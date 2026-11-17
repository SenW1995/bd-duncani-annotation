#!/usr/bin/env python3
"""Paper 1（Scientific Data）图 5：结构层与结构注释

(a) 全部 3,616 个 Boltz-2 模型的 mean pLDDT 分布（灰柱），叠加进入结构命名的
    784 个模型（蓝线，按各自组内百分比归一）；点线为 50 / 70 置信档界
(b) 58 个采纳基因：mean pLDDT vs Foldseek 比对覆盖度，按邻居名类型着色

数据输入
  --tables/plddt_all_models.tsv          全量模型 pLDDT（srv_56_all_plddt.py 生成）
  --root/V2/v4_structure/structure_confidence_all.tsv   784 个被评估模型与 58 个采纳

用法（WSL bioinfo 环境）
  python3 fig_paper1_structure.py --root <GEN_ROOT> --out <GEN_ROOT>/papers/Paper1_ScientificData/figs
"""
import argparse
import collections
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

C = dict(high="#1f4e79", ok="#5b9bd5", dom="#e08a3c", low="#c9c9c9",
         accept="#3f7f4f", reject="#9a9a9a", alt="#c0504d")


def panel(ax, letter, x=-0.20, y=1.12):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def bucket(reason, tier=""):
    r = reason.strip()
    if tier == "no_structure":
        return "No searchable model"
    if not r:
        return "Accepted"
    if r == "无结构文件":
        return "No searchable model"
    if r.startswith("序列已给出"):
        return "Sequence already named the domain"
    if r.startswith("邻居名无功能"):
        return "Neighbour name uninformative"
    if r.startswith("prob<"):
        return "Probability < 0.8"
    if r.startswith("evalue>"):
        return "E-value > 1e-5"
    if "未命名" in r:
        return "Neighbour itself unnamed"
    if "无命中" in r:
        return "No hit"
    if "已有功能名" in r:
        return "Already had a product name"
    if "置信度不足" in r:
        return "Model confidence insufficient"
    if "域支撑" in r or "结构域不足" in r:
        return "Domain coverage insufficient"
    if "一致性" in r:
        return "Sequence identity < 95%"
    return "Other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tables", default=None,
                    help="含 plddt_all_models.tsv 的目录，默认 <out>/../tables")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    tdir = a.tables or os.path.join(os.path.dirname(os.path.abspath(a.out)), "tables")
    all_models = pd.read_csv(os.path.join(tdir, "plddt_all_models.tsv"), sep="\t")

    rows = list(csv.DictReader(open(os.path.join(
        a.root, "V2/v4_structure/structure_confidence_all.tsv"), encoding="utf-8"),
        delimiter="\t"))
    for r in rows:
        r["plddt"] = float(r["mean_pLDDT"]) if r["mean_pLDDT"] else np.nan
        r["prob_f"] = float(r["prob"]) if r["prob"] not in ("", None) else np.nan

    fig = plt.figure(figsize=(7.2, 2.7))

    # ---------------- (a) pLDDT 分布：全部模型 vs 进入命名的子集 ----------------
    ax = fig.add_axes([0.075, 0.21, 0.395, 0.66])
    bins = np.arange(20, 100, 2.5)
    allv = all_models["mean_pLDDT"].values.astype(float)
    subv = np.array([r["plddt"] for r in rows if not np.isnan(r["plddt"])])
    ax.hist(allv, bins=bins, weights=np.full(allv.size, 100.0 / allv.size),
            color="#d6d6d6", lw=0, label=f"all models (n = {allv.size:,})")
    h_sub, _ = np.histogram(subv, bins=bins, weights=np.full(subv.size, 100.0 / subv.size))
    ax.step(bins[:-1], h_sub, where="post", color=C["high"], lw=1.1,
            label=f"assessed for naming (n = {subv.size})")
    top = max((100.0 / allv.size * np.histogram(allv, bins=bins)[0]).max(),
              (100.0 / subv.size * np.histogram(subv, bins=bins)[0]).max())
    ax.set_ylim(0, top * 1.22)
    for x in (50, 70):
        ax.axvline(x, ls=":", lw=0.7, c="#777777", zorder=0)
    ax.set_xlabel("Mean pLDDT", fontsize=7.2)
    ax.set_ylabel("% of models in group", fontsize=7.2)
    ax.set_xlim(20, 98)
    ax.tick_params(labelsize=6.6)
    leg = ax.legend(frameon=True, fontsize=5.5, loc="upper right", handlelength=1.1,
                    labelspacing=0.25, borderpad=0.3, handletextpad=0.5)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("none")
    leg.get_frame().set_alpha(0.92)
    panel(ax, "a", x=-0.26, y=1.14)

    # ---------------- (b) 采纳基因 ----------------
    ax = fig.add_axes([0.615, 0.21, 0.335, 0.66])
    acc = [r for r in rows if r["decision"] == "accept"]
    styles = {"specific": (C["high"], "o", "specific name"),
              "domain_only": (C["dom"], "^", "domain-level name")}
    for q, (col, mk, lab) in styles.items():
        sub = [r for r in acc if r["name_quality"] == q]
        ax.scatter([r["plddt"] for r in sub],
                   [float(r["aln_cov"]) if r["aln_cov"] else np.nan for r in sub],
                   s=17, color=col, marker=mk, lw=0.3, edgecolor="white",
                   alpha=0.9, label=f"{lab} (n={len(sub)})", zorder=3)
    ax.axvline(70, color="#bbbbbb", ls="--", lw=0.6, zorder=1)
    ax.axvline(50, color="#bbbbbb", ls="--", lw=0.6, zorder=1)
    ax.axhline(100, color="#888888", ls=":", lw=0.6, zorder=1)
    ax.text(71, 12, "70", fontsize=5.2, color="#777777")
    ax.text(51, 12, "50", fontsize=5.2, color="#777777")
    ax.set_xlim(45, 92)
    ax.set_ylim(0, 130)
    ax.set_xlabel("Mean pLDDT of the model", fontsize=7.2)
    ax.set_ylabel("Aligned coverage of the protein (%)", fontsize=7.2)
    ax.tick_params(labelsize=6.6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], linestyle="none", marker="o", markerfacecolor=C["high"],
                      markeredgecolor="none", markersize=5, label="specific name (n=13)"),
               Line2D([], [], linestyle="none", marker="^", markerfacecolor=C["dom"],
                      markeredgecolor="none", markersize=5, label="domain-level name (n=45)")]
    ax.legend(handles=handles, frameon=False, fontsize=5.6, loc="lower right",
              handletextpad=0.9, borderpad=0.1, labelspacing=0.4)
    panel(ax, "b", x=-0.30, y=1.14)

    fig.savefig(os.path.join(a.out, "fig5_structure.png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(a.out, "fig5_structure.pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Fig 5 完成 -> {a.out}")


if __name__ == "__main__":
    main()
