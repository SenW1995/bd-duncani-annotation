#!/usr/bin/env python3
"""Paper 1（Scientific Data）图 1–2

Fig 1  组装：染色体长度与基因密度、三套组装比较、BUSCO、与两套已发表组装的共线性
Fig 2  注释：基因类型、外显子数、UTR 长度、内含子长度、转录支持、每基因 reads、总览

用法（WSL bioinfo 环境）
  python3 fig_paper1_resource.py --root <GEN_ROOT> --out <GEN_ROOT>/papers/figs
"""
import argparse
import collections
import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch, PathPatch, Rectangle

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

C = dict(main="#2f6f9f", alt="#c0504d", grey="#8c8c8c", light="#dbe7f2",
         green="#4f8a5b", orange="#e08a3c", purple="#8064a2", red="#b03a2e")
A = re.compile(r'(\w+) "([^"]*)"')


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def panel(ax, letter, x=-0.16, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def _read_paf(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        r = dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                 t=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8]))
        r["block"] = r["qe"] - r["qs"]
        rows.append(r)
    return rows


def _layout(order, lengths, gap_frac=0.008):
    """按序列长度比例分配横轴，返回 {name: (x0, x1)}"""
    total = sum(lengths[n] for n in order)
    usable = 1 - gap_frac * (len(order) - 1)
    off, cur = {}, 0.0
    for n in order:
        w = usable * lengths[n] / total
        off[n] = (cur, cur + w)
        cur += w + gap_frac
    return off


def _ref_track(paf, small_bin=200_000):
    """参考组装：长序列单列，短序列并进一个 other 块"""
    L = {}
    for r in paf:
        L.setdefault(r["t"], r["tlen"])
    big = sorted([t for t in L if L[t] >= small_bin], key=lambda t: -L[t])
    small = [t for t in L if L[t] < small_bin]
    lens = {t: L[t] for t in big}
    if small:
        lens["other"] = int(sum(L[t] for t in small))
    return big, small, lens


def _draw_synteny(ax, paf_s, paf_w, ours, olen, min_block=2000):
    """三轨丝带图：上 = Singh 2023，中 = 本研究，下 = Wang 2022。"""
    s_big, s_small, s_len = _ref_track(paf_s)
    w_big, w_small, w_len = _ref_track(paf_w)

    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(0, 1.02)
    ax.axis("off")

    y_top, y_mid, y_bot = 0.80, 0.50, 0.19
    h = 0.030

    o_off = _layout(ours, olen)
    s_off = _layout(s_big + (["other"] if s_small else []), s_len)
    w_off = _layout(w_big + (["other"] if w_small else []), w_len)
    HALO = [pe.withStroke(linewidth=2.0, foreground="white")]

    def draw_track(off, y, color, labels, fs=6.2, above=True):
        items = sorted(off.items(), key=lambda kv: kv[1][0])
        k = 0
        for n, (x0, x1) in items:
            ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, fc=color, ec="none"))
            if n not in labels:
                continue
            w = x1 - x0
            if w >= 0.055:
                ax.text((x0 + x1) / 2, y + (h * 1.7 if above else -h * 2.0),
                        labels[n], ha="center", va="bottom" if above else "top",
                        fontsize=fs, color="#222222", path_effects=HALO)
            else:
                dy = (0.030 if k % 2 == 0 else 0.075) * (1 if above else -1)
                ax.text((x0 + x1) / 2, y + dy, labels[n], ha="center",
                        va="bottom" if above else "top", fontsize=fs - 0.4,
                        rotation=90, color="#222222", path_effects=HALO)
                k += 1

    draw_track(s_off, y_top, "#9a9a9a",
               {t: f"chr {i+1}" for i, t in enumerate(s_big)} | {"other": "71 contigs"},
               fs=6.0, above=True)
    draw_track(o_off, y_mid, C["main"],
               {q: (f"chr {i+1}" if i < 5 else q.replace("scaffold_", "scaf "))
                for i, q in enumerate(ours)},
               fs=6.4, above=True)
    draw_track(w_off, y_bot, "#9a9a9a",
               {t: f"chr {i+1}" for i, t in enumerate(w_big)}, fs=6.0, above=False)

    def dominant_flip(paf):
        f = {}
        for t in {r["t"] for r in paf}:
            fwd = sum(r["block"] for r in paf if r["t"] == t and r["strand"] == "+"
                      and r["block"] >= min_block)
            rev = sum(r["block"] for r in paf if r["t"] == t and r["strand"] == "-"
                      and r["block"] >= min_block)
            f[t] = rev > fwd
        return f

    flip_s = dominant_flip(paf_s)
    flip_w = dominant_flip(paf_w)

    def frac(off, n, pos, length, flipped=False):
        x0, x1 = off[n]
        if flipped:
            pos = max(length, 1) - pos
        return x0 + (x1 - x0) * (pos / max(length, 1))

    n_rev = 0

    def ribbon(b, off_ref, ref_len, ref_small, y_ref, color, flip):
        t = b["t"] if b["t"] in off_ref else ("other" if b["t"] in ref_small else None)
        if t is None:
            return
        fl = flip.get(b["t"], False)
        xa = frac(o_off, b["q"], b["qs"], olen[b["q"]])
        xb = frac(o_off, b["q"], b["qe"], olen[b["q"]])
        ya = frac(off_ref, t, b["ts"], ref_len[t], fl)
        yb = frac(off_ref, t, b["te"], ref_len[t], fl)
        verts = [(xa, y_mid), (xa, (y_mid + y_ref) / 2), (ya, (y_mid + y_ref) / 2),
                 (ya, y_ref), (yb, y_ref), (yb, (y_mid + y_ref) / 2),
                 (xb, (y_mid + y_ref) / 2), (xb, y_mid), (xa, y_mid)]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=0.42, lw=0))

    for b in paf_s:
        if b["block"] < min_block:
            continue
        inv = (b["strand"] == "-") != flip_s.get(b["t"], False)
        n_rev += inv
        ribbon(b, s_off, s_len, s_small, y_top, C["orange"] if inv else C["main"], flip_s)
    for b in paf_w:
        if b["block"] < min_block:
            continue
        inv = (b["strand"] == "-") != flip_w.get(b["t"], False)
        n_rev += inv
        ribbon(b, w_off, w_len, w_small, y_bot, C["orange"] if inv else C["main"], flip_w)

    ax.text(-0.012, y_top, "Singh 2023", ha="right", va="center", fontsize=6.6, rotation=90)
    ax.text(-0.012, y_mid, "this study", ha="right", va="center", fontsize=6.6,
            rotation=90, color=C["main"], fontweight="bold")
    ax.text(-0.012, y_bot, "Wang 2022", ha="right", va="center", fontsize=6.6, rotation=90)

    ax.plot([], [], color=C["main"], lw=5, alpha=.45, label="collinear, same strand")
    ax.plot([], [], color=C["orange"], lw=5, alpha=.45, label="collinear, inverted")
    ax.legend(fontsize=6.0, loc="upper right", bbox_to_anchor=(1.0, 1.03),
              ncol=2, frameon=False, handlelength=1.2, handletextpad=0.6,
              columnspacing=1.6)
    return n_rev


def fig1(root, out):
    fa = read_fasta(os.path.join(root, "V2/babesia_duncani_no_apicoplast.fa"))
    lens = sorted(((len(v), k) for k, v in fa.items()), reverse=True)
    gtf = os.path.join(root, "papers/Submission_Package/02_data/annotation/babesia_duncani_v3.10.gtf")

    genes = collections.defaultdict(list)
    offsets = {}
    cum = 0
    for _, name in lens:
        offsets[name] = cum
        cum += len(fa[name])
    for line in open(gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(A.findall(f[8]))
        if a.get("gene_biotype") == "protein_coding":
            genes[f[0]].append((int(f[3]) + int(f[4])) // 2)


    fig = plt.figure(figsize=(7.2, 6.4))

    # ---------------- (a) 滑窗 AT 含量 ----------------
    ax = fig.add_axes([0.095, 0.575, 0.885, 0.395])
    win, step = 50_000, 10_000
    at_lo, at_hi = 50.0, 78.0     # 标尺下限/上限（%）
    xmax = lens[0][0] / 1e6
    ylabels = []
    for i, (L, name) in enumerate(lens):
        y = len(lens) - 1 - i
        seq = fa[name]
        short = name.replace("scaffold_", "scaf ")
        head = f"chr {i+1}" if i < 5 else short
        ylabels.append(head)
        # 端部 10% 用浅灰底标出
        for x0 in (0.0, L * 0.9 / 1e6):
            ax.add_patch(Rectangle((x0, y - 0.05), L * 0.1 / 1e6, 0.94,
                                   fc="#e9e9e9", ec="none", zorder=0))
        if L <= win:
            at = (seq.count("A") + seq.count("T")) / len(seq) * 100
            xs = np.array([L / 2e6])
            ys = np.array([at])
        else:
            xs = (np.arange(0, L - win + 1, step) + win / 2) / 1e6
            ys = np.array([
                (seq[p:p + win].count("A") + seq[p:p + win].count("T")) / win * 100
                for p in np.arange(0, L - win + 1, step)
            ])
        hgt = np.clip((ys - at_lo) / (at_hi - at_lo), 0, 1) * 0.86
        ax.plot(xs, y + hgt, color=C["main"], lw=0.55, zorder=3)
        ax.fill_between(xs, y, y + hgt, color=C["main"], alpha=0.38, lw=0, zorder=2)
        ax.text(L / 1e6 + 0.06, y + 0.30,
                f"{L/1e6:.2f} Mb · {len(fa[name]):,} bp",
                va="center", fontsize=5.7, color="#444444")

    ax.set_yticks(range(len(lens)))
    ax.set_yticklabels(ylabels[::-1], fontsize=6.4)
    ax.set_xlim(0, xmax * 1.24)
    ax.set_ylim(-0.12, len(lens) - 0.28)
    ax.set_xlabel("Position along sequence (Mb)", fontsize=7.2)
    ax.tick_params(axis="x", labelsize=6.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.text(0.995, 0.06, f"AT content, {win//1000}-kb windows ({step//1000}-kb step)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.7,
            color="#555555")
    ax.axhline(0, color="#cccccc", lw=0.4, zorder=1)
    panel(ax, "a", x=-0.20, y=1.14)

    # ---------------- (b) 三轨共线性 ----------------
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.465])
    paf_s = _read_paf(os.path.join(root, "annotation_qc/genome_compare/ours_vs_singh.paf"))
    paf_w = _read_paf(os.path.join(root, "annotation_qc/genome_compare/ours_vs_wang.paf"))
    ours = sorted({r["q"] for r in paf_s}, key=lambda s: int(s.split("_")[1]))
    olen = {q: len(fa[q]) for q in ours}
    _draw_synteny(ax, paf_s, paf_w, ours, olen)
    panel(ax, "b", x=-0.045, y=1.00)

    fig.savefig(os.path.join(out, "fig1_assembly.png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(out, "fig1_assembly.pdf"), bbox_inches="tight")
    plt.close(fig)


def fig2(root, out):
    gtf = os.path.join(root, "papers/Submission_Package/02_data/annotation/babesia_duncani_v3.10.gtf")
    biotype, exon_n, u5, u3, intr = collections.Counter(), collections.Counter(), [], [], []
    last_end, last_chrom = {}, {}
    for line in open(gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(A.findall(f[8]))
        g = a.get("gene_id")
        if f[2] == "gene":
            biotype[a.get("gene_biotype")] += 1
        elif f[2] == "exon":
            exon_n[g] += 1
            s, e = int(f[3]), int(f[4])
            if g in last_end and last_chrom[g] == f[0]:
                intr.append(abs(s - last_end[g]) - 1)
            last_end[g] = e
            last_chrom[g] = f[0]
        elif f[2] == "five_prime_utr":
            u5.append(int(f[4]) - int(f[3]) + 1)
        elif f[2] == "three_prime_utr":
            u3.append(int(f[4]) - int(f[3]) + 1)

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.34,
                          left=0.10, right=0.97, top=0.94, bottom=0.11)

    ax = fig.add_subplot(gs[0, 0])
    vals = list(exon_n.values())
    bins = np.arange(0.5, min(20, max(vals)) + 1.5, 1)
    ax.hist(vals, bins=bins, color=C["main"], edgecolor="white", lw=0.4)
    ax.set_xlabel("Exons per gene", fontsize=7.2)
    ax.set_ylabel("Genes", fontsize=7.2)
    ax.set_xlim(0.4, 15.5)
    ax.tick_params(labelsize=6.6)
    ax.text(0.55, 0.82, f"median {int(np.median(vals))}", transform=ax.transAxes, fontsize=6.5)
    panel(ax, "a", x=-0.24)

    ax = fig.add_subplot(gs[0, 1])
    # 与 introns.tsv 保持一致：只统计 >=20 bp 的内含子（fl_validate.py 的阈值）。
    # 4 bp 的 BdWA1_002g00742 GTAG 是 internal-stop 拆分假象，不是真实内含子。
    iv = np.array([x for x in intr if x >= 20])
    bins = np.logspace(0.8, 4.3, 40)
    ax.hist(iv, bins=bins, color=C["purple"], lw=0)
    ax.set_xscale("log")
    ax.axvline(np.median(iv), color="black", ls="--", lw=0.8)
    ax.set_xlabel("Intron length (bp, log)", fontsize=7.2)
    ax.set_ylabel("Introns", fontsize=7.2)
    ax.tick_params(labelsize=6.6)
    ax.text(0.52, 0.62, f"n = {len(iv):,}\nmedian {int(np.median(iv))} bp",
            transform=ax.transAxes, fontsize=6.2, ha="left")
    panel(ax, "b")

    ax = fig.add_subplot(gs[1, 0])
    coding_ids = set()
    for line in open(gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "gene":
            a0 = dict(A.findall(f[8]))
            if a0.get("gene_biotype") == "protein_coding":
                coding_ids.add(a0["gene_id"])
    sup = [r for r in csv.DictReader(open(os.path.join(
        root, "papers/Submission_Package/02_data/tables/01_gene_support.tsv"), encoding="utf-8"),
        delimiter="\t") if r["gene_id"] in coding_ids]
    n = len(sup)
    anyread = sum(1 for r in sup if int(r["n_reads"]) > 0)
    alle = sum(1 for r in sup if float(r["frac_exons_hit"]) >= 1.0)
    ten = sum(1 for r in sup if float(r["frac_exon_bp_10x"]) >= 1.0)
    bars = [("Any read", anyread, C["main"]),
            ("All exons", alle, C["green"]),
            ("≥10× every exon", ten, C["orange"])]
    yy = np.arange(len(bars))[::-1]
    ax.barh(yy, [b[1] / n * 100 for b in bars], color=[b[2] for b in bars], height=0.55)
    for y0, b in zip(yy, bars):
        ax.text(b[1] / n * 100 + 0.8, y0, f"{b[1]/n*100:.1f}%", va="center", fontsize=6.5)
    ax.set_yticks(yy)
    ax.set_yticklabels([b[0] for b in bars], fontsize=6.5)
    ax.set_xlim(0, 112)
    ax.set_xlabel(f"% of {n:,} coding genes", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.6)
    panel(ax, "c", x=-0.55)

    ax = fig.add_subplot(gs[1, 1])
    nr = np.array([int(r["n_reads"]) for r in sup if int(r["n_reads"]) > 0])
    ax.hist(np.log10(nr + 1), bins=40, color=C["main"], lw=0)
    ax.set_xlabel("log10(reads per gene + 1)", fontsize=7)
    ax.set_ylabel("Genes", fontsize=7.2)
    ax.axvline(np.log10(11), color=C["red"], ls="--", lw=0.8)
    ax.text(np.log10(11) + 0.02, 0.86, "10 reads", color=C["red"],
            fontsize=5.8, ha="left", transform=ax.get_xaxis_transform())
    ax.tick_params(labelsize=6.6)
    ax.text(0.97, 0.84, f"n = {len(nr):,}", transform=ax.transAxes, fontsize=6.3, ha="right")
    panel(ax, "d", x=-0.32)

    fig.savefig(os.path.join(out, "fig2_annotation.png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(out, "fig2_annotation.pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig1(a.root, a.out)
    fig2(a.root, a.out)
    print("Paper1 图完成 ->", a.out)


if __name__ == "__main__":
    main()
