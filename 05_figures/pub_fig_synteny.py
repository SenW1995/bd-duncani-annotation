#!/usr/bin/env python3
"""Fig. 1 三轨共线性图：上 = Singh 2023，中 = 本研究，下 = Wang 2022

做法沿用工作区 genome_synteny.py 的丝带图（贝塞尔曲线连接两条轨道），
扩成三条轨道：参考序列在上、本研究在中间、第二套参考在下。
组装统计（总长 / 序列数 / N50 / 覆盖度）不再画条形图，改为 Table 1。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle

from pub_style import DOUBLE, MM, PALETTE, plt, save

C = PALETTE
D = "<GEN_ROOT>/annotation_qc/genome_compare"
OUT = "<WORK_DIR>/figs_pub"

MIN_BLOCK = 2000          # 只画 ≥2 kb 的共线块
SMALL_BIN = 200_000       # 参考组装里短于 200 kb 的序列合并成一个灰块


def read_paf(path):
    p = pd.read_csv(path, sep="\t", header=None, comment="#",
                    usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                    names=["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te"])
    p["block"] = p.qe - p.qs
    return p


def track_layout(order, lengths, gap_frac=0.008):
    """按序列长度比例分配横轴，返回 {name: (x0, x1)}"""
    total = sum(lengths[n] for n in order)
    gap = gap_frac
    usable = 1 - gap * (len(order) - 1)
    off, cur = {}, 0.0
    for n in order:
        w = usable * lengths[n] / total
        off[n] = (cur, cur + w)
        cur += w + gap
    return off


def main():
    paf_s = read_paf(f"{D}/ours_vs_singh.paf")
    paf_w = read_paf(f"{D}/ours_vs_wang.paf")

    # 本研究：五条染色体 + 两条小 scaffold，按编号排序
    ours = sorted(paf_s.q.unique(), key=lambda s: int(s.split("_")[1]))
    olen = {q: int(paf_s[paf_s.q == q].qlen.iloc[0]) for q in ours}

    def ref_track(paf):
        """参考组装轨道：大序列单列，短序列并进一个 other 块"""
        L = paf.drop_duplicates("t").set_index("t").tlen.astype(int)
        big = [t for t in L.index if L[t] >= SMALL_BIN]
        big.sort(key=lambda t: -L[t])
        small = [t for t in L.index if L[t] < SMALL_BIN]
        lens = {t: L[t] for t in big}
        if small:
            lens["other"] = int(sum(L[t] for t in small))
        return big, small, lens

    s_big, s_small, s_len = ref_track(paf_s)
    w_big, w_small, w_len = ref_track(paf_w)

    fig, ax = plt.subplots(figsize=MM(DOUBLE, 88))
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(0, 1.02)
    ax.axis("off")

    y_top, y_mid, y_bot = 0.80, 0.50, 0.19
    h = 0.028

    o_off = track_layout(ours, olen)
    s_off = track_layout(s_big + (["other"] if s_small else []), s_len)
    w_off = track_layout(w_big + (["other"] if w_small else []), w_len)

    HALO = [pe.withStroke(linewidth=2.0, foreground="white")]   # 标签压在丝带上也读得清

    def draw_track(off, y, color, labels, fs=6.2, above=True):
        """窄序列的标签竖排并错层，避免相互压住"""
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

    # 轨道底色：参考用灰、本研究用深蓝，突出中心
    draw_track(s_off, y_top, "#9a9a9a",
               {t: f"chr {i+1}" for i, t in enumerate(s_big)} | {"other": "71 contigs"},
               fs=6.0, above=True)
    draw_track(o_off, y_mid, C["blue"],
               {q: (f"chr {i+1}" if i < 5 else q.replace("scaffold_", "scaf ").replace("_", ""))
                for i, q in enumerate(ours)},
               fs=6.4, above=True)
    draw_track(w_off, y_bot, "#9a9a9a",
               {t: f"chr {i+1}" for i, t in enumerate(w_big)}, fs=6.0, above=False)

    # 参考序列按"主导链方向"整体翻转：否则同一条染色体在参考里是反向的，
    # 丝带会满屏交叉，读不出共线性
    def dominant_flip(paf):
        f = {}
        for t, g in paf[paf.block >= MIN_BLOCK].groupby("t"):
            fwd = g[g.strand == "+"].block.sum()
            rev = g[g.strand == "-"].block.sum()
            f[t] = rev > fwd
        return f

    flip_s = dominant_flip(paf_s)
    flip_w = dominant_flip(paf_w)

    def frac(off, n, pos, length, flipped=False):
        x0, x1 = off[n]
        if flipped:
            pos = max(length, 1) - pos
        return x0 + (x1 - x0) * (pos / max(length, 1))

    def ribbon(b, off_ref, ref_len, ref_small, y_ref, color, flip):
        t = b["t"] if b["t"] in off_ref else ("other" if b["t"] in ref_small else None)
        if t is None:
            return
        fl = flip.get(b["t"], False)
        xa = frac(o_off, b["q"], b["qs"], olen[b["q"]])
        xb = frac(o_off, b["q"], b["qe"], olen[b["q"]])
        ya = frac(off_ref, t, b["ts"], ref_len[t], fl)
        yb = frac(off_ref, t, b["te"], ref_len[t], fl)
        verts = [(xa, y_mid), (xa, (y_mid + y_ref) / 2), (ya, (y_mid + y_ref) / 2), (ya, y_ref),
                 (yb, y_ref), (yb, (y_mid + y_ref) / 2), (xb, (y_mid + y_ref) / 2), (xb, y_mid),
                 (xa, y_mid)]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=0.42, lw=0))

    n_rev = 0
    for _, b in paf_s[paf_s.block >= MIN_BLOCK].iterrows():
        inv = (b.strand == "-") != flip_s.get(b["t"], False)
        n_rev += inv
        ribbon(b, s_off, s_len, s_small, y_top, C["orange"] if inv else C["blue"], flip_s)
    for _, b in paf_w[paf_w.block >= MIN_BLOCK].iterrows():
        inv = (b.strand == "-") != flip_w.get(b["t"], False)
        n_rev += inv
        ribbon(b, w_off, w_len, w_small, y_bot, C["orange"] if inv else C["blue"], flip_w)

    # 轨道名与覆盖度
    ax.text(-0.012, y_top, "Singh 2023", ha="right", va="center", fontsize=6.6, rotation=90)
    ax.text(-0.012, y_mid, "this study", ha="right", va="center", fontsize=6.6, rotation=90,
            color=C["blue"], fontweight="bold")
    ax.text(-0.012, y_bot, "Wang 2022", ha="right", va="center", fontsize=6.6, rotation=90)

    ax.plot([], [], color=C["blue"], lw=5, alpha=.45, label="collinear, same strand")
    ax.plot([], [], color=C["orange"], lw=5, alpha=.45, label="collinear, inverted")
    ax.legend(fontsize=6.0, loc="upper right", bbox_to_anchor=(1.0, 1.02),
              ncol=2, frameon=False, handlelength=1.2, handletextpad=0.6,
              columnspacing=1.6)

    fig.tight_layout()
    save(fig, "fig1_synteny_three_track", OUT)
    print("blocks drawn", int((paf_s.block >= MIN_BLOCK).sum()), int((paf_w.block >= MIN_BLOCK).sum()),
          "inverted", int(n_rev))


if __name__ == "__main__":
    main()
