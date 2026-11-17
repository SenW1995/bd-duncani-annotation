"""论文级图表统一风格。

用法：
    from pub_style import *
    fig, ax = plt.subplots(figsize=MM(85, 60))
    ...
    save(fig, "fig1_overview")
输出 PDF（矢量）+ PNG 600 dpi，尺寸按期刊常用栏宽。
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ---------- 尺寸：期刊常用 ----------
MM = lambda w, h=None: (w / 25.4, (h if h else w * 0.75) / 25.4)
SINGLE = 85      # 单栏 mm
DOUBLE = 170     # 双栏 mm

# ---------- 配色：色盲友好（Okabe-Ito 变体） ----------
PALETTE = {
    "blue":    "#0072B2",
    "orange":  "#E69F00",
    "green":   "#009E73",
    "red":     "#D55E00",
    "purple":  "#CC79A7",
    "sky":     "#56B4E9",
    "yellow":  "#F0E442",
    "grey":    "#999999",
    "dark":    "#333333",
}
SEQ = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"], PALETTE["red"],
       PALETTE["purple"], PALETTE["sky"], PALETTE["yellow"], PALETTE["grey"]]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,       # 字体嵌入，期刊要求
    "ps.fonttype": 42,
    "axes.prop_cycle": plt.cycler(color=SEQ),
})

OUTDIR = os.environ.get("PUB_FIG_DIR", "./figs_pub")


def panel(ax, letter: str, x: float = -0.16, y: float = 1.14):
    """面板标号：放在标题上方，避免与标题、坐标轴压字"""
    if not letter:
        return
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="bottom", ha="left")


# 可选：本机绘图辅助脚本（非复现必需）；如需启用，设置环境变量 SKILL_SCRIPTS
SKILL_SCRIPTS = os.environ.get("SKILL_SCRIPTS", "")


def equalize_panels(fig) -> None:
    """单行多面板：重排成等宽、等间距的绘图区

    tight_layout 会按每个面板的刻度标签宽度各自收缩，导致宽度不等；
    这里取最窄者作为统一宽度，在原有总跨度内等间距重排。
    """
    axes = [a for a in fig.axes if a.get_visible()]
    if len(axes) < 2:
        return
    if len({round(a.get_position().y0, 3) for a in axes}) != 1:
        return                                    # 多行布局不在此处理
    boxes = [a.get_position() for a in axes]
    wmax = max(b.width for b in boxes)
    # 色条等窄轴不参与等宽：否则会被一起压扁，标签反而更挤
    keep = [i for i, b in enumerate(boxes) if b.width >= 0.4 * wmax]
    axes = [axes[i] for i in keep]
    boxes = [boxes[i] for i in keep]
    n = len(boxes)
    if n < 2:
        return
    left = min(b.x0 for b in boxes)
    right = max(b.x1 for b in boxes)
    w = min(b.width for b in boxes)
    if n > 1 and (right - left - n * w) / (n - 1) < 0.008:
        w = (right - left) / n                    # 间距太挤时直接等分
    gap = (right - left - n * w) / (n - 1) if n > 1 else 0
    for k, a in enumerate(axes):
        y0, h = a.get_position().y0, a.get_position().height
        a.set_position([left + k * (w + gap), y0, w, h])


def align_gate(fig, name: str, outdir: str) -> None:
    """多面板图必须过等宽/对齐门（nature-figure 的阻断式检查）"""
    if len(fig.axes) < 2:
        return
    import sys
    if SKILL_SCRIPTS not in sys.path:
        sys.path.insert(0, SKILL_SCRIPTS)
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(
            fig, json_out=os.path.join(outdir, f"{name}.alignment.json"))
        print(f"  alignment ok: {name}")
    except Exception as exc:                      # 报告但不静默：失败会打印全文
        print(f"  ALIGNMENT FAILED: {name}: {exc}")


def save(fig, name: str, outdir: str | None = None, gate: bool = True):
    d = outdir or OUTDIR
    os.makedirs(d, exist_ok=True)
    if gate:
        equalize_panels(fig)
        align_gate(fig, name, d)
    for ext in ("pdf", "png"):
        # bbox_inches="tight" 保证面板标号不被裁掉
        fig.savefig(os.path.join(d, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"saved {name}.pdf/.png")
