#!/usr/bin/env python3
"""任务3：未注释转录区的性质 + 转录起始位点（TSS）分析。

TSS 判定：取每条 primary read 在转录方向上的 5′ 端（+ 链取 reference_start，
- 链取 reference_end），按基因聚合后以 ±20 bp 聚类成 TSS 峰。

通读判定：某基因有 reads 的 5′ 端落在基因起始位点上游 >200 bp，且同链上游
存在另一个基因 —— 这些 reads 多半是从上游基因的启动子起始、越过终止信号读下来的。
"""
from __future__ import annotations

import argparse
import json
import os
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam

WIN = 2000


def parse_gtf(path):
    genes = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        gid = a.get("gene_id")
        if gid:
            genes[gid] = dict(gene_id=gid, chrom=f[0], start=int(f[3]), end=int(f[4]),
                              strand=f[6], biotype=a.get("gene_biotype", ""),
                              product=a.get("product", ""))
    return pd.DataFrame(genes.values())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--cluster", type=int, default=50)
    ap.add_argument("--upstream", type=int, default=200)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    O = args.outdir

    g = parse_gtf(args.gtf)
    g["tss"] = np.where(g.strand == "+", g.start, g.end)

    # 基因按 1 kb 窗口建索引，read 只查自己跨过的窗口
    g = g.reset_index(drop=True)
    win_idx = defaultdict(list)
    gs, ge = g.start.to_numpy(), g.end.to_numpy()
    gc, gst = g.chrom.to_numpy(), g.strand.to_numpy()
    gid_arr = g.gene_id.to_numpy()
    for i in range(len(g)):
        for w in range(gs[i] // WIN, ge[i] // WIN + 1):
            win_idx[(gc[i], w)].append(i)

    pos = defaultdict(list)   # gene_id -> [5' 坐标]
    bam = pysam.AlignmentFile(args.bam, "rb")
    n = 0
    for rec in bam.fetch(until_eof=True):
        if rec.is_secondary or rec.is_supplementary or rec.mapping_quality < args.min_mapq:
            continue
        c = rec.reference_name
        s, e = rec.reference_start + 1, rec.reference_end
        cand = set()
        for w in range(s // WIN, e // WIN + 1):
            cand.update(win_idx.get((c, w), ()))
        if not cand:
            continue
        five = e if rec.is_reverse else s
        strand = "-" if rec.is_reverse else "+"
        hit = False
        for i in cand:
            if gst[i] != strand or ge[i] < s or gs[i] > e:
                continue
            pos[gid_arr[i]].append(five)
            hit = True
        if hit:
            n += 1
    bam.close()
    print(f"assigned reads: {n}, genes with reads: {len(pos)}", flush=True)

    # ---------- TSS 聚类 ----------
    rows = []
    gset = g.set_index("gene_id")
    up_same = {}   # gene -> 上游同链基因
    for c, d in g.groupby("chrom"):
        d = d.sort_values("start")
        plus = d[d.strand == "+"].reset_index(drop=True)
        for i, r in plus.iterrows():
            prev = plus[plus.end < r.start]
            if len(prev):
                p = prev.iloc[-1]
                if r.start - p.end < 2000:
                    up_same[r.gene_id] = p.gene_id

    for gid, ps in pos.items():
        if gid not in gset.index:
            continue
        row = gset.loc[gid]
        ps = np.array(sorted(ps))
        # 聚类
        cl, cur = [], [ps[0]]
        for p in ps[1:]:
            if p - cur[-1] <= args.cluster:
                cur.append(p)
            else:
                cl.append(cur)
                cur = [p]
        cl.append(cur)
        total = len(ps)
        # 只保留占该基因 reads >=3% 的簇，过滤 5′ 端异质性造成的碎片
        cl = [c for c in cl if len(c) >= max(3, 0.03 * total)]
        cl = sorted(cl, key=len, reverse=True)
        if not cl:
            cl = [sorted(ps)]
        top = cl[0]
        tss1 = int(np.median(top))
        # 上游 >200 bp 的 5′ 端比例（相对基因自身 TSS）
        up_frac = float((ps < row.tss - args.upstream).mean()) if row.strand == "+" \
            else float((ps > row.tss + args.upstream).mean())
        rec = dict(gene_id=gid, chrom=row.chrom, strand=row.strand, biotype=row.biotype,
                   gene_tss=int(row.tss), n_reads=total, n_tss=len(cl),
                   tss1=tss1, tss1_n=len(top), tss1_frac=round(len(top) / total, 3),
                   tss1_offset=tss1 - int(row.tss),
                   upstream_frac=round(up_frac, 3),
                   upstream_gene=up_same.get(gid, ""))
        if len(cl) > 1:
            sec = cl[1]
            rec.update(tss2=int(np.median(sec)), tss2_n=len(sec),
                       tss2_frac=round(len(sec) / total, 3),
                       tss2_offset=int(np.median(sec)) - int(row.tss))
        rows.append(rec)
    tss = pd.DataFrame(rows)

    # ---------- 分类 ----------
    def cls(r):
        if r.upstream_frac >= 0.3 and r.upstream_gene:
            return "上游基因通读"
        if r.upstream_frac >= 0.3:
            return "上游起始_无同链上游基因"
        if r.tss1_frac >= 0.8:
            return "单一主起点(>=80%)"
        if r.tss1_frac >= 0.5:
            return "一个主起点+次要起点"
        return "多个起点_无明显主峰"

    tss["class"] = tss.apply(cls, axis=1)
    tss = tss.sort_values("n_reads", ascending=False)
    tss.to_csv(f"{O}/tss_per_gene.tsv", sep="\t", index=False)

    pc = tss[tss.biotype == "protein_coding"]
    summary = dict(
        genes_with_reads=int(len(tss)),
        pc_genes_with_reads=int(len(pc)),
        single_tss=int((pc.n_tss == 1).sum()),
        two_tss=int((pc.n_tss == 2).sum()),
        multi_tss=int((pc.n_tss >= 3).sum()),
        primary_tss_within_100bp=int((pc.tss1_offset.abs() <= 100).sum()),
        primary_tss_within_500bp=int((pc.tss1_offset.abs() <= 500).sum()),
        readthrough=int((pc["class"] == "上游基因通读").sum()),
        upstream_no_gene=int((pc["class"] == "上游起始_无同链上游基因").sum()),
        median_primary_frac=round(float(pc.tss1_frac.median()), 3),
        median_n_tss=float(pc.n_tss.median()),
    )
    json.dump(summary, open(f"{O}/tss_summary.json", "w"), indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print()
    print(pc["class"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
