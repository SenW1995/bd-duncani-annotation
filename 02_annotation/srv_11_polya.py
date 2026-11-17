#!/usr/bin/env python3
"""多聚腺苷化位点（PAS）与可变多聚腺苷化（APA）分析。

原理：全长 cDNA reads 的 3′ 端对应 polyA 加尾位点。
  1. 每条 read 取转录方向的 3′ 端，按基因聚合
  2. ±25 bp 聚类成 PAS，保留占该基因 reads ≥5% 的峰
  3. 每个 PAS 上游 60 bp 搜 AATAAA 类信号，记录相对位置
  4. 按 PAS 数目与主峰占比分类 APA

输出
  pas_sites.tsv      每个 PAS 的坐标、reads 数、motif
  apa_per_gene.tsv   每基因的 PAS 数、主峰占比、3′UTR 长度变化
  apa_summary.json
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

WIN = 1000
BIN = 25
MIN_FRAC = 0.10          # 一个 PAS 至少占该基因 10% 的 reads
MIN_READS = 10
# 只用真正的 polyA 信号；TATAAA 是 TATA box，必须排除
SIGNALS = ["AATAAA", "ATTAAA", "AGTAAA", "AATACA", "AATATA", "CATAAA", "GATAAA"]
WIN_LO, WIN_HI = 10, 45   # 信号到 PAS 的距离窗口（polyA 信号典型在 -20 附近）


def parse_gtf(path):
    genes = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        if a.get("gene_id"):
            genes[a["gene_id"]] = dict(gene_id=a["gene_id"], chrom=f[0], start=int(f[3]),
                                       end=int(f[4]), strand=f[6],
                                       biotype=a.get("gene_biotype", ""))
    return pd.DataFrame(genes.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-mapq", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    g = parse_gtf(a.gtf).reset_index(drop=True)
    win = defaultdict(list)
    gs, ge = g.start.to_numpy(), g.end.to_numpy()
    for i in range(len(g)):
        for w in range(gs[i] // WIN, ge[i] // WIN + 1):
            win[(g.chrom.iat[i], w)].append(i)
    gid, gst = g.gene_id.to_numpy(), g.strand.to_numpy()
    print(f"基因 {len(g)}", flush=True)

    three = defaultdict(list)
    bam = pysam.AlignmentFile(a.bam, "rb")
    n = 0
    for rec in bam.fetch(until_eof=True):
        if rec.is_secondary or rec.is_supplementary or rec.mapping_quality < a.min_mapq:
            continue
        s, e = rec.reference_start + 1, rec.reference_end
        cand = set()
        for w in range(s // WIN, e // WIN + 1):
            cand.update(win.get((rec.reference_name, w), ()))
        if not cand:
            continue
        strand = "-" if rec.is_reverse else "+"
        t3 = s if rec.is_reverse else e          # 转录方向的 3′ 端
        for i in cand:
            if gst[i] != strand or ge[i] < s or gs[i] > e:
                continue
            # 3′ 端应落在基因体内或下游 500 bp 内
            if strand == "+" and not (gs[i] <= t3 <= ge[i] + 500):
                continue
            if strand == "-" and not (gs[i] - 500 <= t3 <= ge[i]):
                continue
            three[gid[i]].append(t3)
            n += 1
    bam.close()
    print(f"归属 reads {n}", flush=True)

    fa = pysam.FastaFile(a.genome)
    gmap = g.set_index("gene_id")

    def upstream(chrom, pos, strand, L=60):
        if strand == "+":
            s0 = max(pos - L, 1)
            seq = fa.fetch(chrom, s0 - 1, pos).upper()
        else:
            seq = fa.fetch(chrom, pos - 1, pos - 1 + L).upper()
            seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        return seq

    site_rows, gene_rows = [], []
    for gene, pos in three.items():
        if len(pos) < MIN_READS:
            continue
        r = gmap.loc[gene]
        v = np.array(sorted(pos))
        cnt = Counter((v // BIN).tolist())
        peaks = [(k, c) for k, c in cnt.items() if c >= max(MIN_READS, MIN_FRAC * len(v))]
        peaks.sort(key=lambda x: -x[1])
        if not peaks:
            continue
        total = len(v)
        info = []
        for k, c in peaks:
            sel = v[(v >= k * BIN) & (v < (k + 1) * BIN)]
            p = int(np.median(sel))
            seq = upstream(r.chrom, p, r.strand)
            hit, off = None, None
            H = len(seq)                       # seq 末端就是 PAS
            best = None
            for sig in SIGNALS:
                for k in range(WIN_LO, WIN_HI + 1):
                    st = H - k
                    if st >= 0 and seq[st:st + 6] == sig:
                        # 越靠近 -20 越优先
                        score = abs(k - 20)
                        if best is None or score < best[0]:
                            best = (score, sig, k)
                        break
            if best:
                hit, off = best[1], best[2]
            info.append(dict(pas=p, n=c, frac=round(c / total, 3), motif=hit, motif_pos=off))
            site_rows.append(dict(gene_id=gene, chrom=r.chrom, pas=p, strand=r.strand,
                                  n_reads=c, frac=round(c / total, 3),
                                  motif=hit, motif_pos=off))
        info.sort(key=lambda x: x["pas"] if r.strand == "+" else -x["pas"])
        major = max(info, key=lambda x: x["n"])
        gene_rows.append(dict(
            gene_id=gene, biotype=r.biotype, strand=r.strand, n_reads=total,
            n_pas=len(info), major_pas=major["pas"], major_frac=major["frac"],
            major_motif=major["motif"],
            pas_span=(max(x["pas"] for x in info) - min(x["pas"] for x in info)),
            apa=("multiple" if len(info) >= 2 else "single"),
        ))

    sites = pd.DataFrame(site_rows)
    genes = pd.DataFrame(gene_rows)
    sites.to_csv(f"{a.outdir}/pas_sites.tsv", sep="\t", index=False)
    genes.to_csv(f"{a.outdir}/apa_per_gene.tsv", sep="\t", index=False)

    pc = genes[genes.biotype == "protein_coding"]
    summary = dict(
        genes_analysed=int(len(genes)),
        genes_pc=int(len(pc)),
        pas_total=int(len(sites)),
        pas_per_gene_median=float(genes.n_pas.median()),
        apa_genes=int((pc.apa == "multiple").sum()),
        apa_pct=round(100 * (pc.apa == "multiple").mean(), 1),
        motif_found_pct=round(100 * sites.motif.notna().mean(), 1),
        canonical_pct=round(100 * sites.motif.isin(["AATAAA", "ATTAAA"]).mean(), 1),
        motif_top=dict(Counter(sites.motif.dropna()).most_common(6)),
        major_frac_median=round(float(pc.major_frac.median()), 3),
        pas_span_median=int(pc.pas_span.median()),
    )
    json.dump(summary, open(f"{a.outdir}/apa_summary.json", "w"), indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
