#!/usr/bin/env python3
"""用全长转录组重建每个基因的 UTR，输出带 five_prime_utr / three_prime_utr 的 GTF。

与本地 v2.1 的区别：v2.1 用的是深度分位数，这里用"剪接链与该基因完全匹配的
全长 reads"的真实 5′/3′ 端众数，边界更准。

用法: srv_05_utr.py --bam WT4.bam --gtf v2.1.gtf --out v3.0.gtf --table utr.tsv
"""
from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam

WIN = 1000
BIN = 25          # 5′/3′ 端聚类窗口
MIN_READS = 3     # 至少几条全长 reads 支持才改边界


def parse_gtf(path):
    genes, tx = {}, defaultdict(lambda: dict(exons=[], cds=[]))
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        gid = a.get("gene_id")
        if not gid:
            continue
        if f[2] == "gene":
            genes[gid] = dict(gene_id=gid, chrom=f[0], start=int(f[3]), end=int(f[4]),
                              strand=f[6], biotype=a.get("gene_biotype", ""), attrs=f[8])
        elif f[2] == "exon":
            tx[gid]["exons"].append((int(f[3]), int(f[4])))
        elif f[2] == "CDS":
            tx[gid]["cds"].append((int(f[3]), int(f[4])))
    return genes, tx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--min-mapq", type=int, default=10)
    a = ap.parse_args()

    genes, tx = parse_gtf(a.gtf)
    g = pd.DataFrame(genes.values()).reset_index(drop=True)
    print(f"基因 {len(g)}", flush=True)

    # 每个基因的内含子集合（用于判断 read 是否属于该基因）
    introns = {}
    for gid, t in tx.items():
        ex = sorted(set(t["exons"]))
        if len(ex) < 2:
            introns[gid] = frozenset()
            continue
        introns[gid] = frozenset((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))

    # 内含子 -> 基因 的哈希表
    jmap = defaultdict(set)
    for gid, ij in introns.items():
        for (js, je) in ij:
            jmap[(genes[gid]["chrom"], js, je)].add(gid)
    jexact = {k: frozenset(v) for k, v in jmap.items()}

    # 窗口索引
    win = defaultdict(list)
    gs, ge = g.start.to_numpy(), g.end.to_numpy()
    gid_arr, gst = g.gene_id.to_numpy(), g.strand.to_numpy()
    for i in range(len(g)):
        for w in range(gs[i] // WIN, ge[i] // WIN + 1):
            win[(g.chrom.iat[i], w)].append(i)

    five = defaultdict(list)
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
        # read 的剪接位点
        rj = set()
        rpos = s
        for op, ln in rec.cigartuples:
            if op == 3:
                rj.add((rpos + 1, rpos + ln))
                rpos += ln
            else:
                rpos += ln
        f5, t3 = (e, s) if rec.is_reverse else (s, e)
        # 先把 read 的剪接位点映射到基因（±5 bp），避免逐基因比对内含子集合
        js_genes = None
        for js, je in rj:
            key = (rec.reference_name, js, je)
            gset = jexact.get(key)
            if gset is None:          # 精确命中失败才做 ±5 展开
                gset = set()
                for dx in range(-5, 6):
                    for dy in range(-5, 6):
                        gset.update(jmap.get((rec.reference_name, js + dx, je + dy), ()))
            js_genes = set(gset) if js_genes is None else (js_genes & gset)
        for i in cand:
            if gst[i] != strand or ge[i] < s or gs[i] > e:
                continue
            gid = gid_arr[i]
            if introns.get(gid) and js_genes is not None and gid not in js_genes:
                continue
            five[gid].append(f5)
            three[gid].append(t3)
            n += 1
    bam.close()
    print(f"归属到基因的 reads: {n}", flush=True)

    def peak(vals, strand, is5):
        """取众数所在的 BIN 区间，返回该区间的中位位置"""
        if not vals:
            return None
        v = np.array(sorted(vals))
        cnt = Counter((v // BIN).tolist())
        k = max(cnt.items(), key=lambda x: (x[1], -abs(x[0])))[0]
        sel = v[(v >= k * BIN) & (v < (k + 1) * BIN)]
        return int(np.median(sel)) if len(sel) else None

    rows = []
    for gid, gg in genes.items():
        f5, t3 = five.get(gid, []), three.get(gid, [])
        rec = dict(gene_id=gid, chrom=gg["chrom"], strand=gg["strand"],
                   biotype=gg["biotype"], gene_start=gg["start"], gene_end=gg["end"],
                   n_support=len(f5), tss=None, tes=None, utr5=0, utr3=0, used=False)
        if len(f5) >= MIN_READS:
            tss = peak(f5, gg["strand"], True)
            tes = peak(t3, gg["strand"], False)
            rec.update(tss=tss, tes=tes)
            if tss is not None:
                rec["utr5"] = (gg["start"] - tss) if gg["strand"] == "+" else (tss - gg["end"])
            if tes is not None:
                rec["utr3"] = (tes - gg["end"]) if gg["strand"] == "+" else (gg["start"] - tes)
            rec["used"] = True
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.to_csv(a.table, sep="\t", index=False)
    print("5′UTR 中位 %d，3′UTR 中位 %d，可用基因 %d"
          % (t.utr5.median(), t.utr3.median(), int(t.used.sum())), flush=True)

    # ---------- 写出 GTF ----------
    out = []
    for gid, gg in genes.items():
        r = t[t.gene_id == gid].iloc[0]
        ex = sorted(set(tx[gid]["exons"]))
        cds = sorted(set(tx[gid]["cds"]))
        if not ex:
            continue
        strand = gg["strand"]
        start, end = gg["start"], gg["end"]
        ex = list(ex)
        # 用全长 reads 的真实端点外延首/末外显子（上限 1500 bp，避免越界到邻居）
        if r.used and r.tss is not None:
            d = r.utr5
            if strand == "+" and 0 < d <= 1500:
                ex[0] = (int(min(ex[0][0], r.tss)), ex[0][1])
                start = int(min(start, r.tss))
            elif strand == "-" and 0 < d <= 1500:
                ex[-1] = (ex[-1][0], int(max(ex[-1][1], r.tss)))
                end = int(max(end, r.tss))
        if r.used and r.tes is not None:
            d = r.utr3
            if strand == "+" and 0 < d <= 1500:
                ex[-1] = (ex[-1][0], int(max(ex[-1][1], r.tes)))
                end = int(max(end, r.tes))
            elif strand == "-" and 0 < d <= 1500:
                ex[0] = (int(min(ex[0][0], r.tes)), ex[0][1])
                start = int(min(start, r.tes))

        attr_g = gg["attrs"]
        if "qc_note" not in attr_g and r.used:
            attr_g = attr_g.rstrip().rstrip(";") + '; utr_source "full_length_reads";'
        out.append((gg["chrom"], start, end, strand, f'gene\t{attr_g}'))

        tid = gid
        tattr = attr_g
        out.append((gg["chrom"], start, end, strand, f'transcript\t{tattr}'))
        for s_, e_ in ex:
            out.append((gg["chrom"], s_, e_, strand, f'exon\tgene_id "{gid}"; transcript_id "{tid}";'))
        for s_, e_ in cds:
            out.append((gg["chrom"], s_, e_, strand, f'CDS\tgene_id "{gid}"; transcript_id "{tid}";'))

        # 显式 UTR feature
        if cds:
            cds_lo, cds_hi = min(c[0] for c in cds), max(c[1] for c in cds)
            if strand == "+":
                if cds_lo > start:
                    out.append((gg["chrom"], start, cds_lo - 1, strand,
                                f'five_prime_utr\tgene_id "{gid}"; transcript_id "{tid}";'))
                if end > cds_hi:
                    out.append((gg["chrom"], cds_hi + 1, end, strand,
                                f'three_prime_utr\tgene_id "{gid}"; transcript_id "{tid}";'))
            else:
                if end > cds_hi:
                    out.append((gg["chrom"], cds_hi + 1, end, strand,
                                f'five_prime_utr\tgene_id "{gid}"; transcript_id "{tid}";'))
                if cds_lo > start:
                    out.append((gg["chrom"], start, cds_lo - 1, strand,
                                f'three_prime_utr\tgene_id "{gid}"; transcript_id "{tid}";'))

    order = {c: i for i, c in enumerate(sorted({x[0] for x in out}))}
    out.sort(key=lambda x: (order[x[0]], x[1], x[2]))
    with open(a.out, "w") as fh:
        fh.write("#!genome-build Babesia_duncani\n#!annotation v3.0 (UTR rebuilt from full-length RNA)\n")
        for chrom, s_, e_, strand, body in out:
            feat, attr = body.split("\t", 1)
            s_, e_ = int(round(float(s_))), int(round(float(e_)))
            fh.write(f"{chrom}\tbd_v3\t{feat}\t{s_}\t{e_}\t.\t{strand}\t.\t{attr}\n")
    print("写出", a.out, len(out), "行", flush=True)


if __name__ == "__main__":
    main()
