#!/usr/bin/env python3
"""用全长转录组（Nanopore cDNA）验证基因组注释。

一次 pysam 遍历完成：基因/外显子 reads 计数、剪接位点、全长比例；
再用 samtools depth 做覆盖度与未注释转录区（内存占用低，避免 OOM）。

用法: fl_validate.py --bam x.mq.bam --gtf anno.gtf --genome genome.fa --outdir out/

输出
  01_gene_support.tsv   每个基因的 reads 支持、外显子命中率、外显子平均深度
  02_introns.tsv        注释内含子及 reads 支持
  03_junctions.tsv      观测到的剪接位点（含是否注释）
  04_novel_loci.tsv     有转录本覆盖但不在任何基因体内的区段
  summary.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import pysam

WIN = 1000  # 分窗索引宽度


def parse_gtf(path: str):
    genes, exons = {}, []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] not in ("gene", "exon"):
            continue
        attr = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        gid = attr.get("gene_id")
        if not gid:
            continue
        if f[2] == "gene":
            genes[gid] = dict(gene_id=gid, chrom=f[0], start=int(f[3]), end=int(f[4]),
                              strand=f[6], biotype=attr.get("gene_biotype", ""),
                              product=attr.get("product", ""))
        else:
            exons.append((gid, f[0], int(f[3]), int(f[4])))
    ex = pd.DataFrame(exons, columns=["gene_id", "chrom", "start", "end"]) \
        .drop_duplicates(["gene_id", "start", "end"]).reset_index(drop=True)
    ex["ex_id"] = ex.gene_id + ":" + ex.start.astype(str) + "-" + ex.end.astype(str)
    ex["len"] = ex.end - ex.start + 1
    return genes, ex


class WinIndex:
    """按 1 kb 窗口建索引，查询时只检查落在同一窗口的区间。"""

    def __init__(self, df: pd.DataFrame):
        self.rows = df.reset_index(drop=True)
        self.idx = defaultdict(list)
        starts = self.rows.start.to_numpy()
        ends = self.rows.end.to_numpy()
        chroms = self.rows.chrom.to_numpy()
        for i in range(len(self.rows)):
            for w in range(starts[i] // WIN, ends[i] // WIN + 1):
                self.idx[(chroms[i], w)].append(i)
        self.starts, self.ends = starts, ends

    def query(self, chrom: str, s: int, e: int) -> np.ndarray:
        cand = []
        for w in range(s // WIN, e // WIN + 1):
            cand.extend(self.idx.get((chrom, w), ()))
        if not cand:
            return np.zeros(0, dtype=int)
        cand = np.unique(np.asarray(cand, dtype=int))
        ov = np.minimum(self.ends[cand], e) - np.maximum(self.starts[cand], s) + 1
        return cand[ov > 0]


def sh(cmd: list[str], stdout_path: str):
    with open(stdout_path, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-junc-reads", type=int, default=3)
    ap.add_argument("--min-novel-cov", type=int, default=50)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    O = args.outdir

    genes, ex = parse_gtf(args.gtf)
    gdf = pd.DataFrame(genes.values())
    print(f"genes: {len(gdf)}, exons: {len(ex)}", flush=True)

    gidx = WinIndex(gdf)
    eidx = WinIndex(ex)

    # ---------- 单次遍历 ----------
    gene_cnt, exon_cnt, junc = Counter(), Counter(), Counter()
    jstrand = {}
    n_read = n_primary = n_fl = 0
    qfrac_sum = 0.0
    cache = f"{O}/00_junctions_raw.tsv"
    bam = pysam.AlignmentFile(args.bam, "rb")
    for rec in ([] if os.path.exists(cache) else bam.fetch(until_eof=True)):
        n_read += 1
        if rec.is_secondary or rec.is_supplementary or rec.mapping_quality < args.min_mapq:
            continue
        n_primary += 1
        if rec.query_length:
            qf = rec.query_alignment_length / rec.query_length
            qfrac_sum += qf
            if qf >= 0.9:
                n_fl += 1
        chrom = rec.reference_name
        s = rec.reference_start + 1
        e = rec.reference_end
        for gi in gidx.query(chrom, s, e):
            gene_cnt[gidx.rows.gene_id.iat[gi]] += 1
        for ei in eidx.query(chrom, s, e):
            exon_cnt[eidx.rows.ex_id.iat[ei]] += 1
        rpos = s
        for op, ln in rec.cigartuples:
            if op == 3:
                key = (chrom, rpos + 1, rpos + ln)
                junc[key] += 1
                jstrand[key] = "-" if rec.is_reverse else "+"
                rpos += ln
            else:
                rpos += ln
    bam.close()
    if not os.path.exists(cache):
        pd.DataFrame([dict(chrom=k[0], start=k[1], end=k[2], strand=jstrand[k],
                           n_reads=v) for k, v in junc.items()]) \
            .to_csv(cache, sep="\t", index=False)
    else:
        _c = pd.read_csv(cache, sep="\t")
        junc = Counter({(r.chrom, int(r.start), int(r.end)): int(r.n_reads) for _, r in _c.iterrows()})
        jstrand = {(r.chrom, int(r.start), int(r.end)): r.strand for _, r in _c.iterrows()}
    print(f"reads: {n_read}, primary: {n_primary}, junctions: {len(jstrand)}", flush=True)

    # ---------- 每个基因的支持度 ----------
    gs = gdf.copy()
    gs["n_reads"] = gs.gene_id.map(gene_cnt).fillna(0).astype(int)
    exn = ex.assign(n_reads=ex.ex_id.map(exon_cnt).fillna(0).astype(int))
    ex_stat = exn.groupby("gene_id").agg(
        n_exons=("ex_id", "size"),
        n_exons_hit=("n_reads", lambda s: int((s > 0).sum()))).reset_index()
    ex_stat["frac_exons_hit"] = (ex_stat.n_exons_hit / ex_stat.n_exons).round(3)
    gs = gs.merge(ex_stat, on="gene_id", how="left")
    gs["n_exons"] = gs.n_exons.fillna(0).astype(int)
    gs["frac_exons_hit"] = gs.frac_exons_hit.fillna(0)
    gs["span"] = gs.end - gs.start + 1

    # ---------- 深度 ----------
    dep_path = f"{O}/tmp_depth.tsv"
    if not (os.path.exists(dep_path) and os.path.getsize(dep_path) > 0):
        sh(["samtools", "depth", "-a", "-Q", str(args.min_mapq), "-q", "0", args.bam], dep_path)
    dep = pd.read_csv(dep_path, sep="\t", header=None, names=["chrom", "pos", "depth"])
    dep_by_chrom = {c: g.depth.to_numpy(dtype=np.int32) for c, g in dep.groupby("chrom", sort=False)}
    fa = pysam.FastaFile(args.genome)
    glen = {c: fa.get_reference_length(c) for c in fa.references}
    fa.close()

    ex_depth = []
    for _, r in exn.iterrows():
        arr = dep_by_chrom.get(r.chrom)
        if arr is None:
            ex_depth.append(0.0)
            continue
        seg = arr[max(r.start - 1, 0):min(r.end, len(arr))]
        ex_depth.append(float(seg.mean()) if len(seg) else 0.0)
    exn["mean_depth"] = ex_depth
    agg = exn.groupby("gene_id").apply(
        lambda d: pd.Series(dict(
            mean_exon_depth=float((d.len * d.mean_depth).sum() / d.len.sum()),
            frac_exon_bp_10x=float((d.mean_depth >= 10).mul(d.len).sum() / d.len.sum()))),
        include_groups=False).reset_index()
    gs = gs.merge(agg, on="gene_id", how="left")
    gs = gs.sort_values(["n_reads", "frac_exons_hit"], ascending=False)
    gs.to_csv(f"{O}/01_gene_support.tsv", sep="\t", index=False)

    # ---------- 剪接位点 ----------
    rows = []
    for gid, g in ex.groupby("gene_id"):
        g = g.sort_values("start")
        if len(g) < 2:
            continue
        s, e = g.start.to_numpy(), g.end.to_numpy()
        chrom, strand = g.chrom.iloc[0], genes[gid]["strand"]
        for i in range(len(s) - 1):
            if s[i + 1] - e[i] - 1 >= 20:
                rows.append(dict(gene_id=gid, chrom=chrom, strand=strand,
                                 start=int(e[i]) + 1, end=int(s[i + 1]) - 1))
    intr = pd.DataFrame(rows)

    jdf = pd.DataFrame([dict(chrom=k[0], start=k[1], end=k[2], strand=jstrand[k], n_reads=v)
                        for k, v in junc.items() if v >= args.min_junc_reads])
    jset = {(r.chrom, r.start, r.end) for _, r in jdf.iterrows()} if len(jdf) else set()

    def supported(c, s, e):
        for ds in range(-5, 6):
            for de in range(-5, 6):
                if (c, s + ds, e + de) in jset:
                    return True
        return False

    if not intr.empty:
        intr["supported"] = [supported(r.chrom, r.start, r.end) for _, r in intr.iterrows()]
        intr = intr.merge(gs[["gene_id", "n_reads"]].rename(columns={"n_reads": "gene_reads"}),
                          on="gene_id", how="left")
        intr.to_csv(f"{O}/02_introns.tsv", sep="\t", index=False)

    if not jdf.empty:
        ann_exact = {(r.chrom, r.start, r.end): r.gene_id for _, r in intr.iterrows()}
        gtag = []
        for _, r in jdf.iterrows():
            hit = ""
            if (r.chrom, r.start, r.end) in ann_exact:
                hit = ann_exact[(r.chrom, r.start, r.end)]
            else:
                for ds in range(-5, 6):
                    for de in range(-5, 6):
                        g = ann_exact.get((r.chrom, r.start + ds, r.end + de))
                        if g:
                            hit = g
                            break
                    if hit:
                        break
            gtag.append(hit)
        jdf["gene_id"] = gtag
        jdf["annotated"] = jdf.gene_id != ""
        jdf = jdf.sort_values("n_reads", ascending=False)
    jdf.to_csv(f"{O}/03_junctions.tsv", sep="\t", index=False)

    # ---------- 未注释转录区 ----------
    novel = []
    for c, arr in dep_by_chrom.items():
        if c not in glen:
            continue
        mask = np.zeros(glen[c], dtype=bool)
        n = min(len(arr), glen[c])
        mask[:n] = arr[:n] >= args.min_novel_cov
        sub = gs[gs.chrom == c]
        for _, r in sub.iterrows():
            mask[max(r.start - 1, 0):min(r.end, len(mask))] = False
        d = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
        starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
        for s, e in zip(starts, ends):
            if e - s >= 300:
                novel.append(dict(chrom=c, start=int(s) + 1, end=int(e),
                                  span=int(e - s), mean_depth=round(float(arr[s:e].mean()), 1)))
    nv = pd.DataFrame(novel).sort_values("span", ascending=False) if novel else \
        pd.DataFrame(columns=["chrom", "start", "end", "span", "mean_depth"])
    nv.to_csv(f"{O}/04_novel_loci.tsv", sep="\t", index=False)

    # ---------- 汇总 ----------
    pc = gs[gs.biotype == "protein_coding"]
    summary = dict(
        genes_total=int(len(gs)), genes_pc=int(len(pc)),
        reads_scanned=int(n_read), reads_primary=int(n_primary),
        reads_query_aligned_90pct=int(n_fl),
        reads_query_aligned_90pct_pct=round(100 * n_fl / max(n_primary, 1), 1),
        mean_query_aligned_frac=round(qfrac_sum / max(n_primary, 1), 3),
        pc_with_reads=int((pc.n_reads > 0).sum()),
        pc_with_reads_pct=round(100 * (pc.n_reads > 0).mean(), 1),
        pc_well_supported=int(((pc.n_reads >= 10) & (pc.frac_exons_hit >= 0.9)).sum()),
        pc_well_supported_pct=round(100 * ((pc.n_reads >= 10) & (pc.frac_exons_hit >= 0.9)).mean(), 1),
        pc_no_reads=int((pc.n_reads == 0).sum()),
        introns_total=int(len(intr)),
        introns_supported=int(intr.supported.sum()) if len(intr) else 0,
        introns_supported_pct=round(100 * intr.supported.mean(), 1) if len(intr) else None,
        junctions_obs=int(len(jdf)),
        junctions_annotated=int(jdf.annotated.sum()) if len(jdf) else 0,
        junctions_novel=int((~jdf.annotated).sum()) if len(jdf) else 0,
        novel_loci=int(len(nv)), novel_loci_bp=int(nv.span.sum()) if len(nv) else 0,
    )
    json.dump(summary, open(f"{O}/summary.json", "w"), indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if os.path.exists(dep_path):
        os.remove(dep_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
