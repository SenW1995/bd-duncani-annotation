#!/usr/bin/env python3
"""验证修订版注释：重新提取 CDS、翻译、统计内部终止与结构问题。"""
from __future__ import annotations

import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq

gtf = sys.argv[1]
gen = sys.argv[2]
out = sys.argv[3]

fa = pysam.FastaFile(gen)
cds = defaultdict(list)
bio, chrom_, strand_ = {}, {}, {}
for line in open(gtf):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    gid = a.get("gene_id")
    if f[2] == "CDS":
        cds[gid].append((int(f[3]), int(f[4])))
    elif f[2] == "gene":
        bio[gid] = a.get("gene_biotype", "")
        chrom_[gid] = f[0]
        strand_[gid] = f[6]

rows = []
for gid, segs in cds.items():
    segs = sorted(segs, key=lambda x: x[0])
    nt = "".join(fa.fetch(chrom_[gid], s - 1, e).upper() for s, e in segs)
    if strand_.get(gid) == "-":
        nt = str(Seq(nt).reverse_complement())
    nt = nt[: len(nt) // 3 * 3]
    aa = str(Seq(nt).translate())
    core = aa[:-1] if aa.endswith("*") else aa
    rows.append(dict(gene_id=gid, biotype=bio.get(gid, ""), nt_len=len(nt), aa=len(core),
                     internal_stops=core.count("*"), has_stop=aa.endswith("*")))
d = pd.DataFrame(rows)
d.to_csv(out, sep="\t", index=False)

pc = d[d.biotype == "protein_coding"]
print(f"基因数 {len(d)}  蛋白编码 {len(pc)}")
print(f"蛋白编码中带内部终止: {int((pc.internal_stops > 0).sum())}")
print(f"蛋白编码长度中位数: {int(pc.aa.median())} aa")
print(f"有终止密码子结尾的蛋白编码基因: {int(pc.has_stop.sum())} ({100*pc.has_stop.mean():.1f}%)")
print()
print("biotype 分布:")
print(d.biotype.value_counts().to_string())
