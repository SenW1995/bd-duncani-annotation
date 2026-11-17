#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 5 个读框断裂的 coding 模型转成假基因（v3.12 → v3.13）。

这 5 个模型 CDS 长度都不是 3 的倍数，翻译出 10–95 个内部终止，table2asn 只报
SuspiciousFrame，因此没有进入此前 51 个 InternalStop 的修复名单。

处理：
  1. gene / transcript 的 biotype 改为 pseudogene
  2. 加 qc_note "frameshift_unresolved_pseudogene"
  3. 写 pseudogene_type（由 BLASTP 最长 ORF 决定：有旁系同源 unprocessed，否则 unitary）
CDS 坐标保留，与已转的 42 个一致。

用法
  python srv_73_fix_frameshift.py -i v3.12.gtf --types frameshift_types.tsv -o v3.13.gtf
"""
import argparse
import collections
import re

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
TARGETS = {"BdWA1_001g00004", "BdWA1_001g00007", "BdWA1_001g00020",
           "BdWA1_001g00024", "BdWA1_002g00092"}


def set_attr(s, key, val):
    if re.search(key + r'\s+"[^"]*"', s):
        return re.sub(key + r'\s+"[^"]*"', f'{key} "{val}"', s, count=1)
    return s.rstrip(";") + f'; {key} "{val}";'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("--types", required=True, help="gene<TAB>type")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()

    types = {}
    for line in open(a.types, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0].startswith("BdWA1"):
            types[f[0]] = f[1]

    n = collections.Counter()
    with open(a.output, "w", encoding="utf-8") as out:
        for line in open(a.input, encoding="utf-8"):
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                out.write(line)
                continue
            at = dict(ATTR_RE.findall(f[8]))
            g = at.get("gene_id")
            if g not in TARGETS:
                out.write(line)
                continue
            s = f[8]
            if f[2] in ("gene", "transcript"):
                # 与既有假基因写法一致：gene 与 transcript 都用 gene_biotype
                s = set_attr(s, "gene_biotype", "pseudogene")
                s = re.sub(r'transcript_biotype\s+"[^"]*";\s*', "", s)
                s = set_attr(s, "pseudogene_type", types.get(g, "unitary"))
                s = set_attr(s, "qc_note", "frameshift_unresolved_pseudogene")
                n[f[2]] += 1
            f[8] = s
            out.write("\t".join(f) + "\n")
    print("改写特征数:", dict(n))
    print("写出:", a.output)


if __name__ == "__main__":
    main()
