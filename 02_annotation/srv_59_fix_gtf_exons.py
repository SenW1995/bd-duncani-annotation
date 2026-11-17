#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重建 GTF 的 exon 记录：exon = merge(CDS ∪ 5'UTR ∪ 3'UTR)。

背景：v3.9/v3.10 GTF 里 404 个蛋白编码基因的 exon 记录是"基因跨度"而非真实外显子，
内含子被算进了 exon（个别基因 exon 段数还少于 CDS 段数）。CDS 与 UTR 记录本身正确，
因此用它们的并集重建 exon 即可。

非编码基因（lncRNA / tRNA / rRNA）没有 CDS，保留原有 exon。

用法
  python srv_59_fix_gtf_exons.py -i in.gtf -o out.gtf
"""
import argparse
import collections
import re
import sys

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_attrs(s):
    return dict(ATTR_RE.findall(s))


def merge_intervals(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:          # 相邻或重叠即合并
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()

    lines = []                    # (原始行, 解析后的字段)
    genes = collections.OrderedDict()   # gene_id -> dict
    for line in open(a.input, encoding="utf-8"):
        if not line.strip() or line.startswith("#"):
            lines.append((line, None))
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            lines.append((line, None))
            continue
        at = parse_attrs(f[8])
        gid = at.get("gene_id")
        feat = f[2]
        rec = dict(seqid=f[0], source=f[1], feat=feat, start=int(f[3]), end=int(f[4]),
                   score=f[5], strand=f[6], frame=f[7], attrs=f[8], raw=line, at=at)
        lines.append((line, rec))
        if gid is None:
            continue
        g = genes.setdefault(gid, dict(biotype=at.get("gene_biotype", ""),
                                       CDS=[], UTR=[], exon=[], order=len(genes)))
        if feat == "CDS":
            g["CDS"].append((rec["start"], rec["end"]))
        elif feat in ("five_prime_utr", "three_prime_utr", "UTR", "five_prime_UTR", "three_prime_UTR"):
            g["UTR"].append((rec["start"], rec["end"]))
        elif feat == "exon":
            g["exon"].append((rec["start"], rec["end"]))

    # 计算每个基因应有的 exon
    new_exon = {}
    n_fixed = 0
    for gid, g in genes.items():
        if not g["CDS"] and not g["UTR"]:
            continue                       # 非编码：保留原 exon
        want = merge_intervals(g["CDS"] + g["UTR"])
        have = sorted(g["exon"])
        if want != have:
            n_fixed += 1
        new_exon[gid] = want

    # 写新 GTF：exon 行按新结构重排；其余行原样
    out = []
    emitted = set()
    # 先按基因归组，保证 gene → transcript → exon/其他 的顺序
    by_gene = collections.OrderedDict()
    for line, rec in lines:
        if rec is None:
            out.append(line)
            continue
        gid = rec["at"].get("gene_id")
        by_gene.setdefault(gid, []).append((line, rec))

    n_added = n_removed = 0
    for gid, recs in by_gene.items():
        if gid in new_exon:
            # 原 exon 全部丢弃，按新结构重写
            kept = [(l, r) for l, r in recs if r["feat"] != "exon"]
            n_removed += len(recs) - len(kept)
            for l, r in kept:
                out.append(l)
            if gid in emitted:
                continue
            emitted.add(gid)
            tmpl = next(r for l, r in recs if r["feat"] == "exon")
            for st, en in new_exon[gid]:
                f = [tmpl["seqid"], tmpl["source"], "exon", str(st), str(en),
                     ".", tmpl["strand"], ".", tmpl["attrs"]]
                out.append("\t".join(f) + "\n")
                n_added += 1
        else:
            for l, r in recs:
                out.append(l)

    with open(a.output, "w", encoding="utf-8") as fh:
        fh.writelines(out)

    print("基因数:", len(genes))
    print("exon 结构被修正的基因:", n_fixed)
    print("原 exon 行删除:", n_removed, "| 新 exon 行写入:", n_added)
    # 统计
    tot_exon = sum(len(v) for v in new_exon.values())
    print("修正后（蛋白编码+含 UTR 基因）exon 总数:", tot_exon)


if __name__ == "__main__":
    main()
