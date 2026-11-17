#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""汇总最终注释与派生表的关键数字，供论文与表格核对。

用法
  python summarize_resource_numbers.py --tables <dir> --out <json>
其中 <dir> 需包含 01_gene_support.tsv、tss/tss_per_gene.tsv、polya/apa_per_gene.tsv、families/families.tsv
（也接受扁平的旧目录结构）。
"""
import argparse
import collections
import json
import os
import statistics as st


def find(root, name):
    p = os.path.join(root, name)
    if os.path.exists(p):
        return p
    for sub in ("tss", "polya", "families", "."):
        q = os.path.join(root, sub, name)
        if os.path.exists(q):
            return q
    return None


def read(p):
    rows = [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8", errors="replace") if l.strip()]
    h = {k: i for i, k in enumerate(rows[0])}
    return h, rows[1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = {}

    p = find(a.tables, "01_gene_support.tsv")
    if p:
        h, rows = read(p)
        pc = [r for r in rows if r[h["biotype"]] == "protein_coding"]
        R["support"] = dict(
            genes_rows=len(rows), pc_rows=len(pc),
            pc_with_reads=sum(1 for r in pc if int(r[h["n_reads"]]) > 0),
            pc_no_reads=sum(1 for r in pc if int(r[h["n_reads"]]) == 0),
            pc_all_exons=sum(1 for r in pc if float(r[h["frac_exons_hit"]]) == 1.0),
            pc_all_exons_10x=sum(1 for r in pc if float(r[h["frac_exons_hit"]]) == 1.0
                                 and float(r[h["mean_exon_depth"]]) >= 10),
            pc_depth10_frac=round(100 * sum(1 for r in pc if float(r[h["mean_exon_depth"]]) >= 10) / len(pc), 1),
        )

    p = find(a.tables, "tss_per_gene.tsv")
    if p:
        h, rows = read(p)
        pc = [r for r in rows if r[h["biotype"]] == "protein_coding"]
        cls = collections.Counter(r[h["class"]] for r in pc)
        R["tss"] = dict(
            pc_genes=len(pc),
            classes=dict(cls),
            two_or_more=sum(1 for r in pc if int(r[h["n_tss"]]) >= 2),
            max_tss=max(int(r[h["n_tss"]]) for r in pc),
            median_primary_frac=round(st.median(float(r[h["tss1_frac"]]) for r in pc
                                                if r[h["tss1_frac"]]), 3),
            with_upstream_gene=sum(1 for r in pc if r[h["upstream_gene"]]),
        )

    p = find(a.tables, "apa_per_gene.tsv")
    if p:
        h, rows = read(p)
        pc = [r for r in rows if r[h["biotype"]] == "protein_coding"]
        R["pas"] = dict(
            pc_genes=len(pc),
            with_signal=sum(1 for r in pc if r[h["major_motif"]]),
            multiple=sum(1 for r in pc if int(r[h["n_pas"]]) >= 2),
            max_pas=max(int(r[h["n_pas"]]) for r in pc),
            median_major_frac=round(st.median(float(r[h["major_frac"]]) for r in pc
                                              if r[h["major_frac"]]), 3),
            motif_top=dict(collections.Counter(r[h["major_motif"]] for r in pc
                                               if r[h["major_motif"]]).most_common(6)),
        )

    p = find(a.tables, "families.tsv")
    if p:
        h, rows = read(p)
        fam = collections.Counter(r[h["family_id"]] for r in rows)
        R["families"] = dict(
            genes=len(rows), families=len(fam),
            multi=sum(1 for v in fam.values() if v >= 2),
            multi_members=sum(v for v in fam.values() if v >= 2),
            largest=max(fam.values()),
        )

    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(R, fh, ensure_ascii=False, indent=2)
    print(json.dumps(R, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
