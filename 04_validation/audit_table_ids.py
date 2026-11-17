#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核对提交表格中的基因编号，输出可映射 / 已删除 / 无法解析三类清单。

用法
  python audit_table_ids.py --gtf <gtf> --map <v38_id_map.tsv> --tables <tsv...> --outdir <dir>
"""
import argparse
import collections
import json
import os
import re

GID = re.compile(r'^BdWA1_(?:nc[0-9]{3}[0-9]{5}|[0-9]{3}(?:pg|g|t|r)[0-9]{5})$')
IDCOL = ("gene_id", "pdb_id", "lncRNA", "lncrna_id", "partner_gene")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--tables", nargs="*", default=[])
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    gtf = {}
    for line in open(a.gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        at = dict(re.findall(r'(\w+)\s+"([^"]*)"', f[8]))
        g = at.get("gene_id")
        gtf[g] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                      biotype=at.get("gene_biotype", ""))

    old2new, removed = {}, set()
    for line in open(a.map, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or f[0] == "old_id":
            continue
        if f[1]:
            old2new[f[0]] = f[1]
        else:
            removed.add(f[0])

    summary, detail = [], []
    for tf in a.tables:
        if not os.path.exists(tf):
            continue
        rows = [l.rstrip("\n").split("\t") for l in open(tf, encoding="utf-8", errors="replace") if l.strip()]
        hdr = rows[0]
        cols = []
        for key in IDCOL:
            if key in hdr:
                cols.append(hdr.index(key))
        if not cols:
            continue
        seen, stat = set(), collections.Counter()
        for r in rows[1:]:
            for ci in cols:
                if ci >= len(r):
                    continue
                g = r[ci]
                if not GID.match(g) or g in gtf or g in seen:
                    continue
                seen.add(g)
                if g in old2new:
                    stat["remap"] += 1
                    detail.append([os.path.basename(tf), g, "remap", old2new[g], ""])
                elif g in removed:
                    stat["removed"] += 1
                    detail.append([os.path.basename(tf), g, "removed", "", ""])
                else:
                    # 按坐标在同 scaffold 找重叠基因（用于判断是否早期版本号）
                    stat["unknown"] += 1
                    detail.append([os.path.basename(tf), g, "unknown", "", ""])
        summary.append(dict(table=os.path.basename(tf), n_bad=len(seen), **stat))

    with open(os.path.join(a.outdir, "table_id_detail.tsv"), "w", encoding="utf-8") as fh:
        fh.write("table\told_id\tstatus\tnew_id\tnote\n")
        for r in detail:
            fh.write("\t".join(r) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
