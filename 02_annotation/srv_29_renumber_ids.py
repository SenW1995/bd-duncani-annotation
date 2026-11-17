#!/usr/bin/env python3
"""把非标准基因号改成统一的 BdWA1_00XgNNNNN 连续编号（v3.5 定稿）

规则（沿用原注释）
  BdWA1_{scaffold 三位}g{该 scaffold 内五位流水号}
  新号从该 scaffold 现有最大号往后接，按文件顺序（即坐标顺序）分配。
  非标准号包括：lncRNA 的 BdWA1_nc00100001 系列、内含子来源基因的 ..._hid1 系列。

同时输出 old->new 对照表，供下游表（家族、lncRNA 靶基因）改名用。
"""
from __future__ import annotations

import re

IN = "<GEN_ROOT>/V2/v3/babesia_duncani_v3.5_pre.gtf"
OUT = "<GEN_ROOT>/V2/v3/babesia_duncani_v3.5.gtf"
MAP = "<GEN_ROOT>/V2/v3/v35_id_map.tsv"

STD = re.compile(r"^BdWA1_(\d{3})g(\d{5})$")


def main():
    lines = open(IN, encoding="utf-8").read().splitlines(keepends=True)

    # 1) 每个 scaffold 现有最大流水号
    maxidx = {}
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        m = STD.match(a.get("gene_id", ""))
        if m:
            scaf, idx = int(m.group(1)), int(m.group(2))
            maxidx[scaf] = max(maxidx.get(scaf, 0), idx)

    # 2) 非标准号 -> 新号（按出现顺序）
    remap, order = {}, []
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        gid = a.get("gene_id", "")
        if not gid or STD.match(gid) or gid in remap:
            continue
        scaf = int(f[0].split("_")[1])
        maxidx[scaf] = maxidx.get(scaf, 0) + 1
        new = f"BdWA1_{scaf:03d}g{maxidx[scaf]:05d}"
        remap[gid] = new
        order.append((gid, new, f[0], int(f[3]), int(f[4]), a.get("gene_biotype", "")))

    # 3) 全文件替换（长的先替换，避免 _hid1 被父基因号截断）
    keys = sorted(remap, key=len, reverse=True)
    out = []
    for ln in lines:
        for k in keys:
            if k in ln:
                ln = re.sub(re.escape(k) + r"(?![0-9A-Za-z_])", remap[k], ln)
        out.append(ln)
    open(OUT, "w", encoding="utf-8", newline="\n").write("".join(out))

    with open(MAP, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("old_id\tnew_id\tscaffold\tstart\tend\tbiotype\n")
        for r in sorted(order, key=lambda x: x[1]):
            fh.write("\t".join(map(str, r)) + "\n")

    print(f"非标准号 {len(remap)} 个已重编号 -> {OUT}")
    for gid, new, *_ in order[:3]:
        print(f"  {gid} -> {new}")


if __name__ == "__main__":
    main()
