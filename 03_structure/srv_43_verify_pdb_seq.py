#!/usr/bin/env python3
"""核验每个结构 PDB 与其注释蛋白序列是否一致（防止 ID 错配）

做法
  1. 从 PDB 的 CA 原子按残基号取序列（3 字母转 1 字母）
  2. 与 families_v32/all.faa 中同名（旧 ID）序列逐位比较
  3. 输出：完全一致 / 长度不同 / 有替换 / 找不到序列
目的：后续基于结构改注释前，先确认"这个结构确实属于这个基因"。
"""
import os
import sys
import csv

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O"}


def read_faa(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def pdb_seq(path):
    out, seen = [], set()
    for line in open(path):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        resi = line[22:26].strip()
        if resi in seen:
            continue
        seen.add(resi)
        out.append(AA3.get(line[17:20].strip().upper(), "X"))
    return "".join(out)


def main():
    faa, dirs, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
    seqs = read_faa(faa)
    rows = []
    for d in dirs:
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".pdb"):
                continue
            gid = fn[:-4]
            p = os.path.join(d, fn)
            sq = pdb_seq(p)
            ref = seqs.get(gid)
            if ref is None:
                rows.append(dict(pdb_id=gid, status="no_ref_sequence", pdb_len=len(sq),
                                 ref_len=0, identity=""))
                continue
            if sq == ref:
                status = "exact"
            elif len(sq) == len(ref):
                diff = sum(1 for a, b in zip(sq, ref) if a != b)
                status = f"substitutions:{diff}"
            else:
                status = "length_mismatch"
            same = sum(1 for a, b in zip(sq, ref) if a == b)
            ident = same / max(1, len(ref)) * 100
            rows.append(dict(pdb_id=gid, status=status, pdb_len=len(sq),
                             ref_len=len(ref), identity=f"{ident:.1f}"))
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb_id", "status", "pdb_len", "ref_len", "identity"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    import collections
    c = collections.Counter(r["status"].split(":")[0] for r in rows)
    print("核验", len(rows), "个结构:", dict(c))
    bad = [r for r in rows if r["status"] not in ("exact",)]
    for r in bad[:10]:
        print("  异常:", r)


if __name__ == "__main__":
    main()
