#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""featureCounts 输出 → counts / TPM / logCPM 表。

用法
  python counts_to_tpm.py --counts counts.txt --out gene_expression.tsv
"""
import argparse
import math


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", default="culture,invivo")
    a = ap.parse_args()
    lab = a.labels.split(",")

    rows = []
    for line in open(a.counts, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "Geneid":
            continue
        g, length = f[0], int(float(f[5]))
        c = [int(float(x)) for x in f[6:6 + len(lab)]]
        rows.append((g, length, c))

    lib = [sum(r[2][i] for r in rows) for i in range(len(lab))]
    with open(a.out, "w", encoding="utf-8") as fh:
        cols = ["gene_id", "length"]
        cols += [f"count_{x}" for x in lab]
        cols += [f"tpm_{x}" for x in lab]
        cols += [f"logCPM_{x}" for x in lab]
        fh.write("\t".join(cols) + "\n")
        for g, length, c in rows:
            tpm = [c[i] / (length / 1000.0) / lib[i] * 1e6 if lib[i] else 0.0 for i in range(len(lab))]
            lc = [math.log2(c[i] / lib[i] * 1e6 + 1) if lib[i] else 0.0 for i in range(len(lab))]
            fh.write("\t".join([g, str(length)] + [str(x) for x in c]
                               + [f"{x:.2f}" for x in tpm] + [f"{x:.3f}" for x in lc]) + "\n")
    print("genes:", len(rows), "library sizes:", lib)


if __name__ == "__main__":
    main()
