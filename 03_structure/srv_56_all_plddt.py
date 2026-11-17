#!/usr/bin/env python3
"""统计全部 Boltz-2 模型的 pLDDT（逐模型 mean over CA atoms）

用途：Paper 1 需要"全部预测蛋白"的 pLDDT 分布，而此前只有用于结构检索的
784 个模型（structure_confidence_all.tsv）和 1,255 个代表模型的汇总
（plddt_all.tsv）。本脚本直接解析 structures/*.pdb 的 B-factor 列（Boltz/
AlphaFold 约定为 pLDDT，0–100），给出全量口径。

用法（WSL bioinfo）
  python3 srv_56_all_plddt.py --pdb-dir <WORK_DIR>/pdb_all --out <TABLE_DIR>
"""
import argparse
import glob
import json
import os
import sys

import numpy as np


def model_plddt(path):
    """返回 (n_ca, mean, frac_ge70, frac_ge90, frac_ge50)"""
    vals = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                continue
    if not vals:
        return None
    a = np.asarray(vals, dtype=float)
    return (int(a.size), float(a.mean()),
            float((a >= 70).mean()), float((a >= 90).mean()), float((a >= 50).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True, help="输出目录")
    ap.add_argument("--pattern", default="*.pdb")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.pdb_dir, a.pattern)))
    if not files:
        sys.exit(f"no PDB under {a.pdb_dir}")

    rows = []
    for i, p in enumerate(files, 1):
        r = model_plddt(p)
        if r is None:
            continue
        gid = os.path.basename(p).replace(".pdb", "")
        rows.append(dict(pdb_id=gid, gene_id=gid, n_res=r[0], mean_pLDDT=round(r[1], 2),
                         frac_ge70=round(r[2], 4), frac_ge90=round(r[3], 4),
                         frac_ge50=round(r[4], 4)))
        if i % 500 == 0:
            print(f"  {i}/{len(files)}", flush=True)

    os.makedirs(a.out, exist_ok=True)
    tsv = os.path.join(a.out, "plddt_all_models.tsv")
    with open(tsv, "w", encoding="utf-8") as fh:
        cols = ["pdb_id", "gene_id", "n_res", "mean_pLDDT", "frac_ge70", "frac_ge90", "frac_ge50"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    m = np.array([r["mean_pLDDT"] for r in rows])
    s = dict(
        n_models=int(m.size),
        mean=round(float(m.mean()), 2),
        median=round(float(np.median(m)), 2),
        q1=round(float(np.percentile(m, 25)), 2),
        q3=round(float(np.percentile(m, 75)), 2),
        min=round(float(m.min()), 2),
        max=round(float(m.max()), 2),
        n_ge70=int((m >= 70).sum()),
        n_50_70=int(((m >= 50) & (m < 70)).sum()),
        n_lt50=int((m < 50).sum()),
        pct_ge70=round(float((m >= 70).mean() * 100), 1),
        pct_50_70=round(float(((m >= 50) & (m < 70)).mean() * 100), 1),
        pct_lt50=round(float((m < 50).mean() * 100), 1),
        n_res_median=int(np.median([r["n_res"] for r in rows])),
    )
    with open(os.path.join(a.out, "plddt_all_models_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(s, fh, indent=2, ensure_ascii=False)

    print(json.dumps(s, indent=2, ensure_ascii=False))
    print(f"wrote {tsv}")


if __name__ == "__main__":
    main()
