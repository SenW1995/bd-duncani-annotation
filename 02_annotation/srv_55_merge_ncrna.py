#!/usr/bin/env python
"""srv_55_merge_ncrna.py — 把 tRNAscan-SE 与 barrnap 结果并入注释，生成 v3.10 GTF

来源
  ~/ncrna/trna.out   tRNAscan-SE 2.0.13 默认模式（48 个，含 1 个 Isotype 未定）
  ~/ncrna/rrna.gff   barrnap 1.10.5 --kingdom fun（18S×2、28S×2、5S×3）

编号
  tRNA: BdWA1_00XtNNNNN（按染色体与坐标连续）
  rRNA: BdWA1_00XrNNNNN
每个位点写 gene/transcript/exon 三行，gene_biotype 分别为 tRNA / rRNA。

用法：python srv_55_merge_ncrna.py --gen-dir .../GEN --ncrna-dir ~/ncrna
输出：GEN/V2/v4_structure/babesia_duncani_v3.10.gtf
      GEN/papers/Paper1_ScientificData/tables/ncrna_added_summary.json
"""
from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

import pandas as pd

SRC_T = "tRNAscan-SE"
SRC_R = "barrnap"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen-dir", default=os.environ.get("GEN_ROOT", "."))
    ap.add_argument("--ncrna-dir", default=os.environ.get("NCRNA_DIR", "./ncrna"))
    a = ap.parse_args()
    out_gtf = f"{a.gen_dir}/V2/v4_structure/babesia_duncani_v3.10.gtf"

    trna = pd.read_csv(f"{a.ncrna_dir}/trna.out", sep=r"\s+", skiprows=3,
                       names=["chrom", "num", "start", "end", "aa", "codon", "intron_begin",
                              "intron_end", "score", "isotype", "isotype_score", "note"],
                       dtype=str)
    trna = trna.dropna(subset=["start"])
    trna["start"] = trna.start.astype(int)
    trna["end"] = trna.end.astype(int)

    rrna = []
    for line in open(f"{a.ncrna_dir}/rrna.gff"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        rrna.append(dict(chrom=p[0], start=int(p[3]), end=int(p[4]), strand=p[6],
                         name=attrs.get("Name", "rRNA"), product=attrs.get("product", "rRNA")))
    rr = pd.DataFrame(rrna)

    def prefix(chrom):
        d = "".join(ch for ch in chrom if ch.isdigit())
        return d.zfill(2) if d else chrom[-2:]

    # 染色体顺序（与现有命名一致）
    order = ["scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5",
             "scaffold_6", "scaffold_8"]
    rows = []
    for kind, df in (("tRNA", trna), ("rRNA", rr)):
        df = df.copy()
        df["order"] = df.chrom.map({c: i for i, c in enumerate(order)})
        df = df.sort_values(["order", "start"])
        tag = "t" if kind == "tRNA" else "r"
        counters = defaultdict(int)
        for _, r in df.iterrows():
            pf = prefix(r["chrom"])
            counters[pf] += 1
            gid = f"BdWA1_{pf}{tag}{counters[pf]:05d}"
            if kind == "tRNA":
                prod = ("tRNA-Undet" if r["aa"] == "Undet" else f"tRNA-{r['aa']}")
                strand = "+" if r["start"] <= r["end"] else "-"
                s, e = sorted((r["start"], r["end"]))
                note = "" if r["aa"] != "Undet" else f" isotype_cm {r['isotype']}"
            else:
                prod = r["product"]
                strand = r["strand"]
                s, e = r["start"], r["end"]
                note = ""
            rows.append(dict(kind=kind, chrom=r.chrom, start=s, end=e, strand=strand,
                             gene_id=gid, product=prod, note=note,
                             score=float(r["score"]) if kind == "tRNA" else None))

    added = pd.DataFrame(rows)
    print(added.groupby("kind", dropna=False)["product"].value_counts().to_string())

    with open(out_gtf, "w", encoding="utf-8") as fh:
        for line in open(f"{a.gen_dir}/V2/v4_structure/babesia_duncani_v3.9.gtf"):
            fh.write(line)
        fh.write("#\n# tRNA and rRNA features added in v3.10 (tRNAscan-SE 2.0.13, barrnap 1.10.5)\n#\n")
        for _, r in added.iterrows():
            bt = r["kind"]
            src = SRC_T if r["kind"] == "tRNA" else SRC_R
            common = f'gene_biotype "{bt}"; product "{r["product"]}";'
            if r["note"]:
                common += f' note "{r["note"].strip()}";'
            if r["kind"] == "tRNA":
                common += f' trna_score "{r["score"]:.1f}";'
            for feat in ("gene", "transcript", "exon"):
                tid = f' transcript_id "{r["gene_id"]}.t1";' if feat != "gene" else ""
                gid = f' gene_id "{r["gene_id"]}";'
                fh.write(f"{r['chrom']}\t{src}\t{feat}\t{r['start']}\t{r['end']}\t.\t{r['strand']}\t.\t"
                         f"{gid}{tid}{common}\n")

    summ = {
        "tRNA_added": int((added["kind"] == "tRNA").sum()) if len(added) else 0,
        "tRNA_undet": int((added["product"] == "tRNA-Undet").sum()) if len(added) else 0,
        "rRNA_added": int((added["kind"] == "rRNA").sum()) if len(added) else 0,
        "rRNA_types": added[added["kind"] == "rRNA"]["name" if "name" in added.columns else "product"].value_counts().astype(int).to_dict() if len(added) else {},
        "genes_v3_9": 4039,
        "genes_v3_10": 4039 + len(added),
    }
    os.makedirs(f"{a.gen_dir}/papers/Paper1_ScientificData/tables", exist_ok=True)
    with open(f"{a.gen_dir}/papers/Paper1_ScientificData/tables/ncrna_added_summary.json", "w") as fh:
        json.dump(summ, fh, indent=2)
    added.to_csv(f"{a.gen_dir}/papers/Paper1_ScientificData/tables/ncrna_added.tsv", sep="\t", index=False)
    print(json.dumps(summ, indent=2))
    print("wrote", out_gtf)


if __name__ == "__main__":
    main()
