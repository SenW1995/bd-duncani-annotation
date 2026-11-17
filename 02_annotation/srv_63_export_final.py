#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从最终 GTF + 基因组生成 CDS/蛋白 FASTA，并汇总统计。

- CDS：所有有 CDS 特征的基因（含假基因），header 里标 biotype
- 蛋白：仅 protein_coding，翻译后去掉末尾终止符；假基因不翻译
- 统计：基因/biotype/外显子/内含子/长度中位

用法
  python srv_63_export_final.py --gtf v3.12.gtf --genome genome.fa --outdir out
"""
import argparse
import collections
import json
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="babesia_duncani_WA1_v3.12")
    a = ap.parse_args()

    genome = {r.id.split()[0]: str(r.seq).upper() for r in SeqIO.parse(a.genome, "fasta")}

    genes = collections.OrderedDict()
    for line in open(a.gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(ATTR_RE.findall(f[8]))
        g = at.get("gene_id")
        if not g:
            continue
        d = genes.setdefault(g, dict(chrom=f[0], strand=f[6], biotype="", CDS=[], exon=[],
                                     product=""))
        if f[2] == "gene":
            d["biotype"] = at.get("gene_biotype", "")
            d["product"] = at.get("product", "")
        elif f[2] == "CDS":
            d["CDS"].append((int(f[3]), int(f[4])))
        elif f[2] == "exon":
            d["exon"].append((int(f[3]), int(f[4])))

    os.makedirs(a.outdir, exist_ok=True)
    cds_path = os.path.join(a.outdir, f"{a.prefix}.cds.fna")
    prot_path = os.path.join(a.outdir, f"{a.prefix}.protein.faa")
    tx_path = os.path.join(a.outdir, f"{a.prefix}.transcript.fa")

    n_cds = n_prot = n_tx = 0
    with open(cds_path, "w", encoding="utf-8") as fc, \
            open(prot_path, "w", encoding="utf-8") as fp, \
            open(tx_path, "w", encoding="utf-8") as ft:
        for g, d in genes.items():
            if d["exon"]:
                esegs = sorted(d["exon"], reverse=(d["strand"] == "-"))
                eparts = [genome[d["chrom"]][s - 1:e] for s, e in esegs]
                if d["strand"] == "-":
                    eparts = [str(Seq(p).reverse_complement()) for p in eparts]
                tx = "".join(eparts)
                ft.write(f">{g} biotype={d['biotype']} len={len(tx)} product={d['product']}\n")
                for i in range(0, len(tx), 70):
                    ft.write(tx[i:i + 70] + "\n")
                n_tx += 1
            if not d["CDS"]:
                continue
            segs = sorted(d["CDS"], reverse=(d["strand"] == "-"))
            parts = [genome[d["chrom"]][s - 1:e] for s, e in segs]
            if d["strand"] == "-":
                parts = [str(Seq(p).reverse_complement()) for p in parts]
            cds = "".join(parts)
            fc.write(f">{g} biotype={d['biotype']} len={len(cds)} product={d['product']}\n")
            for i in range(0, len(cds), 70):
                fc.write(cds[i:i + 70] + "\n")
            n_cds += 1
            if d["biotype"] == "protein_coding":
                p = str(Seq(cds[:len(cds) // 3 * 3]).translate(table=1)).rstrip("*")
                fp.write(f">{g} len={len(p)} product={d['product']}\n")
                for i in range(0, len(p), 70):
                    fp.write(p[i:i + 70] + "\n")
                n_prot += 1

    # 统计
    bt = collections.Counter(d["biotype"] for d in genes.values())
    n_exon = sum(len(d["exon"]) for d in genes.values())
    n_gene = len(genes)
    coding = [d for d in genes.values() if d["biotype"] == "protein_coding"]
    introns = sum(max(0, len(d["exon"]) - 1) for d in coding)
    multi = sum(1 for d in coding if len(d["exon"]) > 1)
    stats = dict(
        genes=n_gene,
        by_biotype=dict(bt),
        exons=n_exon,
        introns_coding=introns,
        coding_genes_with_intron=multi,
        cds_written=n_cds,
        proteins_written=n_prot,
        transcripts_written=n_tx,
    )
    with open(os.path.join(a.outdir, f"{a.prefix}.stats.json"), "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=2)

    print(json.dumps(stats, ensure_ascii=False, indent=2))
    print("写出:", cds_path)
    print("写出:", prot_path)
    print("写出:", tx_path)


if __name__ == "__main__":
    main()
