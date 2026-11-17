#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""论文数字 / 结构命名 / 提交表格 与最终注释（GTF）的一致性审计。

纯标准库实现（不依赖 Biopython）。只读，不修改注释文件。

用法
  python audit_paper_consistency.py --gtf <gtf> --genome <fa> --tables <tsv...> --outdir <dir>
"""
import argparse
import collections
import json
import os
import re

ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
GID = re.compile(r'BdWA1_(?:nc[0-9]{3}[0-9]{5}|[0-9]{3}(?:pg|g|t|r)[0-9]{5})')

# NCBI 标准遗传密码（table 1），碱基顺序 T,C,A,G
CODON = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
BIDX = {"T": 0, "C": 1, "A": 2, "G": 3}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(s):
    out = []
    for i in range(0, len(s) - 2, 3):
        c = s[i:i + 3]
        try:
            out.append(CODON[16 * BIDX[c[0]] + 4 * BIDX[c[1]] + BIDX[c[2]]])
        except KeyError:
            out.append("X")
    return "".join(out)


def read_fasta(path):
    d, name, buf = {}, None, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name:
                d[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        d[name] = "".join(buf).upper()
    return d


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


def q(v, p):
    if not v:
        return None
    v = sorted(v)
    k = (len(v) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(v) - 1)
    return round(v[lo] + (v[hi] - v[lo]) * (k - lo), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--tables", nargs="*", default=[])
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    genes, order = collections.OrderedDict(), []
    for line in open(a.gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(ATTR.findall(f[8]))
        gid = at.get("gene_id")
        if not gid:
            continue
        if gid not in genes:
            genes[gid] = dict(chrom=f[0], strand=f[6], start=int(f[3]), end=int(f[4]),
                              biotype="", product="", exons=[], cds=[], u5=[], u3=[],
                              tx_id="", pseudo_type="")
            order.append(gid)
        d = genes[gid]
        t = f[2]
        if t == "gene":
            d["biotype"] = at.get("gene_biotype", "")
            d["product"] = at.get("product", "")
            d["pseudo_type"] = at.get("pseudogene_type", "")
        elif t == "transcript":
            d["tx_id"] = at.get("transcript_id", "")
        elif t == "exon":
            d["exons"].append((int(f[3]), int(f[4])))
        elif t == "CDS":
            d["cds"].append((int(f[3]), int(f[4])))
        elif t == "five_prime_utr":
            d["u5"].append((int(f[3]), int(f[4])))
        elif t == "three_prime_utr":
            d["u3"].append((int(f[3]), int(f[4])))

    bt = collections.Counter(d["biotype"] for d in genes.values())
    coding = [g for g, d in genes.items() if d["biotype"] == "protein_coding"]
    pseudo = [g for g, d in genes.items() if d["biotype"] == "pseudogene"]
    lnc = [g for g, d in genes.items() if d["biotype"] == "lncRNA"]
    cds_genes = [g for g, d in genes.items() if d["cds"]]

    n_exon = sum(len(d["exons"]) for d in genes.values())
    n_cds_feat = sum(len(d["cds"]) for d in genes.values())
    introns_all = sum(max(0, len(d["exons"]) - 1) for d in genes.values())
    introns_coding = sum(max(0, len(genes[g]["exons"]) - 1) for g in coding)
    introns_lnc = sum(max(0, len(genes[g]["exons"]) - 1) for g in lnc)

    genome = read_fasta(a.genome)

    ilen, gt_ag, noncanon = [], 0, 0
    for g, d in genes.items():
        if len(d["exons"]) < 2:
            continue
        es = sorted(d["exons"])
        seq = genome.get(d["chrom"], "")
        for i in range(len(es) - 1):
            s, e = es[i][1], es[i + 1][0]
            L = e - s - 1
            ilen.append(L)
            if seq and L >= 4:
                don, acc = seq[s:s + 2], seq[e - 3:e - 1]
                if d["strand"] == "-":
                    don, acc = revcomp(acc), revcomp(don)
                if don == "GT" and acc == "AG":
                    gt_ag += 1
                else:
                    noncanon += 1

    u5_len = [sum(e - s + 1 for s, e in genes[g]["u5"]) for g in coding if genes[g]["u5"]]
    u3_len = [sum(e - s + 1 for s, e in genes[g]["u3"]) for g in coding if genes[g]["u3"]]

    internal_stop, no_stop = [], []
    for g in cds_genes:
        d = genes[g]
        segs = sorted(d["cds"], reverse=(d["strand"] == "-"))
        parts = [genome[d["chrom"]][s - 1:e] for s, e in segs]
        if d["strand"] == "-":
            parts = [revcomp(p) for p in parts]
        cds = "".join(parts)
        prot = translate(cds[:len(cds) // 3 * 3])
        ns = prot[:-1].count("*")
        if ns:
            internal_stop.append((g, ns))
        if not prot.endswith("*"):
            no_stop.append(g)

    described = [g for g in coding if genes[g]["product"]
                 and not genes[g]["product"].lower().startswith("hypothetical")]
    struct_named = [g for g in coding if genes[g]["product"].startswith("Protein structure like")]

    series = collections.Counter()
    for g in genes:
        m = re.match(r'BdWA1_(nc)([0-9]{3})([0-9]{5})$', g)
        if m:
            series["nc"] += 1
            continue
        m = re.match(r'BdWA1_([0-9]{3})(pg|g|t|r)([0-9]{5})$', g)
        if m:
            series[m.group(2)] += 1
        else:
            series["UNMATCHED"] += 1

    # 各 scaffold 内编号是否连续
    by_scaf = collections.defaultdict(list)
    for g in genes:
        m = re.match(r'BdWA1_(?:nc)?([0-9]{3})(pg|g|t|r)([0-9]{5})$', g)
        if m:
            by_scaf[(m.group(1), m.group(2))].append(int(m.group(3)))
    gaps = {f"{k[0]}_{k[1]}": sorted(set(range(1, max(v) + 1)) - set(v))
            for k, v in by_scaf.items() if set(range(1, max(v) + 1)) - set(v)}

    stats = dict(
        n_genes=len(genes),
        by_biotype=dict(bt),
        exons=n_exon,
        cds_features=n_cds_feat,
        genes_with_cds=len(cds_genes),
        introns_total=introns_all,
        introns_coding=introns_coding,
        introns_lncrna=introns_lnc,
        intron_median=q(ilen, 0.5),
        intron_iqr=[q(ilen, 0.25), q(ilen, 0.75)],
        intron_min=min(ilen), intron_max=max(ilen),
        intron_n_lt30=sum(1 for x in ilen if x < 30),
        intron_gt_ag_pct=round(pct(gt_ag, len(ilen)), 1),
        intron_noncannon=noncanon,
        n_introns_measured=len(ilen),
        utr5_genes=len(u5_len), utr5_rows=sum(len(d["u5"]) for d in genes.values()),
        utr3_genes=len(u3_len), utr3_rows=sum(len(d["u3"]) for d in genes.values()),
        utr5_median=q(u5_len, 0.5), utr5_q1=q(u5_len, 0.25), utr5_q3=q(u5_len, 0.75),
        utr3_median=q(u3_len, 0.5), utr3_q1=q(u3_len, 0.25), utr3_q3=q(u3_len, 0.75),
        utr_median_ratio=round((q(u3_len, 0.5) or 0) / (q(u5_len, 0.5) or 1), 2),
        utr5_lt50_pct=round(pct(sum(1 for x in u5_len if x < 50), len(u5_len)), 1),
        internal_stop_genes=len(internal_stop),
        internal_stop_detail=internal_stop[:10],
        cds_without_terminal_stop=len(no_stop),
        described=len(described), described_pct=round(pct(len(described), len(coding)), 1),
        structure_named=len(struct_named),
        id_series=dict(series),
        id_gaps=gaps,
        duplicate_ids=[g for g, c in collections.Counter(order).items() if c > 1],
        genes_without_product=[g for g, d in genes.items() if not d["product"]][:10],
        genes_without_transcript_id=[g for g, d in genes.items() if not d["tx_id"]][:10],
        pseudo_types=dict(collections.Counter(genes[g]["pseudo_type"] or "NA" for g in pseudo)),
        pseudo_with_cds=sum(1 for g in pseudo if genes[g]["cds"]),
        genome_seqs=len(genome),
        genome_bp=sum(len(v) for v in genome.values()),
    )

    # ---- 结构表 product 一致性 ----
    struct_check = dict(checked=0, match=0, mismatch=[], not_in_gtf=[], gtf_struct_named=len(struct_named))
    cand = [x for x in a.tables if x.endswith("structure_annotations.tsv")]
    if cand:
        rows = [l.rstrip("\n").split("\t") for l in open(cand[0], encoding="utf-8") if l.strip()]
        hdr = rows[0]
        ig = hdr.index("gene_id") if "gene_id" in hdr else 0
        ip = hdr.index("new_product") if "new_product" in hdr else None
        tbl = set()
        for r in rows[1:]:
            g = r[ig]
            tbl.add(g)
            struct_check["checked"] += 1
            gt = genes.get(g, {}).get("product")
            if g not in genes:
                struct_check["not_in_gtf"].append(g)
            elif ip is not None and gt == r[ip]:
                struct_check["match"] += 1
            else:
                struct_check["mismatch"].append([g, r[ip] if ip is not None else None, gt])
        struct_check["in_gtf_not_in_table"] = sorted(set(struct_named) - tbl)

    # ---- 表格基因编号核对 ----
    id_check = []
    all_ids = set(genes)
    for tf in a.tables:
        if not os.path.exists(tf):
            continue
        txt = open(tf, encoding="utf-8", errors="replace").read()
        found = set(GID.findall(txt))
        bad = sorted(x for x in found if x not in all_ids)
        id_check.append(dict(table=os.path.basename(tf), n_ids=len(found),
                             n_bad=len(bad), examples=bad[:8]))

    out = dict(stats=stats, struct=struct_check, table_ids=id_check)
    with open(os.path.join(a.outdir, "audit_result.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
