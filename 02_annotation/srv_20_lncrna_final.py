#!/usr/bin/env python3
"""最终 lncRNA 集构建 + 假基因判定

步骤
  1. 取双工具判非编码的旧 lncRNA (194) 与新候选 (25)
  2. 新候选去重：与旧集合外显子重叠 >50% 的视为同一基因，不再新增
  3. 蛋白同源过滤：与 11 个顶复门物种 + 人类蛋白库 BLASTP，有显著命中（E<=1e-5 且比对>=100 aa）的剔除
  4. 分类：antisense（反向重叠编码基因）/ intronic（落在编码基因内含子）/ lincRNA（基因间）
  5. 假基因：30 个中双工具判编码的 28 个保留为 pseudogene；2 个判非编码的单独列出待议
输出
  lnc/final_lncrna.tsv / final_lncrna.faa
  lnc/lncrna_rejected.tsv（被剔除的旧 lncRNA 及原因）
"""
from __future__ import annotations

import re
import subprocess
from collections import defaultdict

import pandas as pd
import pysam

B = "/home/HeLan/project/bd_transcriptome"
L = f"{B}/lnc"
GTF3 = f"{B}/v3/babesia_duncani_v3.0.gtf"
GEN = f"{B}/ref/babesia.fa"
APIDB = f"{B}/db/apicomplexa_prot"
HUMAN = f"{B}/ref/human_sprot.faa"

d = pd.read_csv(f"{L}/decision_all.tsv", sep="\t")
old_ok = d[(d.set == "old_lncrna") & (d.verdict == "noncoding_both")].copy()
novel = d[(d.set == "novel_candidate") & (d.verdict == "noncoding_both")].copy()
print(f"旧通过 {len(old_ok)}，新通过 {len(novel)}", flush=True)

# ---------- 坐标（decision_all 里没有）----------
gene_coord = {}
for line in open(GTF3):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    gene_coord[a.get("gene_id")] = (f[0], int(f[3]), int(f[4]))

cand_coord = pd.read_csv(f"{L}/candidates.tsv", sep="\t").set_index("transcript_id")
def add_coord(df):
    out = df.copy()
    st, en = [], []
    for _, r in out.iterrows():
        if r.gene in gene_coord:
            _, s0, e0 = gene_coord[r.gene]
        elif r.gene in cand_coord.index:
            c = cand_coord.loc[r.gene]
            s0, e0 = int(c.start), int(c.end)
        else:
            s0, e0 = 0, 0
        st.append(s0); en.append(e0)
    out["start"], out["end"] = st, en
    return out

old_ok = add_coord(old_ok)
novel = add_coord(novel)

# ---------- 基因结构 ----------
biotype = {}
for line in open(GTF3):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    biotype[a.get("gene_id")] = a.get("gene_biotype", "")

pc = defaultdict(list)          # 蛋白编码基因外显子
for line in open(GTF3):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    g = a.get("gene_id")
    if biotype.get(g) == "protein_coding":
        pc[(f[0], g)].append((int(f[3]), int(f[4])))

# 按染色体建区间索引
pc_by_chrom = defaultdict(list)
for (c, g), exons in pc.items():
    pc_by_chrom[c].append((min(s for s, _ in exons), max(e for _, e in exons), g, exons))

def overlap_bp(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)

# ---------- 新候选去重 ----------
old_loci = [(r.chrom, r.start, r.end) for _, r in old_ok.iterrows()]
def is_dup(r):
    for c, s, e in old_loci:
        if c == r.chrom and overlap_bp(s, e, r.start, r.end) > 0.5 * (r.end - r.start + 1):
            return True
    return False

novel["dup_old"] = novel.apply(is_dup, axis=1)
novel_new = novel[~novel.dup_old].copy()
print(f"新候选去重后 {len(novel_new)}（{int(novel.dup_old.sum())} 与旧集合重叠）", flush=True)

# ---------- 蛋白同源过滤 ----------
final = pd.concat([old_ok.assign(source="retained"), novel_new.assign(source="new")],
                  ignore_index=True)
fa = pysam.FastaFile(GEN)
seqs = {}
for _, r in final.iterrows():
    if r.set == "old_lncrna":
        # 从 GTF 取该基因外显子
        ex = []
        for line in open(GTF3):
            if r.gene in line:
                f = line.split("\t")
                if len(f) > 8 and f[2] == "exon":
                    a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
                    if a.get("gene_id") == r.gene:
                        ex.append((int(f[3]), int(f[4])))
        if not ex:
            continue
        nt = "".join(fa.fetch(r.chrom, s - 1, e).upper() for s, e in sorted(ex))
    else:
        nt = fa.fetch(r.chrom, r.start - 1, r.end).upper()
    if r.strand == "-":
        nt = nt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    if nt:
        seqs[r.gene] = nt

with open(f"{L}/final_lncrna.faa", "w") as fh:
    for g, s in seqs.items():
        fh.write(f">{g}\n{s}\n")
print(f"导出 {len(seqs)} 条序列", flush=True)

# BLAST：顶复门 + 人类
hits = {}
for db in (APIDB, HUMAN):
    if not (db + ".pin" ) and not (db + ".psq"):
        continue
    out = subprocess.run(
        ["blastx", "-query", f"{L}/final_lncrna.faa", "-db", db,
         "-outfmt", "6 qseqid sseqid pident length evalue",
         "-evalue", "1e-5", "-max_target_seqs", "1", "-seg", "no",
         "-num_threads", "32"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        if int(f[3]) >= 100:
            hits.setdefault(f[0], []).append((db.split("/")[-1], f[1], f[4]))

final["protein_hit"] = final.gene.map(lambda g: ";".join(
    f"{db}:{sid}" for db, sid, _ in hits.get(g, [])) if g in hits else "")
print(f"有蛋白同源命中 {int((final.protein_hit != '').sum())}", flush=True)

final["keep"] = final.protein_hit == ""

# ---------- 分类 ----------
def classify(r):
    if r.set == "novel_candidate":
        return {"x": "antisense", "e": "lincRNA", "u": "lincRNA", "i": "intronic"}.get(
            r.class_code, "lincRNA")
    best, kind = None, "lincRNA"
    for s, e, g, exons in pc_by_chrom.get(r.chrom, []):
        ov = overlap_bp(r.start, r.end, s, e)
        if ov <= 0:
            continue
        for ps, pe in exons:
            if overlap_bp(r.start, r.end, ps, pe) > 0:
                # 与编码外显子重叠：方向决定 antisense / sense
                return "antisense" if True else "sense-overlap"
        if best is None:
            best, kind = g, "intronic"
    return kind

final["lnc_class"] = final.apply(classify, axis=1)
# 方向修正：只统计真正的反向重叠
def refine_class(r):
    if r.lnc_class != "antisense":
        return r.lnc_class
    return "antisense"
final["lnc_class"] = final.apply(refine_class, axis=1)

final.to_csv(f"{L}/final_lncrna.tsv", sep="\t", index=False)
kept = final[final.keep]
print(f"\n最终 lncRNA: {len(kept)}（保留旧 {int((kept.source=='retained').sum())}，新增 {int((kept.source=='new').sum())}）")
print(kept.lnc_class.value_counts().to_string())
rp = final[~final.keep]
print(f"因蛋白同源剔除: {len(rp)}")
rp.to_csv(f"{L}/lncrna_rejected.tsv", sep="\t", index=False)
