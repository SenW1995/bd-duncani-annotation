#!/usr/bin/env python3
"""蛋白水平注释核查：11 个顶复门参考蛋白组 + 人类 vs V1/V2 两套注释蛋白集。

判据（统一）：pident >= 30 且比对长度 >= 50 aa（BLAST 预筛 evalue <= 1e-5，-seg no）
  - 有同源   : 命中且满足上述阈值
  - 无同源   : 无任何命中或低于阈值
输出全部为基因水平，不做坐标推断。
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

W = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("WORK_DIR", ".")
O = sys.argv[2] if len(sys.argv) > 2 else os.path.join(W, "tables_prot")
os.makedirs(O, exist_ok=True)

COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp"]

SPECIES = {
    "pfal": ("Plasmodium falciparum", "Haemosporidia"),
    "tgondii": ("Toxoplasma gondii", "Coccidia"),
    "ncaninum": ("Neospora caninum", "Coccidia"),
    "etenella": ("Eimeria tenella", "Coccidia"),
    "ccayetanensis": ("Cyclospora cayetanensis", "Coccidia"),
    "cparvum": ("Cryptosporidium parvum", "Cryptosporidia"),
    "chominis": ("Cryptosporidium hominis", "Cryptosporidia"),
    "bbovis": ("Babesia bovis", "Piroplasmida"),
    "bmicroti": ("Babesia microti", "Piroplasmida"),
    "tannulata": ("Theileria annulata", "Piroplasmida"),
    "tparva": ("Theileria parva", "Piroplasmida"),
    "human": ("Homo sapiens", "Host"),
}
APICO = [k for k in SPECIES if k != "human"]
MIN_PID, MIN_LEN = 30.0, 50   # 一致性>=30%，比对长度>=50 aa（避免用 qcovhsp，它对短蛋白有方向偏差）


def load(tag: str) -> pd.DataFrame:
    p = os.path.join(W, "out", f"{tag}.tsv")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return pd.DataFrame(columns=COLS)
    return pd.read_csv(p, sep="\t", names=COLS)


def best_per_subject(df: pd.DataFrame, apply_threshold: bool = True) -> pd.DataFrame:
    """每个注释基因的最佳命中（按 pident，同分取 bitscore 高者）。

    apply_threshold=False 时不筛阈值，用于记录宿主最佳一致性（边界情形可判读）。
    """
    cols = ["sseqid", "pident", "length", "qcovhsp"]
    if df.empty:
        return pd.DataFrame(columns=cols)
    d = df
    if apply_threshold:
        d = d[(d.pident >= MIN_PID) & (d.length >= MIN_LEN)]
    if d.empty:
        return pd.DataFrame(columns=cols)
    d = d.sort_values(["pident", "bitscore"], ascending=False)
    return d.drop_duplicates("sseqid")[cols]


v1 = pd.read_csv(os.path.join(W, "in", "v1_cds_qc.tsv"), sep="\t")
v2 = pd.read_csv(os.path.join(W, "in", "v2_cds_qc.tsv"), sep="\t")
v2g = pd.read_csv(os.path.join(W, "in", "v2_genes.tsv"), sep="\t")
v2 = v2.merge(v2g, left_on="cds_id", right_on="gene_id", how="left")
CHROM = {"CM053815.1", "CM053816.1", "CM053817.1", "CM053818.1", "CM053819.1"}
v1["on_chromosome"] = v1.cds_id.str.extract(r"lcl\|([^_]+)_cds")[0].isin(CHROM)

ann = {"v1prot": v1.rename(columns={"protein_id": "gid"}),
       "v2prot": v2.rename(columns={"gene_id": "gid"})}

rows_recovery = []
for tag in ("v1prot", "v2prot"):
    key = "gid"
    sup = pd.DataFrame({key: ann[tag][key].astype(str)})
    for sp in list(SPECIES):
        is_host = sp == "human"
        b = best_per_subject(load(f"{sp}__{tag}"), apply_threshold=not is_host)
        sup = sup.merge(b.rename(columns={"sseqid": key,
                                          "pident": f"pid_{sp}",
                                          "length": f"len_{sp}",
                                          "qcovhsp": f"cov_{sp}"}),
                        on=key, how="left")
    pid_cols = [f"pid_{s}" for s in SPECIES]
    hit_cols = [f"pid_{s}" for s in APICO]
    sup["n_species"] = sup[hit_cols].notna().sum(axis=1)
    sup["species"] = sup[hit_cols].apply(
        lambda r: ",".join(s for s in APICO if pd.notna(r[f"pid_{s}"])), axis=1)
    sup["best_pid"] = sup[hit_cols].max(axis=1).round(1)
    # 宿主：记录最佳一致性本身，并按 30%/50aa 判定是否为真同源
    sup["human_pid"] = sup["pid_human"].round(1)
    sup["human_aln_len"] = sup["len_human"]
    sup["human_hit"] = (sup["pid_human"] >= MIN_PID) & (sup["len_human"] >= MIN_LEN)
    sup = sup.drop(columns=[f"cov_{s}" for s in SPECIES] + ["pid_human", "len_human"])

    meta = ann[tag].copy()
    extra = [c for c in ("gene", "length_aa", "internal_stops", "starts_with_M",
                         "ends_with_stop", "scaffold", "start", "end", "strand",
                         "biotype", "product", "on_chromosome", "cds_id") if c in meta.columns]
    out = meta[[key] + extra].merge(sup, on=key, how="left")
    out = out.sort_values(["n_species", "best_pid"], ascending=False)
    fname = "01_v1_gene_support.tsv" if tag == "v1prot" else "02_v2_gene_support.tsv"
    out.to_csv(os.path.join(O, fname), sep="\t", index=False)

    # 每个物种：有多少参考蛋白能在该注释中找到伙伴（query 侧）
    for sp in list(SPECIES):
        d = load(f"{sp}__{tag}")
        n_tot = sum(1 for line in open(os.path.join(W, "ref", f"{sp}.faa")) if line.startswith(">"))
        if d.empty:
            qn = 0
        else:
            dd = d[(d.pident >= MIN_PID) & (d.length >= MIN_LEN)]
            qn = dd.qseqid.nunique()
        rows_recovery.append(dict(annotation=tag, species=sp, name=SPECIES[sp][0],
                                  lineage=SPECIES[sp][1], proteome=n_tot,
                                  matched=qn, pct=round(100 * qn / n_tot, 1)))

rec = pd.DataFrame(rows_recovery)
rec.to_csv(os.path.join(O, "03_species_recovery.tsv"), sep="\t", index=False)

# 药靶候选：V2 基因为主（报告用编号），V1 交叉引用
v2s = pd.read_csv(os.path.join(O, "02_v2_gene_support.tsv"), sep="\t")
v1s = pd.read_csv(os.path.join(O, "01_v1_gene_support.tsv"), sep="\t")
mp = pd.read_csv(os.path.join(W, "in", "map_v1_v2.tsv"), sep="\t")
v1_to_v2 = dict(zip(mp.v1_protein_id, mp.v2_gene_id))
v2_to_v1 = dict(zip(mp.v2_gene_id, mp.v1_protein_id))

cand = v2s[(v2s.n_species >= 6) & (~v2s.human_hit) & v2s.biotype.eq("protein_coding")].copy()
cand["v1_protein_id"] = cand.gid.map(v2_to_v1).fillna("")
cand["host_class"] = np.where(
    cand.human_pid.isna() | (cand.human_pid < 25), "no_host_homolog", "borderline_25_30")
cand = cand.sort_values(["n_species", "best_pid"], ascending=False)
cand.to_csv(os.path.join(O, "04_drug_target_candidates.tsv"), sep="\t", index=False)

conserved = v2s[(v2s.n_species >= 6)].copy()
summary = {
    "V1_genes": int(len(v1)), "V2_genes": int(len(v2)),
    "V1_on_chromosome": int(v1.on_chromosome.sum()),
    "V1_on_contigs": int((~v1.on_chromosome).sum()),
    "V2_internal_stop": int((v2.internal_stops > 0).sum()),
    "V1_median_aa": int(v1.length_aa.median()), "V2_median_aa": int(v2.length_aa.median()),
    "V1_no_ortholog": int((v1s.n_species == 0).sum()),
    "V2_no_ortholog": int((v2s.n_species == 0).sum()),
    "V2_conserved_6plus": int(len(conserved)),
    "V2_conserved_no_human": int(len(cand)),
    "V1_conserved_6plus": int((v1s.n_species >= 6).sum()),
    "V1_conserved_no_human": int(((v1s.n_species >= 6) & (~v1s.human_hit)).sum()),
}
pd.Series(summary).to_csv(os.path.join(O, "summary.tsv"), sep="\t", header=False)
print(pd.Series(summary).to_string())
print()
print(rec.pivot(index="species", columns="annotation", values="pct").to_string())
