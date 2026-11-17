#!/usr/bin/env python3
"""把可信的结构邻居转成注释：product 前加 "Protein structure like ..."

三级结构置信度（见 srv_46 的判定）
  global_high  整条平均 pLDDT >= 70          → 整体可信
  global_ok    整条平均 pLDDT 50-70          → 整体基本可信
  domain_ok    整条平均 < 50，但含 >=40 aa   → 仅结构域可信
               且域内均值 >=70 的高置信片段
  low          以上都不满足                  → 不采用

其余入选条件
  序列对得上：结构序列与 v3.8 翻译蛋白一致性 >= 95%
  比对可信  ：Foldseek top1 prob >= 0.8、evalue <= 1e-5，邻居描述有实际名字
  domain_ok 额外要求：高置信结构域残基数 >= 比对覆盖长度的一半，避免"整条
             只有一段乱码碰上了"的情况
只改功能未知的基因；已有名称不动。
product_source：整体可信写 structure-based，结构域级写 structure-domain-based。

输出
  <outdir>/babesia_duncani_v3.9.gtf
  <outdir>/structure_annotations.tsv        采用者明细
  <outdir>/structure_confidence_all.tsv     全部被查结构的判定表（含拒因）
  <outdir>/annotation_summary.tsv           统计
"""
import argparse
import collections
import csv
import os
import re

A = re.compile(r'(\w+) "([^"]*)"')
UNKNOWN = {"hypothetical protein", "uncharacterized protein", "unknown function",
           "protein of unknown function", ""}

# 名字信息量分级：generic 只说明"有信号肽/是膜蛋白"，不含功能信息
GENERIC = re.compile(
    r"^(signal peptide[ -]containing protein|membrane protein.*|integral membrane protein|"
    r"multi-pass transmembrane protein|checkpoint protein|transmembrane protein.*|"
    r"transmembrane and coiled-coil domain-containing protein)$", re.I)
REVIEW = re.compile(r"viral|xylem|plant|arabidopsis|toxoplasma", re.I)


def seq_domain_class(attrs):
    """判断"序列层面能不能拿到这个域/功能"。

    named_domain  PFAM 里有已知域（非 DUF）→ 序列已给出域
    described    只有 DUF，但 pfam 描述写了功能名（如 GET4/CCDC47/MINDY）→ 序列已给出功能
    duf_only     只有 DUF 且描述为 unknown function → 序列只给出"未知功能域"
    ko_only      无 PFAM 但 eggNOG 给了 KEGG KO → 序列已归入直系同源组
    none         序列层面没有任何域注释
    """
    doms = [d for d in (attrs.get("PFAM", "") or "").split("|") if d]
    named = [d for d in doms if not d.startswith("DUF")]
    if named:
        return "named_domain", "|".join(named)
    desc = attrs.get("pfam", "") or ""
    if doms and desc and "unknown function" not in desc.lower():
        return "described", desc
    if doms:
        return "duf_only", "|".join(doms)
    if attrs.get("KEGG_ko"):
        return "ko_only", attrs["KEGG_ko"]
    return "none", ""


def name_quality(desc):
    if GENERIC.match(desc.strip()):
        return "generic"
    if "domain-containing protein" in desc.lower() or "domain containing protein" in desc.lower():
        return "domain_only"
    return "specific"


def is_unknown(p):
    p = (p or "").strip()
    return p in UNKNOWN or p.startswith("DUF") or p.lower().startswith("uncharacterized")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--foldseek", required=True)
    ap.add_argument("--foldseek-extra", default="")
    ap.add_argument("--domain", required=True, help="srv_46 输出的 domain_plddt.tsv")
    ap.add_argument("--seqcheck", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-plddt", type=float, default=50)
    ap.add_argument("--min-prob", type=float, default=0.8)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    ap.add_argument("--min-identity", type=float, default=95)
    ap.add_argument("--domain-cover", type=float, default=0.5,
                    help="结构域级：高置信域残基 / 比对覆盖长度 的最低比例")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    to_old = {r["new_id"]: r["old_id"] for r in
              csv.DictReader(open(a.map, encoding="utf-8"), delimiter="\t")}
    dom = {r["pdb_id"]: r for r in
           csv.DictReader(open(a.domain, encoding="utf-8"), delimiter="\t")}
    ident = {r["pdb_id"]: float(r["identity"]) for r in
             csv.DictReader(open(a.seqcheck, encoding="utf-8"), delimiter="\t")
             if r.get("identity")}

    gtf_prod, gtf_attr = {}, {}
    for line in open(a.gtf, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "gene":
            attrs = dict(A.findall(f[8]))
            gid = attrs.get("gene_id", "")
            gtf_prod[gid] = attrs.get("product", "")
            gtf_attr[gid] = attrs

    fs = {}
    for path in [a.foldseek] + ([a.foldseek_extra] if a.foldseek_extra else []):
        if not path or not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
            if r.get("pdb_id"):
                fs.setdefault(r["pdb_id"], r)

    new2old = {r["new_id"]: r["old_id"] for r in
               csv.DictReader(open(a.map, encoding="utf-8"), delimiter="\t")
               if r.get("new_id")}
    old2gene = {}
    for g, o in to_old.items():
        old2gene.setdefault(o, g)

    changes, decisions = [], []
    reasons = collections.Counter()
    prefilter = collections.Counter()
    excluded = []
    for pid, r in fs.items():
        # 结构文件已按 v3.9 基因号命名，优先用结果表里的 gene_id；兼容旧表回退映射
        g = (r.get("gene_id") or "").strip() or old2gene.get(pid)
        # 入口条件：序列层面必须给不出功能名或域。已有名字、或 PFAM/KO 已给出域的基因
        # 不属于结构命名的候选池，在建池阶段就排除，不进入判定表。
        if not g:
            prefilter["无基因号"] += 1
            continue
        if not is_unknown(gtf_prod.get(g, "")):
            prefilter["序列已有功能名"] += 1
            excluded.append(dict(gene_id=g, pdb_id=pid,
                                 product=gtf_prod.get(g, ""),
                                 seq_domain_class=seq_domain_class(gtf_attr.get(g, {}))[0],
                                 seq_domain_evidence=seq_domain_class(gtf_attr.get(g, {}))[1],
                                 excluded_because="product name already assigned"))
            continue
        seq_cls, seq_ev = seq_domain_class(gtf_attr.get(g, {}))
        if seq_cls not in ("none", "duf_only"):
            prefilter["序列已给出域/功能"] += 1
            excluded.append(dict(gene_id=g, pdb_id=pid,
                                 product=gtf_prod.get(g, ""),
                                 seq_domain_class=seq_cls,
                                 seq_domain_evidence=seq_ev,
                                 excluded_because="sequence gave domain or function"))
            continue
        d = dom.get(pid, {})
        rec = dict(gene_id=g or "", pdb_id=pid,
                   pdb_id_legacy=new2old.get(pid, pid),
                   old_product=gtf_prod.get(g, "") if g else "",
                   tier=d.get("tier", "no_structure"),
                   mean_pLDDT=d.get("mean_pLDDT", ""),
                   frac_ge70=d.get("frac_ge70", ""),
                   dom_frac=d.get("dom_frac", ""),
                   longest_dom_len=d.get("longest_dom_len", ""),
                   longest_dom_mean=d.get("longest_dom_mean", ""),
                   n_aa=r.get("length_aa", ""),
                   aln_cov=r.get("top1_cov", ""),
                   domain_support="",
                   neighbor=(r.get("top1_desc") or "").strip(),
                   name_quality="", review_flag="",
                   seq_domain_class="", seq_domain_evidence="",
                   db=r.get("top1_db", ""), prob=r.get("top1_prob", ""),
                   evalue=r.get("top1_evalue", ""),
                   identity=ident.get(pid, ""), decision="reject", reason="",
                   new_product="")
        rec["seq_domain_class"], rec["seq_domain_evidence"] = seq_cls, seq_ev
        rec["name_quality"] = name_quality(rec["neighbor"]) if rec["neighbor"] else ""
        rec_quality = rec["name_quality"]
        reason = ""
        if r.get("status") != "ok":
            reason = "Foldseek 无命中"
        else:
            try:
                prob = float(r["top1_prob"])
                ev = float(r["top1_evalue"])
            except (ValueError, KeyError):
                reason = "缺乏概率/期望值"
            else:
                desc = (r.get("top1_desc") or "").strip()
                idn = ident.get(pid, 100.0)
                if idn < a.min_identity:
                    reason = f"序列一致性<{a.min_identity:.0f}%"
                elif prob < a.min_prob:
                    reason = f"prob<{a.min_prob}"
                elif ev > a.max_evalue:
                    reason = f"evalue>{a.max_evalue:g}"
                elif is_unknown(desc):
                    reason = "邻居本身未命名"
                elif seq_cls not in ("none", "duf_only"):
                    rec["seq_gate"] = "1"
                    reason = f"序列已给出域/功能（{seq_cls}）"
                elif rec_quality == "generic":
                    rec["seq_gate"] = "1"
                    reason = "邻居名无功能信息（序列可预测）"
                else:
                    rec["seq_gate"] = "1"
                    tier = d.get("tier", "")
                    if tier in ("global_high", "global_ok"):
                        pass
                    elif tier == "domain_ok":
                        n_aa = float(r.get("length_aa") or 0) or float(d.get("n_res") or 0)
                        cov = float(r.get("top1_cov") or 0) / 100.0
                        aln_len = max(1.0, cov * n_aa)
                        sup = float(d.get("dom_res_total") or 0) / aln_len
                        rec["domain_support"] = round(sup, 2)
                        if sup < a.domain_cover:
                            reason = "高置信结构域不足以支撑该比对"
                    else:
                        reason = "结构置信度不足"
                    if not reason:
                        tier_label = "structure-based" if tier.startswith("global") else \
                                     "structure-domain-based"
                        rec["decision"] = "accept"
                        rec["tier_label"] = tier_label
                        rec["new_product"] = f"Protein structure like {desc}"
                        changes.append(dict(
                            gene_id=g, pdb_id=pid,
                            pdb_id_legacy=new2old.get(pid, pid),
                            old_product=gtf_prod.get(g, ""),
                            new_product=rec["new_product"], tier=tier,
                            tier_label=tier_label,
                            pLDDT=d.get("mean_pLDDT", ""),
                            dom_frac=d.get("dom_frac", ""),
                            longest_dom_len=d.get("longest_dom_len", ""),
                            longest_dom_mean=d.get("longest_dom_mean", ""),
                            neighbor=desc, name_quality=name_quality(desc),
                            seq_domain_class=seq_cls, seq_domain_evidence=seq_ev,
                            review_flag="1" if REVIEW.search(desc) else "",
                            neighbor_id=r.get("top1_target", ""),
                            db=r.get("top1_db", ""), prob=prob, evalue=ev,
                            seqid=r.get("top1_seqid", ""), cov=r.get("top1_cov", ""),
                            identity=idn))
        rec["name_quality"] = name_quality(rec["neighbor"]) if rec["neighbor"] else ""
        rec["review_flag"] = "1" if REVIEW.search(rec["neighbor"] or "") else ""
        rec["reason"] = reason
        if reason:
            reasons[reason] += 1
        decisions.append(rec)

    # 写入 GTF
    by_gene = {c["gene_id"]: c for c in changes}
    out_lines, applied = [], collections.Counter()
    for line in open(a.gtf, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            out_lines.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            out_lines.append(line)
            continue
        attrs = dict(A.findall(f[8]))
        c = by_gene.get(attrs.get("gene_id", ""))
        if c:
            old = attrs.get("product", "")
            f[8] = f[8].replace(f'product "{old}"', f'product "{c["new_product"]}"')
            f[8] = f[8].replace(f'description "{attrs.get("description", old)}"',
                                f'description "{c["new_product"]}"')
            f[8] += f'; product_source "{c["tier_label"]}"'
            applied[c["tier_label"]] += 1
        out_lines.append("\t".join(f) + "\n")

    # 规范化：transcript 行补 transcript_id（v3.6 起一直缺，gffread/NCBI 都要求）
    fixed = 0
    for i, line in enumerate(out_lines):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        attrs = dict(A.findall(f[8]))
        if attrs.get("transcript_id"):
            continue
        gid = attrs.get("gene_id", "")
        if not gid:
            continue
        f[8] = f'gene_id "{gid}"; transcript_id "{gid}";' + \
               f[8].split(f'gene_id "{gid}";', 1)[-1]
        out_lines[i] = "\t".join(f) + "\n"
        fixed += 1

    gtf_out = os.path.join(a.outdir, "babesia_duncani_v3.9.gtf")
    open(gtf_out, "w", encoding="utf-8", newline="\n").write("".join(out_lines))
    print(f"transcript 行补 transcript_id：{fixed}")

    def dump(name, rows):
        if not rows:
            return
        keys = list(rows[0])
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(os.path.join(a.outdir, name), "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    audit = [r for r in decisions if r.get("seq_gate") == "1"]
    audit.sort(key=lambda r: (r["decision"] != "accept", r["gene_id"]))
    dump("sequence_domain_audit.tsv", audit)
    if excluded:
        dump("structure_entry_excluded.tsv", excluded)
    dump("structure_annotations.tsv", changes)
    dump("structure_confidence_all.tsv", decisions)

    # 分层判定汇总（论文补充表/图注用）
    tier_note = {
        "global_high": "整条平均 pLDDT>=70，整体可信",
        "global_ok": "整条平均 pLDDT 50-70，整体基本可信",
        "domain_ok": "整条<50，但含>=40 aa 的高置信结构域",
        "low": "整条<50 且无可信结构域",
        "no_structure": "无结构文件",
    }
    tier_stat = collections.defaultdict(lambda: collections.Counter())
    for d in decisions:
        s = tier_stat[d["tier"]]
        s["n"] += 1
        s[d["decision"]] += 1
    with open(os.path.join(a.outdir, "structure_tier_summary.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("tier\tn_structures\taccepted\trejected\tdefinition\n")
        for k in ("global_high", "global_ok", "domain_ok", "low", "no_structure"):
            if k not in tier_stat:
                continue
            v = tier_stat[k]
            fh.write(f"{k}\t{v['n']}\t{v['accept']}\t{v['reject']}\t{tier_note[k]}\n")

    summary = [("采用（整体可信 structure-based）", applied["structure-based"]),
               ("采用（结构域级 structure-domain-based）", applied["structure-domain-based"]),
               ("合计采用", sum(applied.values()))]
    for k, v in collections.Counter(c["name_quality"] for c in changes).most_common():
        summary.append((f"  其中邻居名信息量={k}", v))
    summary.append(("  需人工复核（远缘物种名）",
                    sum(1 for c in changes if c["review_flag"] == "1")))
    summary.append(("入口前排除（序列已有功能名）", prefilter["序列已有功能名"]))
    summary.append(("入口前排除（序列已给出域/功能）", prefilter["序列已给出域/功能"]))
    summary.append(("进入结构判定的候选", len(decisions)))
    summary.append(("被拒", sum(reasons.values())))
    for k, v in reasons.most_common():
        summary.append((f"  拒因：{k}", v))
    with open(os.path.join(a.outdir, "annotation_summary.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        for k, v in summary:
            fh.write(f"{k}\t{v}\n")
    print("结构注释写入", sum(applied.values()), "个基因",
          dict(applied), "->", gtf_out)
    for k, v in reasons.most_common(6):
        print(f"  拒因 {k}: {v}")


if __name__ == "__main__":
    main()
