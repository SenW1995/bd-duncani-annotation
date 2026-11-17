#!/usr/bin/env python3
"""把 GTF 的 product 属性清洗成 GenBank 能接受的短名词短语。

NCBI 的 discrepancy 报告对 /product 有一票否决项（FATAL），本脚本处理其中
可以在不改动基因模型的前提下修好的部分：

  contains '. '            整段 UniProt 功能描述被当成产物名
  starts with '-'          名字前多了连字符
  'ase' 截断               "cysteine ase" 实为 "cysteine protease"
  unbalanced brackets      括号没配对
  contains 'genome'        KEGG 通路名误当产物名
  equals '16S ribosomal RNA' 等  rRNA 基因被写成 mRNA（由 gtf_to_gff3 处理类型）

不改任何坐标、不改 biotype，只改 product 字符串。
"""
import argparse
import re
import sys

# 逐条人工确认过的替换
EXACT = {
    "With signal peptide and GMGPP repeat at": "GMGPP repeat protein",
    "-choline ethanolamine kinase": "Choline ethanolamine kinase",
    "Translation initiation factor (SUI1": "Translation initiation factor SUI1",
    "Protein structure like Xylem cysteine ase 1-like protein, putative":
        "Protein structure like Xylem cysteine proteinase 1-like protein, putative",
    "Modulation by host of viral RNA genome replication": "hypothetical protein",
    # tRNAscan 未定异亮，但同源模型给的是 Lys（note: isotype_cm Lys）
    "tRNA-Undet": "tRNA-Lys",
}

SUFFIX_WORDS = ("protein", "family", "subunit", "domain", "enzyme", "subfamily")


def clean(p, biotype=""):
    """返回清洗后的 product。"""
    if not p:
        return "lncRNA" if biotype == "lncRNA" else ("hypothetical protein" if biotype else "")
    s = " ".join(p.split())

    if s in EXACT:
        s = EXACT[s]

    # 以 "Belongs to the X superfamily. Y family" 开头的，取后半段的家族名
    if s.startswith("Belongs to"):
        parts = [x.strip() for x in s.split(". ") if x.strip()]
        if len(parts) >= 2 and len(parts[-1]) >= 8:
            s = parts[-1]
        else:
            s = re.sub(r"^Belongs to (the )?", "", parts[0])
            s = re.sub(r"\s*\(TC [^)]*\)", "", s)

    # 整段描述只留第一句
    if ". " in s:
        s = s.split(". ")[0]

    s = s.rstrip(" .;,")

    # 前导符号（GenBank 产物名惯用小写，不做首字母大写）
    s = re.sub(r"^[-–—:;,\s]+", "", s)

    # Pfam/结构域 ID 当名字
    m = re.match(r"^Pfam:(\S+)$", s)
    if m:
        dom = m.group(1).replace("_", " ")
        s = f"{dom} domain-containing protein"
    else:
        s = s.replace("_", " ")

    # 括号配对（必须在截断之后，截断很容易把右括号切掉）
    # cysteine ase -> cysteine protease（只此一处）
    s = s.replace("cysteine ase", "cysteine protease")

    # 断词修复："iron- sulfur" 这类连字符后多了空格
    s = re.sub(r"([A-Za-z])- ([a-z])", r"\1-\2", s)

    # 过长截断到词边界
    if len(s) > 180:
        s = s[:180].rsplit(" ", 1)[0].rstrip(" .,;:-")

    # 截断后可能留下悬空的虚词
    s = re.sub(r"(?:\s+(?:and|or|of|with|the|a|an|at|for|in|into|to|from|by|on|as|"
               r"that|which|where|is|are|was|were|be|been|both|its|their))+$",
               "", s, flags=re.I).rstrip(" .,;:-")

    # 截断后再修一次括号
    if s.count("(") != s.count(")"):
        depth = 0
        out_chars = []
        for ch in s:
            if ch == "(":
                depth += 1
            elif ch == ")":
                if depth == 0:
                    continue
                depth -= 1
            out_chars.append(ch)
        s = "".join(out_chars)
        # 仍有未闭合的左括号：从最后一个左括号处截断
        while "(" in s and s.count("(") != s.count(")"):
            s = s[: s.rfind("(")].rstrip(" ,;:-")
        s = " ".join(s.split())

    if len(s) < 2:
        s = "hypothetical protein"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--report", default=None)
    a = ap.parse_args()

    changed = []
    out = []
    biotype = {}
    for line in open(a.input, encoding="utf-8"):
        if line.startswith("#"):
            out.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            out.append(line)
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', f[8]))
        gid = attrs.get("gene_id", "")
        if f[2] == "gene":
            biotype[gid] = attrs.get("gene_biotype", "")
        # 没有 product 的 lncRNA：NCBI 会默认填 "mRNA"，显式补上 "lncRNA"
        if "product" not in attrs and f[2] in ("gene", "transcript"):
            bt = attrs.get("gene_biotype", biotype.get(gid, ""))
            if bt == "lncRNA":
                f[8] = f[8].rstrip().rstrip(";") + '; product "lncRNA";'
                line = "\t".join(f) + "\n"
        if "product" in attrs:
            old = attrs["product"]
            new = clean(old, attrs.get("gene_biotype", biotype.get(gid, "")))
            if new != old:
                changed.append((gid, f[2], old, new))
                f[8] = re.sub(r'product\s+"[^"]*"', 'product "%s"' % new.replace('"', ""), f[8], count=1)
                line = "\t".join(f) + "\n"
        out.append(line)

    with open(a.output, "w", encoding="utf-8") as fh:
        fh.writelines(out)

    to_hypo = sorted({g for g, _, _, n in changed if n == "hypothetical protein"})
    print(f"改写 {len(changed)} 条 product（涉及 {len({c[0] for c in changed})} 个基因）")
    print(f"其中变为 hypothetical protein 的基因 {len(to_hypo)} 个: {to_hypo[:10]}")

    if a.report:
        with open(a.report, "w", encoding="utf-8") as fh:
            fh.write("gene\tfeature\told\tnew\n")
            for g, feat, o, n in changed:
                fh.write(f"{g}\t{feat}\t{o}\t{n}\n")
        print("明细写入", a.report)


if __name__ == "__main__":
    sys.exit(main())
