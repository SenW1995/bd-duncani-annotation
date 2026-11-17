#!/usr/bin/env python3
"""汇总结构预测结果，生成可检索索引

产物（都在 struct_all/ 下）
  struct_index.tsv    每个结构一行：基因、产物、家族、长度、平均 pLDDT、高置信比例、文件路径
  struct_index.html   单文件可搜索表格（搜索框 + 列排序 + 置信度筛选 + 直接打开 PDB）
"""
from __future__ import annotations

import glob
import json
import os
import re

import pandas as pd

BASE = "/home/HeLan/project/bd_transcriptome"
OUT = f"{BASE}/struct_all"
OLD = f"{BASE}/struct"
GTF = f"{BASE}/v3/babesia_duncani_v3.2.gtf"
FAM = f"{BASE}/families_v32/families.tsv"
TARGETS = f"{OUT}/targets.tsv"


def biotypes():
    bio, prod = {}, {}
    for line in open(GTF):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
        bio[a["gene_id"]] = a.get("gene_biotype", "")
        prod[a["gene_id"]] = a.get("product", "")
    return bio, prod


def plddt_stats(pdb):
    """Boltz 把 pLDDT 写进 B-factor 列（0–100）"""
    vals, n = [], 0
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
                n += 1
            except ValueError:
                pass
    if not vals:
        return None, None, None, 0
    vals.sort()
    mean = sum(vals) / len(vals)
    hi = 100.0 * sum(1 for v in vals if v >= 70) / len(vals)
    vhi = 100.0 * sum(1 for v in vals if v >= 90) / len(vals)
    return mean, hi, vhi, n


def main():
    bio, prod = biotypes()
    fam = pd.read_csv(FAM, sep="\t")
    fam_idx = fam.drop_duplicates("gene_id").set_index("gene_id")
    tgt = pd.read_csv(TARGETS, sep="\t") if os.path.exists(TARGETS) else None
    tgt_len = dict(zip(tgt.gene_id, tgt.length)) if tgt is not None else {}

    rows = []
    for d, tag in ((OUT, "new"), (OLD, "earlier")):
        for p in sorted(glob.glob(f"{d}/*.pdb")):
            gid = os.path.basename(p)[:-4]
            mean, hi, vhi, nres = plddt_stats(p)
            size = int(fam_idx.family_size.get(gid, 1)) if gid in fam_idx.index else 1
            rows.append(dict(
                gene_id=gid,
                product=str(prod.get(gid, ""))[:100],
                family_id=str(fam_idx.family_id.get(gid, "")) if gid in fam_idx.index else "",
                family_size=size,
                representative="yes" if size >= 2 else "single-copy",
                length=tgt_len.get(gid) or nres,   # 先前折的结构不在 targets.tsv 里，用残基数代替
                mean_plddt=round(mean, 1) if mean else "",
                pct_ge70=round(hi, 1) if hi else "",
                pct_ge90=round(vhi, 1) if vhi else "",
                source=tag,
                pdb=os.path.relpath(p, OUT) if tag == "new" else os.path.relpath(p, OUT),
            ))

    df = pd.DataFrame(rows).sort_values(["family_size", "gene_id"], ascending=[False, True])
    df.to_csv(f"{OUT}/struct_index.tsv", sep="\t", index=False)

    data = df.fillna("").to_dict("records")
    html = """<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>B. duncani 结构预测索引</title>
<style>
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:18px;color:#1a1a1a}
h1{font-size:18px;margin:0 0 4px}
p.sub{color:#666;font-size:13px;margin:0 0 14px}
input,select{font-size:13px;padding:5px 8px;border:1px solid #ccc;border-radius:4px}
input[type=search]{width:320px}
table{border-collapse:collapse;margin-top:12px;font-size:12.5px;width:100%}
th,td{border-bottom:1px solid #e6e6e6;padding:5px 8px;text-align:left}
th{background:#fafafa;position:sticky;top:0;cursor:pointer;user-select:none;white-space:nowrap}
th:hover{background:#f0f0f0}
tr:hover td{background:#f7fbff}
td.n{text-align:right;font-variant-numeric:tabular-nums}
a{color:#0b62c4;text-decoration:none}
a:hover{text-decoration:underline}
.c{display:inline-block;min-width:44px;padding:1px 6px;border-radius:3px;color:#fff;font-size:11px}
</style></head><body>
<h1>Babesia duncani 全蛋白组结构（Boltz-2）</h1>
<p class="sub">共 <b id="n"></b> 个结构 · 多拷贝家族只保留一个代表 · 点击表头排序 · 输入基因号或产物名检索</p>
<div>
  <input type="search" id="q" placeholder="搜索基因号 / 产物 / 家族号…">
  <select id="fs">
    <option value="">全部家族规模</option>
    <option value="multi">仅多拷贝代表</option>
    <option value="single">仅单拷贝</option>
  </select>
  <select id="pl">
    <option value="0">全部置信度</option>
    <option value="70">平均 pLDDT ≥ 70</option>
    <option value="80">平均 pLDDT ≥ 80</option>
    <option value="90">平均 pLDDT ≥ 90</option>
  </select>
</div>
<table><thead><tr>
<th data-k="gene_id">gene_id</th><th data-k="product">product</th>
<th data-k="family_id">family</th><th data-k="family_size" class="n">copies</th>
<th data-k="length" class="n">aa</th><th data-k="mean_plddt" class="n">mean pLDDT</th>
<th data-k="pct_ge90" class="n">% ≥90</th><th>PDB</th>
</tr></thead><tbody id="tb"></tbody></table>
<script>
const D = __DATA__;
let sortK = "family_size", asc = false;
function col(v){ if(v>=90) return "#1a7f37"; if(v>=70) return "#4c9be8"; if(v>=50) return "#d9a300"; return "#c0392b"; }
function render(){
  const q = document.getElementById('q').value.toLowerCase();
  const fs = document.getElementById('fs').value;
  const pl = parseFloat(document.getElementById('pl').value);
  let rows = D.filter(r=>{
    if(q && !(r.gene_id+" "+r.product+" "+r.family_id).toLowerCase().includes(q)) return false;
    if(fs==="multi" && !(r.family_size>1)) return false;
    if(fs==="single" && !(r.family_size<=1)) return false;
    if(pl && !(parseFloat(r.mean_plddt)>=pl)) return false;
    return true;
  });
  rows.sort((a,b)=>{
    let x=a[sortK], y=b[sortK];
    const nx=parseFloat(x), ny=parseFloat(y);
    if(!isNaN(nx)&&!isNaN(ny)) return asc?nx-ny:ny-nx;
    return asc? String(x).localeCompare(String(y)) : String(y).localeCompare(String(x));
  });
  document.getElementById('n').textContent = rows.length;
  document.getElementById('tb').innerHTML = rows.map(r=>{
    const c = r.mean_plddt==="" ? "#999" : col(parseFloat(r.mean_plddt));
    return `<tr><td><a href="${r.pdb}">${r.gene_id}</a></td><td>${r.product}</td>
    <td>${r.family_id}</td><td class="n">${r.family_size}</td><td class="n">${r.length}</td>
    <td class="n"><span class="c" style="background:${c}">${r.mean_plddt}</span></td>
    <td class="n">${r.pct_ge90}</td><td><a href="${r.pdb}">open</a></td></tr>`;
  }).join("");
}
document.querySelectorAll('th[data-k]').forEach(th=>th.onclick=()=>{
  const k=th.dataset.k; asc = (sortK===k)? !asc : false; sortK=k; render();
});
document.getElementById('q').oninput=render;
document.getElementById('fs').onchange=render;
document.getElementById('pl').onchange=render;
render();
</script></body></html>"""
    html = html.replace("__DATA__", json.dumps(data, ensure_ascii=False))
    open(f"{OUT}/struct_index.html", "w", encoding="utf-8").write(html)

    print(f"结构 {len(df)} 个（new {int((df.source=='new').sum())} / earlier {int((df.source=='earlier').sum())}）")
    if len(df):
        print("平均 pLDDT 中位", round(pd.to_numeric(df.mean_plddt, errors="coerce").median(), 1))
        print("pLDDT>=70 的比例", f"{100*(pd.to_numeric(df.mean_plddt, errors='coerce')>=70).mean():.1f}%")
    print("写出 struct_index.tsv / struct_index.html")


if __name__ == "__main__":
    main()
