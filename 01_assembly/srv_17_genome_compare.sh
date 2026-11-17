#!/bin/bash
# 基因组对照：我们的组装 vs 两篇已发表论文的组装
#   Singh 2023 (Nat Microbiol)  = GCA_028658345.1（我们的 v1_genome.fa，5 染色体 + 160 contig）
#   Wang 2022 (BMC Biol)        = GWHBECJ00000000（NECAT/Canu ONT，7.86 Mb，7 scaffold）
# 产出
#   ref/wang_genome.fna                  Wang 组装
#   compare/stats.tsv                    三套组装的统计
#   compare/{ours_vs_singh,ours_vs_wang}.paf   共线性比对（dot plot 用）
#   compare/synteny.tsv                  我们每条 chromosome 被两套组装覆盖的比例
set -eu
source /home/HeLan/miniconda/etc/profile.d/conda.sh
conda activate bioinfo

BASE=/home/HeLan/project/bd_transcriptome
REF=$BASE/ref
CMP=$BASE/compare
THREADS=32
mkdir -p "$CMP"

# ---- 1. Wang 组装 ----
if [ ! -s "$REF/wang_genome.fna" ]; then
    wget -q -O "$REF/wang_genome.fasta.gz" \
        "https://download.cncb.ac.cn/gwh/Protists/Babesia_duncani_Babesia_duncani_genome_GWHBECJ00000000/GWHBECJ00000000.genome.fasta.gz"
    zcat "$REF/wang_genome.fasta.gz" > "$REF/wang_genome.fna"
fi
echo "wang genome: $(grep -c '^>' "$REF/wang_genome.fna") seqs"

# ---- 2. 组装统计 ----
python - "$REF/babesia.fa" "$REF/v1_singh_genome.fna" "$REF/wang_genome.fna" "$CMP/stats.tsv" <<'PY'
import gzip, sys, re

def stats(path):
    lens, gc, total = [], 0, 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    s = "".join(cur); lens.append(len(s)); gc += s.count("G") + s.count("C") + s.count("g") + s.count("c")
                    total += len(s); cur = []
            else:
                cur.append(line.strip())
        if cur:
            s = "".join(cur); lens.append(len(s)); gc += s.count("G") + s.count("C")
            total += len(s)
    lens.sort(reverse=True)
    half, n50, l50 = total / 2, 0, 0
    c = 0
    for i, L in enumerate(lens, 1):
        c += L
        if c >= half:
            n50, l50 = L, i
            break
    return dict(path=path, n_seq=len(lens), total=total, gc=round(100 * gc / total, 2),
                max=lens[0], n50=n50, l50=l50, top5=sum(lens[:5]), n_ge1mb=sum(1 for L in lens if L >= 1_000_000))

rows = [stats(p) for p in sys.argv[1:4]]   # 3 个 fasta；argv[4] 是输出
names = ["ours_v2", "singh_v1", "wang_gwh"]
with open(sys.argv[4], "w") as fh:
    fh.write("assembly\tn_seq\ttotal_bp\tGC_pct\tmax_bp\tN50\tL50\ttop5_bp\tn_seq_ge1Mb\n")
    for name, r in zip(names, rows):
        fh.write(f"{name}\t{r['n_seq']}\t{r['total']}\t{r['gc']}\t{r['max']}\t{r['n50']}\t{r['l50']}\t{r['top5']}\t{r['n_ge1mb']}\n")
print(open(sys.argv[4]).read())
PY

# ---- 3. 共线性（asm5）----
for PAIR in "ours_vs_singh $REF/babesia.fa $REF/v1_singh_genome.fna" \
            "ours_vs_wang $REF/babesia.fa $REF/wang_genome.fna"; do
    set -- $PAIR
    TAG=$1; Q=$2; T=$3
    if [ ! -s "$CMP/$TAG.paf" ]; then
        minimap2 -x asm5 -t $THREADS "$T" "$Q" 2>/dev/null > "$CMP/$TAG.paf"
    fi
done

# ---- 4. 覆盖统计 ----
python - "$CMP" <<'PY'
import sys, pandas as pd
cmp_dir = sys.argv[1]
rows = []
for tag in ("ours_vs_singh", "ours_vs_wang"):
    paf = pd.read_csv(f"{cmp_dir}/{tag}.paf", sep="\t", usecols=[0, 1, 2, 3],
                      names=["q", "qlen", "qstart", "qend"], header=None, comment="#")
    cov = 0
    for q, g in paf.groupby("q"):
        iv = sorted(zip(g.qstart, g.qend))
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        cov += sum(e - s for s, e in merged)
        rows.append(dict(comparison=tag, query=q, qlen=int(g.qlen.iloc[0]),
                         covered_bp=sum(e - s for s, e in merged)))
d = pd.DataFrame(rows)
d["pct"] = (100 * d.covered_bp / d.qlen).round(2)
d.to_csv(f"{cmp_dir}/synteny.tsv", sep="\t", index=False)
print(d.to_string(index=False))
PY
echo DONE_COMPARE
