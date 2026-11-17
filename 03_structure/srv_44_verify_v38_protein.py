#!/usr/bin/env python3
"""用 v3.8 注释现场翻译蛋白，与结构 PDB 序列比对（权威一致性判定）

比对对象
  PDB（按 CA 残基号取序列） vs 由 v3.8 GTF 的 CDS 按相/链拼接翻译出的蛋白
输出
  tsv：gene_id、pdb_id、状态（exact/substitutions/length_mismatch/no_pdb）、
       pdb_len、cds_aa_len、identity
"""
import os
import re
import sys
import csv
import collections

CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def translate(cds):
    return "".join(CODON.get(cds[i:i + 3].upper(), "X") for i in range(0, len(cds) - 2, 3))


def main():
    gtf, fa, pdbmap, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    genome = read_fasta(fa)
    cds = collections.defaultdict(list)
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        m = re.search(r'gene_id "([^"]+)"', f[8])
        if m:
            phase = int(f[7]) if f[7].isdigit() else 0
            cds[m.group(1)].append((f[0], int(f[3]), int(f[4]), f[6], phase))
    # pdbmap: gene_id \t pdb_id \t pdb_path
    rows = []
    for line in open(pdbmap):
        g, pid, path = line.rstrip("\n").split("\t")
        segs = sorted(cds.get(g, []), key=lambda x: x[1])   # 一律升序；负链整段反向互补
        if not segs:
            rows.append(dict(gene_id=g, pdb_id=pid, status="no_cds", pdb_len=0, aa_len=0, identity=""))
            continue
        chrom = segs[0][0]
        seq = genome.get(chrom, "")
        if not seq:
            rows.append(dict(gene_id=g, pdb_id=pid, status="no_chrom", pdb_len=0, aa_len=0, identity=""))
            continue
        nt = "".join(seq[s - 1:e] for _, s, e, _, _ in segs)
        if segs[0][3] == "-":
            nt = nt.translate(COMP)[::-1]
        prot = translate(nt)
        if prot.endswith("*"):
            prot = prot[:-1]
        pdb_seq = []
        seen = set()
        for ln in open(path):
            if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                k = ln[22:26].strip()
                if k in seen:
                    continue
                seen.add(k)
                pdb_seq.append(AA3.get(ln[17:20].strip().upper(), "X"))
        pdb_seq = "".join(pdb_seq)
        if pdb_seq == prot:
            st = "exact"
        elif len(pdb_seq) == len(prot):
            st = "substitutions:%d" % sum(1 for a, b in zip(pdb_seq, prot) if a != b)
        else:
            st = "length_mismatch"
        same = sum(1 for a, b in zip(pdb_seq, prot) if a == b)
        rows.append(dict(gene_id=g, pdb_id=pid, status=st, pdb_len=len(pdb_seq),
                         aa_len=len(prot), identity=f"{same / max(1, len(prot)) * 100:.1f}"))
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    c = collections.Counter(r["status"].split(":")[0] for r in rows)
    print("核对", len(rows), "个结构:", dict(c))
    for r in rows:
        if r["status"] != "exact":
            print("  ", r)


if __name__ == "__main__":
    main()
