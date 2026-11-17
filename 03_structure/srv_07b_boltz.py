#!/usr/bin/env python3
"""用 Boltz-2 预测蛋白结构，输出 struct/<id>.pdb

为什么不用 ESMFold：服务器上的 esmfold 环境只有 torch + esm，缺 openfold/omegaconf，
补装要编译 C++ 扩展；而 boltz 环境已装好 boltz2_conf + boltz2_aff 权重并能直接跑，
精度也更好（带 MSA）。

用法
  python srv_07b_boltz.py <faa> <outdir> [--max-len 1200] [--shard 0 --nshard 2]
  --shard/--nshard 用于把任务分到两张卡上并行。

产物
  <outdir>/<id>.pdb        每个输入序列的最佳模型
  <outdir>/boltz_in/       每个序列的 yaml（可复现）
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys

BOLTZ = "/home/HeLan/miniconda/envs/boltz/bin/boltz"
WORK = "/home/HeLan/project/bd_transcriptome/work"


def read_faa(path):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("faa")
    ap.add_argument("outdir")
    ap.add_argument("--max-len", type=int, default=1200)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--msa", action="store_true", help="使用 MMseqs2 在线 MSA 服务")
    a = ap.parse_args()

    seqs = read_faa(a.faa)
    todo = [(k, v) for k, v in seqs.items() if len(v) <= a.max_len]
    todo = [x for i, x in enumerate(todo) if i % a.nshard == a.shard]
    print(f"{a.faa}: {len(seqs)} 条，本分片 {len(todo)} 条", flush=True)
    if not todo:
        return

    tag = os.path.basename(a.faa).replace(".faa", "") + f"_{a.shard}"
    indir = f"{WORK}/boltz_in_{tag}"
    outdir = f"{WORK}/boltz_out_{tag}"
    os.makedirs(indir, exist_ok=True)
    os.makedirs(a.outdir, exist_ok=True)

    for k, v in todo:
        yml = f"{indir}/{k}.yaml"
        if os.path.exists(yml):
            continue
        with open(yml, "w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {v}\n")

    cmd = [BOLTZ, "predict", indir, "--out_dir", outdir,
           "--output_format", "pdb", "--model", "boltz2",
           "--devices", "1", "--accelerator", "gpu", "--no_kernels"]
    if a.msa:
        cmd.append("--use_msa_server")
    print("RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=False)

    n = 0
    for k, _ in todo:
        hits = sorted(glob.glob(f"{outdir}/**/{k}/*model_0.pdb", recursive=True))
        if not hits:
            hits = sorted(glob.glob(f"{outdir}/**/{k}*.pdb", recursive=True))
        if hits:
            shutil.copy(hits[0], f"{a.outdir}/{k}.pdb")
            n += 1
    print(f"复制 {n}/{len(todo)} 个结构到 {a.outdir}", flush=True)


if __name__ == "__main__":
    main()
