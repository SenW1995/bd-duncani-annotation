#!/usr/bin/env python3
"""全蛋白组结构预测（Boltz-2），多拷贝家族只取一个代表

设计
  * 目标集：v3.2 的 protein_coding 基因；多拷贝家族（家族成员 >=2）只保留一个代表，
    代表 = 该家族里最长的序列（结构覆盖度最好）。
  * 断点续跑：以输出目录里是否已有 <gene>.pdb 为准，重复运行自动跳过。
  * 分批：每批 N 个 yaml 调一次 boltz（模型只加载一次，批间可恢复），
    任何一批失败不会丢掉整个队列。
  * 分片：--shard/--nshard 把目标分到多张卡上，互不重叠。

用法
  python srv_22_fold_all.py prepare                 # 生成 targets.tsv / targets.faa 并打印统计
  python srv_22_fold_all.py run --shard 0 --nshard 4 --batch 40
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import sys
import time
import subprocess
import zlib
from collections import defaultdict

import pandas as pd

BASE = "/home/HeLan/project/bd_transcriptome"
GTF = f"{BASE}/v3/babesia_duncani_v3.2.gtf"
FAM = f"{BASE}/families_v32/families.tsv"
FAA = f"{BASE}/families_v32/all.faa"
OLD = f"{BASE}/struct"                      # 先前折的 61 个（药靶 + 人源同源）
OUT = f"{BASE}/struct_all"
WORK = f"{BASE}/work"
BOLTZ = "/home/HeLan/miniconda/envs/boltz/bin/boltz"
MAXLEN = 1200                               # Boltz-2 的实用上限


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


def done_ids():
    """已折过的基因：先前 struct/ 里的 + 本轮 struct_all/ 里的"""
    ids = set()
    for d in (OLD, OUT):
        for p in glob.glob(f"{d}/*.pdb"):
            ids.add(os.path.basename(p)[:-4])
    return ids


def prepare():
    bio, prod = biotypes()
    seqs = read_faa(FAA)
    fam = pd.read_csv(FAM, sep="\t")

    # 每个多拷贝家族取最长序列作代表
    rep_of = {}
    for fid, g in fam.groupby("family_id"):
        if len(g) < 2:
            continue
        cand = [(len(seqs.get(r.gene_id, "")), r.gene_id) for _, r in g.iterrows()]
        cand = [c for c in cand if c[0] > 0]
        if cand:
            rep_of[fid] = max(cand)[1]

    done = done_ids()
    rows, skipped_done, skipped_len = [], 0, 0
    for _, r in fam.iterrows():
        gid = r.gene_id
        if bio.get(gid) != "protein_coding":
            continue
        size = int(r.family_size)
        if size >= 2 and rep_of.get(r.family_id, gid) != gid:
            continue                              # 非代表成员，跳过
        L = len(seqs.get(gid, ""))
        if L == 0:
            continue
        if gid in done:
            skipped_done += 1
            continue
        if L > MAXLEN:
            skipped_len += 1
            continue
        rows.append(dict(gene_id=gid, family_id=r.family_id, family_size=size,
                         length=L, product=str(prod.get(gid, ""))[:120]))

    t = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    t.to_csv(f"{OUT}/targets.tsv", sep="\t", index=False)
    with open(f"{OUT}/targets.faa", "w") as fh:
        for gid in t.gene_id:
            s = seqs[gid]
            fh.write(f">{gid}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    print(f"protein_coding 总数      {int((pd.Series(bio) == 'protein_coding').sum())}")
    print(f"多拷贝家族代表保留      {len(rep_of)}")
    print(f"已在 struct/ 中，跳过   {skipped_done}")
    print(f"超过 {MAXLEN} aa，跳过     {skipped_len}")
    print(f"本次待折                {len(t)}")
    print(f"其中多拷贝代表          {int((t.family_size >= 2).sum())}")
    print(f"长度中位 / 最大         {int(t.length.median())} / {int(t.length.max())}")


def gpu_occupied_by_others(gpu):
    """返回占用说明（别人的进程 / 显存不足），空闲返回 None。"""
    def sh(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=30).stdout
        except Exception:
            return ""
    me = os.environ.get("USER", "")
    apps = sh(["nvidia-smi", "-i", str(gpu), "--query-compute-apps=pid,used_memory",
               "--format=csv,noheader"])
    foreign_mb = 0
    foreign = []
    for line in apps.strip().splitlines():
        parts = [x.strip() for x in line.split(",")]
        if len(parts) < 2 or not parts[0].isdigit():
            continue
        pid, mb = parts[0], int(parts[1].split()[0]) if parts[1].split()[0].isdigit() else 0
        user = sh(["ps", "-o", "user=", "-p", pid]).strip()
        if user and user != me:
            foreign_mb += mb
            foreign.append(f"{user}:{pid}({mb}MB)")
    # 小占用（<2 GB，多半是监控/预览进程）不视为冲突；真跑起来的任务才让路
    if foreign_mb >= 2000:
        return f"GPU{gpu} 上有别人的任务 {' '.join(foreign)}"
    free = sh(["nvidia-smi", "-i", str(gpu), "--query-gpu=memory.free",
               "--format=csv,noheader,nounits"]).strip()
    if free.isdigit() and int(free) < 8000:
        return f"GPU{gpu} 只剩 {int(free)/1024:.1f} GB 显存"
    return None


def run(shard, nshard, batch, guard_gpu=None, guard_exit_minutes=15):
    t = pd.read_csv(f"{OUT}/targets.tsv", sep="\t")
    seqs = read_faa(f"{OUT}/targets.faa")
    # 用 crc32 而不是内置 hash：后者每个进程随机加盐，分片会不一致
    todo = [g for g in t.gene_id if (zlib.crc32(g.encode()) % nshard) == shard]
    os.makedirs(OUT, exist_ok=True)
    tag = f"{os.getpid()}"
    indir = f"{WORK}/boltz_all_in_{tag}"
    tmp_out = f"{WORK}/boltz_all_out_{tag}"
    os.makedirs(indir, exist_ok=True)

    total = len(todo)
    print(f"shard {shard}/{nshard}: {total} 个目标", flush=True)
    waited = 0
    while True:
        pending = [g for g in todo if not os.path.exists(f"{OUT}/{g}.pdb")]
        if not pending:
            print(f"shard {shard} 完成", flush=True)
            return
        if guard_gpu is not None:
            busy = gpu_occupied_by_others(guard_gpu)
            if busy:
                waited += 1
                print(f"shard {shard}: {busy}，让路等待（第 {waited} 次）", flush=True)
                if waited * 5 >= guard_exit_minutes:
                    print(f"shard {shard}: 让路超时，退出释放显存，稍后由外层脚本重启", flush=True)
                    sys.exit(3)
                time.sleep(300)
                continue
            waited = 0
        chunk = pending[:batch]
        # 清掉上一批残留的 yaml，避免重复计算
        for f in glob.glob(f"{indir}/*.yaml"):
            os.remove(f)
        for gid in chunk:
            with open(f"{indir}/{gid}.yaml", "w") as fh:
                fh.write("version: 1\nsequences:\n")
                fh.write(f"  - protein:\n      id: A\n      sequence: {seqs[gid]}\n")
        cmd = [BOLTZ, "predict", indir, "--out_dir", tmp_out,
               "--output_format", "pdb", "--model", "boltz2",
               "--devices", "1", "--accelerator", "gpu", "--no_kernels",
               "--use_msa_server"]
        subprocess.run(cmd, check=False)
        got = 0
        for gid in chunk:
            if os.path.exists(f"{OUT}/{gid}.pdb"):
                got += 1
                continue
            hits = sorted(glob.glob(f"{tmp_out}/**/{gid}/*model_0.pdb", recursive=True))
            if not hits:
                hits = sorted(glob.glob(f"{tmp_out}/**/{gid}*.pdb", recursive=True))
            if hits:
                shutil.copy(hits[0], f"{OUT}/{gid}.pdb")
                got += 1
        left = sum(1 for g in todo if not os.path.exists(f"{OUT}/{g}.pdb"))
        print(f"shard {shard}: 本批 {got}/{len(chunk)}，剩余 {left}/{total}", flush=True)
        if got == 0:
            # 整批为 0：可能是显存不足或 MSA 服务抖动。逐条重试，把真正失败的登记后跳过，
            # 否则整条队列会卡在同一批上无限循环。
            print(f"shard {shard}: 整批失败，逐条重试", flush=True)
            failed = []
            for gid in chunk:
                for f in glob.glob(f"{indir}/*.yaml"):
                    os.remove(f)
                with open(f"{indir}/{gid}.yaml", "w") as fh:
                    fh.write("version: 1\nsequences://n")
                    fh.write(f"  - protein://n      id: A\n      sequence: {seqs[gid]}\n")
                subprocess.run([BOLTZ, "predict", indir, "--out_dir", tmp_out,
                                "--output_format", "pdb", "--model", "boltz2",
                                "--devices", "1", "--accelerator", "gpu", "--no_kernels",
                                "--use_msa_server"], check=False)
                hits = sorted(glob.glob(f"{tmp_out}/**/{gid}/*model_0.pdb", recursive=True))
                if not hits:
                    hits = sorted(glob.glob(f"{tmp_out}/**/{gid}*.pdb", recursive=True))
                if hits:
                    shutil.copy(hits[0], f"{OUT}/{gid}.pdb")
                else:
                    failed.append(gid)
                    with open(f"{OUT}/failed.tsv", "a") as fh:
                        fh.write(f"{gid}\t{len(seqs[gid])}\n")
            print(f"shard {shard}: 逐条重试完成，失败 {len(failed)} 条", flush=True)
            continue


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "run"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--batch", type=int, default=40)
    ap.add_argument("--guard-gpu", type=int, default=None,
                    help="指定卡号后，每批开始前检查该卡是否被别人占用，占用就让路")
    a = ap.parse_args()
    if a.mode == "prepare":
        prepare()
    else:
        run(a.shard, a.nshard, a.batch, a.guard_gpu)
