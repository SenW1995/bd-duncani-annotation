#!/usr/bin/env python3
"""Foldseek 网页 API 批量搜索（未知功能蛋白的结构邻居）

用法
  python3 srv_35_foldseek_api.py --manifest foldseek_manifest.tsv --pdb-dir struct_all \
      --out foldseek_results.tsv [--workers 3] [--database afdb50] [--limit N]

要点
  * 每个蛋白一次提交（POST /api/ticket），轮询 /api/ticket/<id>，完成后取
    /api/result/<id>/0。单条约 10–20 s，3 并发约 5–7 s/条。
  * 提交节流：全局最小间隔 1 s，避免把公共网页服务打爆。
  * 断点续跑：结果表里已有的 pdb_id 直接跳过；每条写完立即 flush。
  * 失败重试 3 次（指数退避），仍失败记 status=failed，不阻塞队列。
  * 只保留前 5 个命中，字段：target（含注释）、prob、evalue、seqId、覆盖度。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import threading
import time
import urllib.error
import urllib.request
import uuid

API = "https://search.foldseek.com/api"
_tls = threading.local()
_submit_lock = threading.Lock()
_last_submit = [0.0]


def _post_file(path, fields, timeout=180):
    boundary = "----wb" + uuid.uuid4().hex
    body = bytearray()
    for k, v in fields.items():
        body += (f"--{boundary}\r\nContent-Disposition: form-data; "
                 f'name="{k}"\r\n\r\n{v}\r\n').encode()
    name = os.path.basename(path)
    body += (f"--{boundary}\r\nContent-Disposition: form-data; "
             f'name="q"; filename="{name}"\r\n'
             "Content-Type: application/octet-stream\r\n\r\n").encode()
    with open(path, "rb") as fh:
        body += fh.read()
    body += f"\r\n--{boundary}--\r\n".encode()
    req = urllib.request.Request(
        API + "/ticket", data=bytes(body),
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _submit(path, databases, mode, timeout=180):
    boundary = "----wb" + uuid.uuid4().hex
    body = bytearray()
    for db in databases:
        body += (f"--{boundary}\r\nContent-Disposition: form-data; "
                 f'name="database[]"\r\n\r\n{db}\r\n').encode()
    body += (f"--{boundary}\r\nContent-Disposition: form-data; "
             f'name="mode"\r\n\r\n{mode}\r\n').encode()
    name = os.path.basename(path)
    body += (f"--{boundary}\r\nContent-Disposition: form-data; "
             f'name="q"; filename="{name}"\r\n'
             "Content-Type: application/octet-stream\r\n\r\n").encode()
    with open(path, "rb") as fh:
        body += fh.read()
    body += f"\r\n--{boundary}--\r\n".encode()
    req = urllib.request.Request(
        API + "/ticket", data=bytes(body),
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _get_json(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.loads(r.read().decode())


def search_one(pdb_path, databases, mode, poll=15, max_wait=600):
    attempt = 0
    while True:
        try:
            with _submit_lock:                       # 提交节流：全局 ≥1 s 一次
                gap = time.time() - _last_submit[0]
                if gap < 5.0:
                    time.sleep(5.0 - gap)
                _last_submit[0] = time.time()
                tid = _submit(pdb_path, databases, mode)["id"]
            t0 = time.time()
            while time.time() - t0 < max_wait:
                st = _get_json(f"{API}/ticket/{tid}").get("status")
                if st == "COMPLETE":
                    return _get_json(f"{API}/result/{tid}/0")
                if st == "ERROR":
                    raise RuntimeError("ticket ERROR")
                time.sleep(poll)
            raise RuntimeError("timeout")
        except Exception as exc:                      # noqa: BLE001
            # 429 是公共服务的配额限制，不是错误：长时间退避后继续，不消耗重试次数
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 429:
                print("  配额 429，等 300 s 再试", flush=True)
                time.sleep(300)
                continue
            attempt += 1
            if attempt > 4:
                raise
            wait = min(90, 10 * 2 ** attempt)
            print(f"  重试 {attempt + 1}/4：{type(exc).__name__} "
                  f"{getattr(exc, 'code', '')}，等 {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def flatten(res):
    """把每个库的命中拉平成一条列表，按概率降序；带 db 标记。"""
    out = []
    for block in res.get("results", []) or []:
        db = block.get("db", "")
        try:
            hits = block["alignments"][0]
        except (KeyError, IndexError, TypeError):
            continue
        for h in hits:
            tgt = str(h.get("target", ""))
            acc, _, desc = tgt.partition(" ")
            qlen, dlen = h.get("qLen") or 1, h.get("dbLen") or 1
            denom = min(qlen, dlen)
            cov = h.get("alnLength", 0) / denom * 100 if denom else 0
            out.append(dict(db=db, target=acc, desc=desc, prob=h.get("prob"),
                            evalue=h.get("eval"), seqid=h.get("seqId"),
                            cov=round(cov, 1), aln=h.get("alnLength")))
    out.sort(key=lambda x: (-(x["prob"] or 0), x["evalue"] if x["evalue"] is not None else 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--database", action="append", default=None,
                    help="可重复给多个库，默认 afdb50+pdb100")
    ap.add_argument("--mode", default="3diaa")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if not a.database:
        a.database = ["afdb50", "pdb100"]

    rows = list(csv.DictReader(open(a.manifest, encoding="utf-8"), delimiter="\t"))
    done = set()
    if os.path.exists(a.out):
        done = {r["pdb_id"] for r in csv.DictReader(open(a.out, encoding="utf-8"), delimiter="\t")}
    todo = [r for r in rows if r["pdb_id"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"清单 {len(rows)} 条，已完成 {len(done)}，本次处理 {len(todo)}", flush=True)

    cols = (["pdb_id", "gene_id", "length_aa", "current_product", "family_id", "family_members",
             "status", "n_hits"]
            + [f"top{i}_{k}" for i in (1, 2, 3, 5)
               for k in ("db", "target", "desc", "prob", "evalue", "seqid", "cov")])
    new = not os.path.exists(a.out)
    fh = open(a.out, "a", encoding="utf-8", newline="")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    if new:
        w.writeheader()
        fh.flush()
    lock = threading.Lock()
    counter = {"ok": 0, "fail": 0}

    def work(r):
        p = os.path.join(a.pdb_dir, r["pdb_id"] + ".pdb")
        rec = {k: r.get(k, "") for k in ("pdb_id", "gene_id", "length_aa",
                                         "current_product", "family_id", "family_members")}
        try:
            hits = flatten(search_one(p, a.database, a.mode))
            rec["status"] = "ok" if hits else "no_hit"
            rec["n_hits"] = len(hits)
            for i in (1, 2, 3):
                if len(hits) >= i:
                    h = hits[i - 1]
                    rec[f"top{i}_db"] = h["db"]
                    rec[f"top{i}_target"], rec[f"top{i}_desc"] = h["target"], h["desc"]
                    rec[f"top{i}_prob"], rec[f"top{i}_evalue"] = h["prob"], h["evalue"]
                    rec[f"top{i}_seqid"], rec[f"top{i}_cov"] = h["seqid"], h["cov"]
            if len(hits) >= 5:
                h = hits[4]
                rec["top5_db"] = h["db"]
                rec["top5_target"], rec["top5_desc"] = h["target"], h["desc"]
                rec["top5_prob"], rec["top5_evalue"] = h["prob"], h["evalue"]
                rec["top5_seqid"], rec["top5_cov"] = h["seqid"], h["cov"]
        except Exception as exc:                      # noqa: BLE001
            code = getattr(exc, "code", "")
            rec["status"] = f"failed:{type(exc).__name__}{code}"
            rec["n_hits"] = 0
        with lock:
            w.writerow(rec)
            fh.flush()
            counter["ok" if rec["status"] in ("ok", "no_hit") else "fail"] += 1
            n = counter["ok"] + counter["fail"]
            if n % 10 == 0:
                print(f"  进度 {n}/{len(todo)}（失败 {counter['fail']}）", flush=True)

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(work, todo))
    fh.close()
    print(f"完成：成功 {counter['ok']}，失败 {counter['fail']}", flush=True)


if __name__ == "__main__":
    main()
