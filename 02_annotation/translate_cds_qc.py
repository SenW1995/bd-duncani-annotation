#!/usr/bin/env python3
"""Translate a CDS fasta (NCBI cds_from_genomic style) and report integrity stats.

Outputs:
  - protein fasta (id = protein_id, or cds index if absent)
  - TSV: cds_id, protein_id, gene, length_aa, internal_stops, starts_with_M, ends_with_stop
"""
import argparse
import re
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

HDR = re.compile(r"\[protein=([^\]]*)\]")
PID = re.compile(r"\[protein_id=([^\]]*)\]")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cds_fasta")
    ap.add_argument("out_faa")
    ap.add_argument("out_tsv")
    ap.add_argument("--table", type=int, default=1, help="NCBI translation table")
    args = ap.parse_args()

    rows = []
    n = 0
    with open(args.out_faa, "w") as fh:
        for rec in SeqIO.parse(args.cds_fasta, "fasta"):
            n += 1
            desc = rec.description
            m_gene = HDR.search(desc)
            m_pid = PID.search(desc)
            gene = m_gene.group(1) if m_gene else "NA"
            pid = m_pid.group(1) if m_pid else f"cds_{n}"
            nt = str(rec.seq).upper()
            nt = nt[: len(nt) - len(nt) % 3]
            try:
                aa = str(Seq(nt).translate(table=args.table))
            except Exception as exc:  # noqa: BLE001
                print(f"WARN translate failed {pid}: {exc}", file=sys.stderr)
                continue
            core = aa[:-1] if aa.endswith("*") else aa
            rows.append(
                dict(
                    cds_id=rec.id,
                    protein_id=pid,
                    gene=gene,
                    nt_len=len(nt),
                    length_aa=len(core),
                    internal_stops=core.count("*"),
                    starts_with_M=int(core.startswith("M")),
                    ends_with_stop=int(aa.endswith("*")),
                )
            )
            fh.write(f">{pid} {gene}\n{core}\n")

    pd.DataFrame(rows).to_csv(args.out_tsv, sep="\t", index=False)
    print(f"translated {len(rows)} CDS -> {args.out_faa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
