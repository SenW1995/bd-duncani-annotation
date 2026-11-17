#!/bin/bash
# 四套组装/注释的 BUSCO 对比：genome 模式 ×4、protein 模式 ×4
set -u
cd /home/HeLan/project/bd_transcriptome || exit 1
source /home/HeLan/miniconda/etc/profile.d/conda.sh
conda activate annotation
mkdir -p busco_cmp
export BUSCO_CONFIG_FILE=""

run_genome () {   # $1 = 标签, $2 = fasta
  busco -i "$2" -m genome -l apicomplexa_odb10 -c 12 -o "g_$1" \
        --download_path work/busco_dl --force --offline \
        --out_path busco_cmp > "busco_cmp/g_$1.log" 2>&1
}

run_prot () {     # $1 = 标签, $2 = faa
  busco -i "$2" -m protein -l apicomplexa_odb10 -c 12 -o "p_$1" \
        --download_path work/busco_dl --force --offline \
        --out_path busco_cmp > "busco_cmp/p_$1.log" 2>&1
}

echo "=== start $(date +%H:%M) ==="
run_genome ours   V2/babesia_duncani_no_apicoplast.fa &
run_genome singh  ref/singh_genome.fna &
wait
echo "=== genome 1/2 done $(date +%H:%M) ==="

run_genome singhd ref/singh_dedup.fna &
run_genome wang   ref/wang_genome.fna &
wait
echo "=== genome done $(date +%H:%M) ==="

run_prot ours   work/v39_proteins.faa &
run_prot singh  ref/singh_protein.faa &
wait
echo "=== protein 1/2 done $(date +%H:%M) ==="

run_prot singhd ref/singh_dedup_protein.faa &
[ -s ref/wang_protein.faa ] && run_prot wang ref/wang_protein.faa
wait
echo "=== all done $(date +%H:%M) ==="

for d in busco_cmp/*/short_summary*.txt; do
  [ -f "$d" ] || continue
  echo "--- $d"
  grep -m1 "C:" "$d"
done
