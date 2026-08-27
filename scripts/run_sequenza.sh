#!/usr/bin/env bash
# Step5_sequenza（Slurm版，可外部传参）：Sequenza 纯度/倍性估计
#   1) gc_wiggle（缺失时自动生成）
#   2) bam2seqz（按染色体）+ 合并（含 SIGPIPE 修复）
#   3) seqz_binning 降采样
#   4) sequenza R 分析 → purity / ploidy
# Output: ${BASE_OUT}/sequenza/<TUMOR>/<TUMOR>.purity_ploidy.txt
#
# 用法：
#   A) 默认跑 PAIR_LIST：
#      bash run_sequenza.sh
#   B) 传入一个 pair_list 文件（两列：normal \t tumor）：
#      bash run_sequenza.sh /path/to/pairs.tsv
#   C) 传入一个/多个 tumor 名（从 PAIR_LIST 第二列筛选补跑）：
#      bash run_sequenza.sh TUMOR001 TUMOR002
#
# Environment-variable override example:
#   BASE_OUT=... PAIR_LIST=... QUEUE=compute CORE=8 MEM=64 WALLTIME=24:00:00 \
#   bash run_sequenza.sh
#
# WGS/WES：Sequenza 对两者通用，无需区分；WGS 的 seqz 文件较大，可调大 BIN_SIZE。
#
# 控制是否提交：SUBMIT=0 只生成子脚本；SUBMIT=1 生成并 sbatch 提交（默认 1）

set -euo pipefail
[ -n "${BASH_VERSION:-}" ] || { echo "请用: bash $0"; exit 2; }

###########################
# 配置参数（允许 env 覆盖）
###########################
ENV_FILE="${ENV_FILE:-}"

# Output root
BASE_OUT="${BASE_OUT:-./results}"

# BQSR BAM directory
BQSR_DIR="${BQSR_DIR:-${BASE_OUT}/BQSR}"
BAM_SUFFIX="${BAM_SUFFIX:-.BQSR.bam}"

# 配对表（两列：normal tumor；可有第3列 gender）
DEFAULT_PAIR_LIST="${DEFAULT_PAIR_LIST:-./config/pairs.example.tsv}"
PAIR_LIST="${PAIR_LIST:-$DEFAULT_PAIR_LIST}"
DEFAULT_GENDER="${DEFAULT_GENDER:-male}"

# Reference resources
REF_ROOT="${REF_ROOT:-/path/to/reference}"
REF_FASTA="${REF_FASTA:-${REF_ROOT}/GRCh38.fasta}"
GC_WIGGLE="${GC_WIGGLE:-${REF_ROOT}/GRCh38.gc50.wig.gz}"
BIN_SIZE="${BIN_SIZE:-50}"

# Sequenza 工具
SEQUENZA_BIN_DIR="${SEQUENZA_BIN_DIR:-/path/to/sequenza/bin}"
SEQUENZA_UTILS="${SEQUENZA_UTILS:-${SEQUENZA_BIN_DIR}/sequenza-utils}"
RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript || true)}"

# 资源
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-8}"
MEM="${MEM:-64}"
WALLTIME="${WALLTIME:-24:00:00}"

# 是否提交
SUBMIT="${SUBMIT:-1}"

###########################
# 入参：pair_list 文件 / tumor 列表
###########################
if [[ $# -gt 0 ]]; then
  if [[ -f "$1" ]]; then
    PAIR_LIST="$1"
  else
    _tmp_pair="$(mktemp)"
    trap 'rm -f "$_tmp_pair"' EXIT
    awk -v OFS="\t" '
      BEGIN{ for(i=1;i<ARGC;i++){ want[ARGV[i]]=1 } }
      NF>=2 && want[$2]{ print $0 }
    ' "$@" < <(awk 'NF>=2{print $0}' "$PAIR_LIST" | tr -d '\r' | sed '1s/^\xEF\xBB\xBF//') > "$_tmp_pair"
    PAIR_LIST="$_tmp_pair"
  fi
fi

###########################
# 检查
###########################
[[ -z "$ENV_FILE" || -s "$ENV_FILE" ]] || { echo "[FATAL] environment file not found: $ENV_FILE" >&2; exit 1; }
[[ -s "$PAIR_LIST" ]]  || { echo "[FATAL] 缺少/为空 PAIR_LIST: $PAIR_LIST" >&2; exit 1; }
[[ -s "$SEQUENZA_UTILS" ]] || { echo "[FATAL] 找不到 sequenza-utils: $SEQUENZA_UTILS" >&2; exit 1; }
[[ -s "$RSCRIPT_BIN" ]] || { echo "[FATAL] 找不到 Rscript: $RSCRIPT_BIN" >&2; exit 1; }
[[ -s "$REF_FASTA" ]]  || { echo "[FATAL] 缺少参考 fasta: $REF_FASTA" >&2; exit 1; }
[[ -d "$BQSR_DIR" ]]   || { echo "[FATAL] 缺少 BQSR 目录: $BQSR_DIR" >&2; exit 1; }

# 输出目录
mkdir -p "${BASE_OUT}/05_sequenza/sh" \
         "${BASE_OUT}/sequenza" \
         "${BASE_OUT}/sequenza_seqz" \
         "${BASE_OUT}/sequenza_logs"

# 配对表是否带 gender（第3列）
PAIR_COLS=$(awk 'NF>=2{print NF; exit}' <(tr -d '\r' < "$PAIR_LIST"))
if [[ "${PAIR_COLS:-2}" -ge 3 ]]; then
  HAS_GENDER=true
else
  HAS_GENDER=false
fi

clean_stream(){ awk 'NF>=2{print $0}' "$PAIR_LIST" | tr -d '\r' | sed '1s/^\xEF\xBB\xBF//'; }

total=$(clean_stream | wc -l | awk '{print $1}')
echo "== Step5 (Slurm): Sequenza 纯度/倍性估计 =="
echo "BASE_OUT     : $BASE_OUT"
echo "BQSR_DIR     : $BQSR_DIR"
echo "PAIR_LIST    : $PAIR_LIST"
echo "REF_FASTA    : $REF_FASTA"
echo "GC_WIGGLE    : $GC_WIGGLE"
echo "HAS_GENDER   : $HAS_GENDER (否则用默认 $DEFAULT_GENDER)"
echo "QUEUE/CORE/MEM/WALLTIME : $QUEUE / $CORE / ${MEM}G / $WALLTIME"
echo "配对总数     : $total"
echo ""

###########################
# 生成子脚本
###########################
gen=0
while IFS=$'\t ' read -r normal tumor gender _rest; do
  [[ -z "${normal:-}" || -z "${tumor:-}" ]] && continue
  [[ "$HAS_GENDER" == true && -n "${gender:-}" ]] || gender="$DEFAULT_GENDER"

  tumor_bam="${BQSR_DIR}/${tumor}${BAM_SUFFIX}"
  normal_bam="${BQSR_DIR}/${normal}${BAM_SUFFIX}"
  if [[ ! -s "$tumor_bam" ]]; then
    echo "[WARN] 找不到 tumor BAM，跳过：$tumor_bam" >&2
    continue
  fi
  if [[ ! -s "$normal_bam" ]]; then
    echo "[WARN] 找不到 normal BAM，跳过：$normal_bam" >&2
    continue
  fi

  shp="${BASE_OUT}/05_sequenza/sh/${tumor}.step5_sequenza.sh"
  cat > "$shp" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J step5_seq___TUMOR__
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/05_sequenza/sh/__TUMOR__.step5_sequenza.out
#SBATCH -e __BASE__/05_sequenza/sh/__TUMOR__.step5_sequenza.err

set -uo pipefail

[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__" || true
export PATH="__SEQUENZA_BIN_DIR__:${PATH}"

SEQUENZA_UTILS="__SEQUENZA_UTILS__"
RSCRIPT="__RSCRIPT_BIN__"
REF_FASTA="__REF_FASTA__"
GC_WIGGLE="__GC_WIGGLE__"
BIN_SIZE="__BIN_SIZE__"

TUMOR="__TUMOR__"
NORMAL="__NORMAL__"
GENDER="__GENDER__"

tumor_bam="__BQSR_DIR__/__TUMOR____BAM_SUFFIX__"
normal_bam="__BQSR_DIR__/__NORMAL____BAM_SUFFIX__"

seqz_dir="__BASE__/sequenza_seqz"
out_dir="__BASE__/sequenza/__TUMOR__"
log_dir="__BASE__/sequenza_logs"
mkdir -p "$seqz_dir" "$out_dir" "$log_dir"

seqz_file="${seqz_dir}/${TUMOR}.seqz.gz"
seqz_binned="${seqz_dir}/${TUMOR}.small.seqz.gz"
purity_file="${out_dir}/${TUMOR}.purity_ploidy.txt"

[[ -s "$tumor_bam" ]]  || { echo "缺少 tumor BAM: $tumor_bam"; exit 2; }
[[ -s "$normal_bam" ]] || { echo "缺少 normal BAM: $normal_bam"; exit 2; }

echo "[$(date)] === Step5 开始: $TUMOR (normal=$NORMAL, gender=$GENDER) ==="
pipeline_start=$(date +%s)

# ========= Step 0: GC wiggle =========
GC_WIGGLE_USE="${GC_WIGGLE%.wig.gz}.tab.wig.gz"
if [[ ! -s "$GC_WIGGLE" ]]; then
  echo "[$(date)] [0] 生成 GC wiggle"
  wig_tmp="${GC_WIGGLE%.gz}"
  rm -f "$wig_tmp" "$GC_WIGGLE"
  "$SEQUENZA_UTILS" gc_wiggle -w "$BIN_SIZE" -f "$REF_FASTA" -o "$wig_tmp" \
    || { echo "gc_wiggle 生成失败"; exit 3; }
  gzip -f "$wig_tmp"
fi
if [[ ! -s "$GC_WIGGLE_USE" ]]; then
  zcat "$GC_WIGGLE" | awk 'BEGIN{OFS="\t"}
      /^variableStep/ {print; next}
      NF==0 {next}
      NF>=2 {print $1,$2; next}
      {print}' | gzip -c > "$GC_WIGGLE_USE"
fi
[[ -s "$GC_WIGGLE_USE" ]] || { echo "GC wiggle 不存在: $GC_WIGGLE_USE"; exit 3; }

# ========= Step 1: bam2seqz（按染色体）+ 合并 =========
echo "[$(date)] [1/3] bam2seqz + 合并"
step1_start=$(date +%s)

chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
if [[ "$GENDER" == "male" ]]; then
  chromosomes="$chromosomes chrX chrY"
else
  chromosomes="$chromosomes chrX"
fi

RUN_BAM2SEQZ=true
if [[ -s "$seqz_file" ]]; then
  seqz_size_mb=$(( $(stat -c%s "$seqz_file" 2>/dev/null || echo 0) / 1048576 ))
  if (( seqz_size_mb >= 10 )); then
    set +o pipefail
    line_count=$(zcat "$seqz_file" 2>/dev/null | head -n 2000 | wc -l 2>/dev/null || echo 0)
    set -o pipefail
    if (( line_count >= 1000 )); then
      echo "  ✓ 已有有效 seqz (${seqz_size_mb}MB)，跳过 bam2seqz"
      RUN_BAM2SEQZ=false
    fi
  fi
  [[ "$RUN_BAM2SEQZ" == true ]] && rm -f "$seqz_file"
fi

if [[ "$RUN_BAM2SEQZ" == true ]]; then
  rm -f ${seqz_dir}/${TUMOR}.chr*.seqz.gz
  failed_chrs=()
  for chr in $chromosomes; do
    chr_seqz="${seqz_dir}/${TUMOR}.${chr}.seqz.gz"
    echo "  处理 $chr ..."
    if ! "$SEQUENZA_UTILS" bam2seqz \
        -n "$normal_bam" -t "$tumor_bam" \
        --fasta "$REF_FASTA" -gc "$GC_WIGGLE_USE" \
        -C "$chr" -o "$chr_seqz" \
        > "${log_dir}/${TUMOR}.bam2seqz.${chr}.log" 2>&1; then
      echo "  [ERROR] $chr 失败"
      failed_chrs+=("$chr")
      rm -f "$chr_seqz"
    fi
  done
  if (( ${#failed_chrs[@]} > 5 )); then
    echo "[FATAL] 过多染色体失败 (${#failed_chrs[@]})"; exit 4
  fi

  # 合并（修复 SIGPIPE：用 || true 忽略管道错误，靠输出文件判断成败）
  echo "  合并染色体文件..."
  shopt -s nullglob
  chr_files=( "${seqz_dir}/${TUMOR}".chr*.seqz.gz )
  shopt -u nullglob
  (( ${#chr_files[@]} > 0 )) || { echo "[FATAL] 无 chr seqz 文件"; exit 4; }

  tmp_seqz="${seqz_dir}/${TUMOR}.merged.tmp.$$.seqz"
  rm -f "$tmp_seqz" "$seqz_file"
  zcat "${chr_files[0]}" 2>/dev/null | head -n 1 > "$tmp_seqz" || true
  [[ -s "$tmp_seqz" ]] || { echo "[FATAL] 表头提取失败"; rm -f "$tmp_seqz"; exit 4; }
  grep -q "chromosome" "$tmp_seqz" || { echo "[FATAL] 表头无效"; rm -f "$tmp_seqz"; exit 4; }

  for f in "${chr_files[@]}"; do
    [[ -s "$f" ]] && { zcat "$f" 2>/dev/null | tail -n +2 >> "$tmp_seqz" || true; }
  done
  tmp_lines=$(wc -l < "$tmp_seqz")
  (( tmp_lines >= 1000 )) || { echo "[FATAL] 合并后行数过少: $tmp_lines"; rm -f "$tmp_seqz"; exit 4; }

  gzip -c "$tmp_seqz" > "$seqz_file" || { echo "[FATAL] gzip 失败"; rm -f "$tmp_seqz" "$seqz_file"; exit 4; }
  rm -f "$tmp_seqz" ${seqz_dir}/${TUMOR}.chr*.seqz.gz
fi
[[ -s "$seqz_file" ]] || { echo "[FATAL] seqz 生成失败"; exit 4; }
step1_time=$(( $(date +%s) - step1_start ))
echo "  ✓ bam2seqz 完成，耗时 ${step1_time}s"

# ========= Step 2: seqz_binning =========
echo "[$(date)] [2/3] seqz_binning"
step2_start=$(date +%s)
if [[ ! -s "$seqz_binned" ]]; then
  "$SEQUENZA_UTILS" seqz_binning -s "$seqz_file" -w "$BIN_SIZE" -o "$seqz_binned" \
      > "${log_dir}/${TUMOR}.seqz_binning.log" 2>&1 \
    || { echo "[FATAL] seqz_binning 失败"; exit 5; }
fi
[[ -s "$seqz_binned" ]] || { echo "[FATAL] binned 文件缺失"; exit 5; }
step2_time=$(( $(date +%s) - step2_start ))
echo "  ✓ seqz_binning 完成，耗时 ${step2_time}s"

# ========= Step 3: Sequenza R 分析 =========
echo "[$(date)] [3/3] Sequenza R 分析"
step3_start=$(date +%s)
if [[ ! -s "$purity_file" ]]; then
  if [[ "$GENDER" == "female" ]]; then is_female="TRUE"; else is_female="FALSE"; fi
  R_script="${out_dir}/${TUMOR}_sequenza.R"
  cat > "$R_script" <<RSCRIPT
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(sequenza))
seqz_file <- "$seqz_binned"
sample_id <- "$TUMOR"
out_dir   <- "$out_dir"
cat("读取 seqz 文件...\n")
seqz_data <- sequenza.extract(seqz_file, verbose = FALSE)
cat("拟合模型...\n")
CP <- sequenza.fit(seqz_data, female = $is_female)
cat("提取结果...\n")
sequenza.results(sequenza.extract = seqz_data, cp.table = CP,
                 sample.id = sample_id, out.dir = out_dir)
confints_file <- file.path(out_dir, paste0(sample_id, "_confints_CP.txt"))
if (file.exists(confints_file)) {
    confints <- read.table(confints_file, header = TRUE, sep = "\t")
    purity_ploidy <- data.frame(purity = confints\$cellularity[1],
                                ploidy = confints\$ploidy.estimate[1])
    write.table(purity_ploidy,
                file.path(out_dir, paste0(sample_id, ".purity_ploidy.txt")),
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat("purity:", confints\$cellularity[1], " ploidy:", confints\$ploidy.estimate[1], "\n")
}
cat("Sequenza 分析完成!\n")
RSCRIPT
  "$RSCRIPT" "$R_script" > "${log_dir}/${TUMOR}.sequenza_R.log" 2>&1 \
    || echo "[WARN] R 分析返回非零，检查 ${log_dir}/${TUMOR}.sequenza_R.log"
fi
step3_time=$(( $(date +%s) - step3_start ))

confints_file="${out_dir}/${TUMOR}_confints_CP.txt"
if [[ ! -s "$purity_file" && -s "$confints_file" ]]; then
  purity=$(tail -n1 "$confints_file" | cut -f1)
  ploidy=$(tail -n1 "$confints_file" | cut -f2)
  printf 'purity\tploidy\n%s\t%s\n' "$purity" "$ploidy" > "$purity_file"
fi
[[ -s "$purity_file" ]] || { echo "[FATAL] 未生成 purity/ploidy 结果"; exit 6; }

purity=$(tail -n1 "$purity_file" | cut -f1)
ploidy=$(tail -n1 "$purity_file" | cut -f2)
total_time=$(( $(date +%s) - pipeline_start ))
echo "  ✓ R 分析完成，耗时 ${step3_time}s"
echo ""
echo "  purity=$purity  ploidy=$ploidy"
echo "  耗时: bam2seqz ${step1_time}s | binning ${step2_time}s | R ${step3_time}s | 总 ${total_time}s"
echo "[$(date)] === Step5 done: $TUMOR ==="
echo "  → purity/ploidy (deliverable): $purity_file"
EOSH

  sed -i \
    -e "s#__QUEUE__#${QUEUE}#g" \
    -e "s#__CORE__#${CORE}#g" \
    -e "s#__MEM__#${MEM}#g" \
    -e "s#__WALLTIME__#${WALLTIME}#g" \
    -e "s#__TUMOR__#${tumor}#g" \
    -e "s#__NORMAL__#${normal}#g" \
    -e "s#__GENDER__#${gender}#g" \
    -e "s#__BASE__#${BASE_OUT}#g" \
    -e "s#__ENV_FILE__#${ENV_FILE}#g" \
    -e "s#__BQSR_DIR__#${BQSR_DIR}#g" \
    -e "s#__BAM_SUFFIX__#${BAM_SUFFIX}#g" \
    -e "s#__SEQUENZA_UTILS__#${SEQUENZA_UTILS}#g" \
    -e "s#__SEQUENZA_BIN_DIR__#${SEQUENZA_BIN_DIR}#g" \
    -e "s#__RSCRIPT_BIN__#${RSCRIPT_BIN}#g" \
    -e "s#__REF_FASTA__#${REF_FASTA}#g" \
    -e "s#__GC_WIGGLE__#${GC_WIGGLE}#g" \
    -e "s#__BIN_SIZE__#${BIN_SIZE}#g" \
    "$shp"

  chmod +x "$shp"
  echo "→ 生成 step5 子脚本: $shp"
  gen=$((gen+1))
done < <(clean_stream)

echo ""
echo "== 生成完成: $gen / $total =="

###########################
# 批量提交
###########################
if [[ "$SUBMIT" == "1" ]]; then
  echo "== 开始批量提交 Step5 (sbatch) =="
  sub=0
  for s in "${BASE_OUT}/05_sequenza/sh/"*.step5_sequenza.sh; do
    [[ -s "$s" ]] || continue
    jid=$(sbatch "$s" | awk '{print $4}') || { echo "[WARN] sbatch 失败 -> $s"; continue; }
    echo "提交: $s -> $jid"
    sub=$((sub+1))
  done
  echo "Step5 提交完成: $sub / $gen"
else
  echo "SUBMIT=0: 只生成不提交（手动 sbatch ${BASE_OUT}/05_sequenza/sh/*.sh）"
fi
