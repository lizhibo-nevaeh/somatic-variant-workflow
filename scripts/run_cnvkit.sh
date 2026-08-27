#!/usr/bin/env bash
# CNVkit driver for paired tumor-normal BAMs.
# Supports WES (targets BED) and WGS (--method wgs).
# If Sequenza purity/ploidy results are available, they are used by cnvkit call.

set -euo pipefail
[ -n "${BASH_VERSION:-}" ] || { echo "Please run with bash: $0"; exit 2; }

ENV_FILE="${ENV_FILE:-}"
BASE_OUT="${BASE_OUT:-./results}"
BQSR_DIR="${BQSR_DIR:-${BASE_OUT}/BQSR}"
BAM_SUFFIX="${BAM_SUFFIX:-.BQSR.bam}"
SEQUENZA_DIR="${SEQUENZA_DIR:-${BASE_OUT}/sequenza}"
DEFAULT_PAIR_LIST="${DEFAULT_PAIR_LIST:-./config/pairs.example.tsv}"
PAIR_LIST="${PAIR_LIST:-$DEFAULT_PAIR_LIST}"
DEFAULT_GENDER="${DEFAULT_GENDER:-unknown}"
REF_ROOT="${REF_ROOT:-/path/to/reference}"
REF_FASTA="${REF_FASTA:-${REF_ROOT}/GRCh38.fasta}"
SEQ_TYPE="${SEQ_TYPE:-WES}"
CNV_TARGETS="${CNV_TARGETS:-${REF_ROOT}/targets.bed}"
CNV_ACCESS="${CNV_ACCESS:-}"
CNVKIT_BIN="${CNVKIT_BIN:-cnvkit.py}"
DEFAULT_PLOIDY="${DEFAULT_PLOIDY:-2}"
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-8}"
MEM="${MEM:-50}"
WALLTIME="${WALLTIME:-24:00:00}"
SUBMIT="${SUBMIT:-1}"

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

if [[ "$SEQ_TYPE" == "WGS" ]]; then
  CNV_TARGETS=""
elif [[ -z "$CNV_TARGETS" || ! -s "$CNV_TARGETS" ]]; then
  echo "[WARN] WES targets were not found; falling back to CNVkit WGS mode." >&2
  CNV_TARGETS=""
fi

[[ -z "$ENV_FILE" || -s "$ENV_FILE" ]] || { echo "[FATAL] environment file not found: $ENV_FILE" >&2; exit 1; }
[[ -s "$PAIR_LIST" ]] || { echo "[FATAL] pair list not found or empty: $PAIR_LIST" >&2; exit 1; }
command -v "$CNVKIT_BIN" >/dev/null 2>&1 || [[ -x "$CNVKIT_BIN" ]] || { echo "[FATAL] CNVkit not found: $CNVKIT_BIN" >&2; exit 1; }
[[ -s "$REF_FASTA" ]] || { echo "[FATAL] reference FASTA not found: $REF_FASTA" >&2; exit 1; }
[[ -d "$BQSR_DIR" ]] || { echo "[FATAL] BQSR directory not found: $BQSR_DIR" >&2; exit 1; }

mkdir -p "${BASE_OUT}/05_5_cnvkit/sh" "${BASE_OUT}/cnv" "${BASE_OUT}/cnv_logs"

PAIR_COLS=$(awk 'NF>=2{print NF; exit}' <(tr -d '\r' < "$PAIR_LIST"))
if [[ "${PAIR_COLS:-2}" -ge 3 ]]; then HAS_GENDER=true; else HAS_GENDER=false; fi
clean_stream(){ awk 'NF>=2{print $0}' "$PAIR_LIST" | tr -d '\r' | sed '1s/^\xEF\xBB\xBF//'; }

total=$(clean_stream | wc -l | awk '{print $1}')
echo "CNVkit pairs: $total"

gen=0
while IFS=$'\t ' read -r normal tumor gender _rest; do
  [[ -z "${normal:-}" || -z "${tumor:-}" ]] && continue
  [[ "$HAS_GENDER" == true && -n "${gender:-}" ]] || gender="$DEFAULT_GENDER"

  if [[ "$gender" == "male" ]]; then sample_sex="male"; elif [[ "$gender" == "female" ]]; then sample_sex="female"; else sample_sex=""; fi

  tumor_bam="${BQSR_DIR}/${tumor}${BAM_SUFFIX}"
  normal_bam="${BQSR_DIR}/${normal}${BAM_SUFFIX}"
  [[ -s "$tumor_bam" ]] || { echo "[WARN] tumor BAM not found: $tumor_bam" >&2; continue; }
  [[ -s "$normal_bam" ]] || { echo "[WARN] normal BAM not found: $normal_bam" >&2; continue; }

  shp="${BASE_OUT}/05_5_cnvkit/sh/${tumor}.step5_5_cnvkit.sh"
  cat > "$shp" <<'EOSH'
#!/usr/bin/env bash
#SBATCH -p __QUEUE__
#SBATCH -J cnvkit___TUMOR__
#SBATCH -c __CORE__
#SBATCH --mem=__MEM__G
#SBATCH -t __WALLTIME__
#SBATCH -o __BASE__/05_5_cnvkit/sh/__TUMOR__.cnvkit.out
#SBATCH -e __BASE__/05_5_cnvkit/sh/__TUMOR__.cnvkit.err

set -euo pipefail
[[ -n "__ENV_FILE__" && -s "__ENV_FILE__" ]] && source "__ENV_FILE__" || true

CNVKIT="__CNVKIT_BIN__"
REF_FASTA="__REF_FASTA__"
CNV_TARGETS="__CNV_TARGETS__"
CNV_ACCESS="__CNV_ACCESS__"
TUMOR="__TUMOR__"
NORMAL="__NORMAL__"
GENDER="__GENDER__"
SAMPLE_SEX="__SAMPLE_SEX__"
DEFAULT_PLOIDY="__DEFAULT_PLOIDY__"
tumor_bam="__BQSR_DIR__/__TUMOR____BAM_SUFFIX__"
normal_bam="__BQSR_DIR__/__NORMAL____BAM_SUFFIX__"
out_dir="__BASE__/cnv/__TUMOR__"
log_dir="__BASE__/cnv_logs"
mkdir -p "$out_dir" "$log_dir"

cnr_file="${out_dir}/${TUMOR}.cnr"
cns_file="${out_dir}/${TUMOR}.cns"
call_cns="${out_dir}/${TUMOR}.call.cns"

purity=""
ploidy="$DEFAULT_PLOIDY"
purity_file="__SEQUENZA_DIR__/__TUMOR__/__TUMOR__.purity_ploidy.txt"
confints_file="__SEQUENZA_DIR__/__TUMOR__/__TUMOR___confints_CP.txt"
if [[ -s "$purity_file" ]]; then
  purity=$(tail -n1 "$purity_file" | cut -f1)
  ploidy=$(tail -n1 "$purity_file" | cut -f2)
elif [[ -s "$confints_file" ]]; then
  purity=$(tail -n1 "$confints_file" | cut -f1)
  ploidy=$(tail -n1 "$confints_file" | cut -f2)
fi
case "$purity" in ""|NA|purity|cellularity) purity="" ;; esac
case "$ploidy" in ""|NA|ploidy) ploidy="$DEFAULT_PLOIDY" ;; esac

export MPLCONFIGDIR="$out_dir"
batch_cmd=("$CNVKIT" batch "$tumor_bam" --normal "$normal_bam" --fasta "$REF_FASTA" --output-dir "$out_dir" -p __CORE__)
if [[ -n "$CNV_TARGETS" && -s "$CNV_TARGETS" ]]; then
  batch_cmd+=(--targets "$CNV_TARGETS")
else
  batch_cmd+=(--method wgs)
fi
[[ -n "$CNV_ACCESS" && -s "$CNV_ACCESS" ]] && batch_cmd+=(--access "$CNV_ACCESS")
"${batch_cmd[@]}" > "${log_dir}/${TUMOR}.cnvkit_batch.log" 2>&1

bam_base=$(basename "$tumor_bam" .bam)
[[ -s "${out_dir}/${bam_base}.cnr" ]] && mv "${out_dir}/${bam_base}.cnr" "$cnr_file"
[[ -s "${out_dir}/${bam_base}.cns" ]] && mv "${out_dir}/${bam_base}.cns" "$cns_file"
[[ -s "$cnr_file" && -s "$cns_file" ]] || { echo "[FATAL] expected CNVkit batch outputs were not generated"; exit 3; }

call_cmd=("$CNVKIT" call "$cns_file" -m clonal --ploidy "$ploidy" -o "$call_cns")
[[ -n "$SAMPLE_SEX" ]] && call_cmd+=(--sample-sex "$SAMPLE_SEX")
[[ -n "$purity" ]] && call_cmd+=(--purity "$purity")
"${call_cmd[@]}" > "${log_dir}/${TUMOR}.cnvkit_call.log" 2>&1
[[ -s "$call_cns" ]] || { echo "[FATAL] CNVkit call output was not generated"; exit 4; }

if [[ -s "$cnr_file" && -s "$cns_file" ]]; then
  "$CNVKIT" scatter "$cnr_file" -s "$cns_file" -o "${out_dir}/${TUMOR}.scatter.png" >/dev/null 2>&1 || true
fi
"$CNVKIT" diagram "$cns_file" -o "${out_dir}/${TUMOR}.diagram.pdf" >/dev/null 2>&1 || true
EOSH

  sed -i \
    -e "s#__QUEUE__#${QUEUE}#g" \
    -e "s#__CORE__#${CORE}#g" \
    -e "s#__MEM__#${MEM}#g" \
    -e "s#__WALLTIME__#${WALLTIME}#g" \
    -e "s#__TUMOR__#${tumor}#g" \
    -e "s#__NORMAL__#${normal}#g" \
    -e "s#__GENDER__#${gender}#g" \
    -e "s#__SAMPLE_SEX__#${sample_sex}#g" \
    -e "s#__BASE__#${BASE_OUT}#g" \
    -e "s#__ENV_FILE__#${ENV_FILE}#g" \
    -e "s#__BQSR_DIR__#${BQSR_DIR}#g" \
    -e "s#__BAM_SUFFIX__#${BAM_SUFFIX}#g" \
    -e "s#__SEQUENZA_DIR__#${SEQUENZA_DIR}#g" \
    -e "s#__CNVKIT_BIN__#${CNVKIT_BIN}#g" \
    -e "s#__DEFAULT_PLOIDY__#${DEFAULT_PLOIDY}#g" \
    -e "s#__REF_FASTA__#${REF_FASTA}#g" \
    -e "s#__CNV_TARGETS__#${CNV_TARGETS}#g" \
    -e "s#__CNV_ACCESS__#${CNV_ACCESS}#g" \
    "$shp"
  chmod +x "$shp"
  gen=$((gen+1))
done < <(clean_stream)

echo "Generated CNVkit jobs: $gen / $total"

if [[ "$SUBMIT" == "1" ]]; then
  for s in "${BASE_OUT}/05_5_cnvkit/sh/"*.step5_5_cnvkit.sh; do
    [[ -s "$s" ]] || continue
    sbatch "$s"
  done
else
  echo "SUBMIT=0: job scripts generated but not submitted."
fi
