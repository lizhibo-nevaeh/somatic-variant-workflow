# CNVkit analysis using BQSR BAMs and optional Sequenza purity/ploidy estimates.

CNV_CFG = config.get("cnv", {})

CNV_PIPELINE_SCRIPT = CNV_CFG.get(
    "pipeline_script",
    f"{SCRIPT_DIR}/run_cnvkit.sh"
)
CNV_BQSR_DIR = CNV_CFG.get("bqsr_dir", BQSR_BAM_DIR)
CNV_SEQUENZA_SOURCE = str(
    CNV_CFG.get("sequenza_source", "workflow")
).lower()

if CNV_SEQUENZA_SOURCE not in {"external", "workflow"}:
    raise ValueError(
        "cnv.sequenza_source must be 'external' or 'workflow', "
        f"got: {CNV_SEQUENZA_SOURCE!r}"
    )

if CNV_SEQUENZA_SOURCE == "workflow":
    CNV_SEQUENZA_DIR = f"{OUT}/sequenza"
else:
    CNV_SEQUENZA_DIR = CNV_CFG.get("sequenza_dir", "")

CNV_TARGETS = CNV_CFG.get("targets", "")
CNV_ACCESS = CNV_CFG.get("access", "")
CNVKIT_BIN = CNV_CFG.get("cnvkit", "cnvkit.py")
CNV_THREADS = int(CNV_CFG.get("threads", 8))
CNV_MEM_MB = int(CNV_CFG.get("mem_mb", 50000))
CNV_RUNTIME = int(CNV_CFG.get("runtime", 1440))
CNV_QUEUE = str(CNV_CFG.get("queue", DEFAULT_PARTITION))


def cnv_normal_bam(wildcards):
    normal = NORMAL_BY_TUMOR[wildcards.tumor]
    return f"{CNV_BQSR_DIR}/{normal}.BQSR.bam"


def cnv_normal_bai(wildcards):
    return cnv_normal_bam(wildcards) + ".bai"


def cnv_tumor_bam(wildcards):
    return f"{CNV_BQSR_DIR}/{wildcards.tumor}.BQSR.bam"


def cnv_tumor_bai(wildcards):
    return cnv_tumor_bam(wildcards) + ".bai"


CNV_GENDER_BY_TUMOR = {}
with open(PAIR_FILE, encoding="utf-8") as _fh:
    for _line in _fh:
        _line = _line.strip()
        if not _line or _line.startswith("#"):
            continue
        _f = _line.split()
        if len(_f) >= 3:
            CNV_GENDER_BY_TUMOR[_f[1]] = _f[2].lower()


def cnv_purity(wildcards):
    return (
        f"{CNV_SEQUENZA_DIR}/"
        f"{wildcards.tumor}/"
        f"{wildcards.tumor}.purity_ploidy.txt"
    )


rule cnvkit_from_sequenza:
    input:
        normal_bam=cnv_normal_bam,
        normal_bai=cnv_normal_bai,
        tumor_bam=cnv_tumor_bam,
        tumor_bai=cnv_tumor_bai,
        purity=cnv_purity
    output:
        call=f"{OUT}/cnv/{{tumor}}/{{tumor}}.call.cns"
    params:
        normal=lambda wildcards: NORMAL_BY_TUMOR[wildcards.tumor],
        gender=lambda wildcards: CNV_GENDER_BY_TUMOR.get(
            wildcards.tumor, "unknown"
        )
    log:
        f"{OUT}/logs/cnv/{{tumor}}.cnv_pipeline.log"
    threads: CNV_THREADS
    resources:
        mem_mb=CNV_MEM_MB,
        runtime=CNV_RUNTIME,
        slurm_partition=CNV_QUEUE
    shell:
        r'''
        set -euo pipefail

        mkdir -p "$(dirname "{log}")"

        pair_file=$(mktemp)
        trap 'rm -f "$pair_file"' EXIT
        printf '%s\t%s\t%s\n' \
          "{params.normal}" "{wildcards.tumor}" "{params.gender}" > "$pair_file"

        BASE_OUT="{OUT}" \
        BQSR_DIR="{CNV_BQSR_DIR}" \
        BAM_SUFFIX=".BQSR.bam" \
        SEQUENZA_DIR="{CNV_SEQUENZA_DIR}" \
        PAIR_LIST="$pair_file" \
        ENV_FILE="{ENV_FILE}" \
        SEQ_TYPE="{SEQ_TYPE}" \
        CNV_TARGETS="{CNV_TARGETS}" \
        CNV_ACCESS="{CNV_ACCESS}" \
        CNVKIT_BIN="{CNVKIT_BIN}" \
        QUEUE="{CNV_QUEUE}" \
        CORE="{threads}" \
        MEM="50" \
        WALLTIME="24:00:00" \
        SUBMIT=0 \
        bash "{CNV_PIPELINE_SCRIPT}" "$pair_file" \
          > "{log}" 2>&1

        child="{OUT}/05_5_cnvkit/sh/{wildcards.tumor}.step5_5_cnvkit.sh"
        test -s "$child"
        bash "$child" >> "{log}" 2>&1
        test -s "{output.call}"
        '''
