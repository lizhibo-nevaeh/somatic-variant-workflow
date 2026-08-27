# BAM-level QC with samtools flagstat and samtools stats.

QC_CFG = config.get("qc", {})
QC_BAM_DIR = str(QC_CFG.get("bam_dir", BQSR_BAM_DIR))
QC_BAM_SUFFIX = str(QC_CFG.get("bam_suffix", ".BQSR.bam"))
QC_THREADS = int(QC_CFG.get("threads", 4))
QC_MEM_MB = int(QC_CFG.get("mem_mb", 8000))
QC_RUNTIME = int(QC_CFG.get("runtime", 120))
QC_QUEUE = str(QC_CFG.get("queue", DEFAULT_PARTITION))


def qc_bam(wildcards):
    return f"{QC_BAM_DIR}/{wildcards.sample}{QC_BAM_SUFFIX}"


rule bam_qc_collect:
    input:
        bam=qc_bam
    output:
        flagstat=f"{OUT}/qc/flagstat/{{sample}}.flagstat",
        stats=f"{OUT}/qc/stats/{{sample}}.stats"
    log:
        f"{OUT}/logs/qc/{{sample}}.bam_qc.log"
    threads: QC_THREADS
    resources:
        mem_mb=QC_MEM_MB,
        runtime=QC_RUNTIME,
        slurm_partition=QC_QUEUE
    shell:
        r'''
        set -euo pipefail

        mkdir -p \
          "$(dirname "{output.flagstat}")" \
          "$(dirname "{output.stats}")" \
          "$(dirname "{log}")"

        [[ -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

        SAMTOOLS=$(command -v samtools || true)
        if [[ -z "$SAMTOOLS" ]]; then
            echo "[FATAL] samtools not found" > "{log}"
            exit 1
        fi

        "$SAMTOOLS" flagstat -@ {threads} "{input.bam}" \
          > "{output.flagstat}" 2> "{log}"

        "$SAMTOOLS" stats -@ {threads} "{input.bam}" \
          > "{output.stats}" 2>> "{log}"

        test -s "{output.flagstat}"
        test -s "{output.stats}"
        '''


rule bam_qc_complete:
    input:
        flagstat=expand(
            f"{OUT}/qc/flagstat/{{sample}}.flagstat",
            sample=SAMPLES
        ),
        stats=expand(
            f"{OUT}/qc/stats/{{sample}}.stats",
            sample=SAMPLES
        )
