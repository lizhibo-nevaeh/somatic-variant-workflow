# 01_7_normal_contam_qc.smk
# Normal contamination QC (report-only)
#
# Reuses the same GATK GetPileupSummaries parameters/resources as
# 02_somatic_calling.smk, but evaluates normal samples independently.
#
# This module:
#   1) runs GetPileupSummaries on each normal BQSR BAM
#   2) runs CalculateContamination on each normal pileup table
#   3) summarizes all normals into qc_report/normal_contamination.tsv
#
# It DOES NOT:
#   - remove/exclude normals
#   - modify PON construction
#   - impose a contamination cutoff

NORMAL_CONTAM_JAVA_MEM_G = int(SOMATIC.get("java_mem_g", 80))

NORMAL_CONTAM_ROOT = f"{OUT}/qc/normal_contamination"
NORMAL_CONTAM_SUMMARY = f"{OUT}/qc_report/normal_contamination.tsv"


rule normal_contam_pileup:
    input:
        bam=lambda wc: bqsr_bam(wc.normal),
        bai=lambda wc: bqsr_bai(wc.normal),
    output:
        table=f"{NORMAL_CONTAM_ROOT}/pileup/{{normal}}.pileups.table"
    log:
        f"{OUT}/logs/qc/normal_contamination/{{normal}}.pileup.log"
    threads: PILEUP_THREADS
    resources:
        mem_mb=PILEUP_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.table})" "$(dirname {log})"
        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

        "{GATK}" \
          --java-options "-Xms{NORMAL_CONTAM_JAVA_MEM_G}G -Xmx{NORMAL_CONTAM_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
          GetPileupSummaries \
          -I "{input.bam}" \
          -L "{INTERVALS}" \
          -V "{GNOMAD_AF}" \
          -O "{output.table}" \
          > "{log}" 2>&1

        test -s "{output.table}"
        """


rule normal_contam_calculate:
    input:
        pileup=f"{NORMAL_CONTAM_ROOT}/pileup/{{normal}}.pileups.table"
    output:
        table=f"{NORMAL_CONTAM_ROOT}/table/{{normal}}.contamination.table"
    log:
        f"{OUT}/logs/qc/normal_contamination/{{normal}}.calculate.log"
    threads: CONTAM_THREADS
    resources:
        mem_mb=CONTAM_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.table})" "$(dirname {log})"
        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

        "{GATK}" \
          --java-options "-Xms{NORMAL_CONTAM_JAVA_MEM_G}G -Xmx{NORMAL_CONTAM_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
          CalculateContamination \
          -I "{input.pileup}" \
          -O "{output.table}" \
          > "{log}" 2>&1

        test -s "{output.table}"
        """


rule normal_contam_summary:
    input:
        tables=expand(
            f"{NORMAL_CONTAM_ROOT}/table/{{normal}}.contamination.table",
            normal=NORMALS
        )
    output:
        summary=NORMAL_CONTAM_SUMMARY
    log:
        f"{OUT}/logs/qc/normal_contamination/summary.log"
    run:
        from pathlib import Path

        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)

        rows = []

        for normal, table_path in zip(NORMALS, input.tables):
            data_row = None

            with open(table_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if line.lower().startswith("sample"):
                        continue
                    data_row = line
                    break

            if data_row is None:
                raise ValueError(
                    f"No contamination result row found for {normal}: {table_path}"
                )

            fields = data_row.split()
            if len(fields) < 3:
                raise ValueError(
                    f"Unexpected CalculateContamination output for {normal}: {data_row}"
                )

            gatk_sample = fields[0]
            contamination = float(fields[1])
            error = float(fields[2])

            rows.append(
                (normal, gatk_sample, contamination, contamination * 100.0, error)
            )

        with open(output.summary, "w", encoding="utf-8") as out_fh:
            out_fh.write(
                "normal_id\tgatk_sample\tcontamination\t"
                "contamination_percent\terror\n"
            )
            for normal, gatk_sample, contamination, contamination_pct, error in rows:
                out_fh.write(
                    f"{normal}\t{gatk_sample}\t{contamination:.10g}\t"
                    f"{contamination_pct:.6f}\t{error:.10g}\n"
                )

        with open(log[0], "w", encoding="utf-8") as log_fh:
            log_fh.write(
                f"normal_contam_summary: normals={len(rows)} "
                f"output={output.summary}\n"
            )
