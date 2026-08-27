# 01_6_pair_qc.smk
# Tumor/Normal identity QC with Somalier
#
# Assumes workflow/rules/00_common.smk has already defined:
#   OUT, TUMORS, NORMAL_BY_TUMOR, BQSR_BAM_DIR, REFS
#
# Design:
#   1) one Somalier extract job per BAM
#   2) one Somalier relate job per expected tumor-normal pair
#   3) one lightweight cohort-level pair_qc.tsv summary
#
# This module does NOT modify PON logic and does NOT calculate normal contamination.

PAIR_QC_CFG = config.get("pair_qc", {})

SOMALIER_BIN = str(PAIR_QC_CFG.get("somalier", "somalier"))
SOMALIER_SITES = str(PAIR_QC_CFG.get("sites", REFS.get("somalier_sites", "")))
SOMALIER_REF = str(PAIR_QC_CFG.get("reference", REF_FASTA))

PAIR_REL_PASS = float(PAIR_QC_CFG.get("relatedness_pass", 0.90))
PAIR_CONC_PASS = float(PAIR_QC_CFG.get("concordance_pass", 0.90))

PAIR_QC_THREADS = int(PAIR_QC_CFG.get("threads", 2))
PAIR_QC_MEM_MB = int(PAIR_QC_CFG.get("mem_mb", 8000))
PAIR_QC_RUNTIME = int(PAIR_QC_CFG.get("runtime", 60))
PAIR_QC_QUEUE = str(PAIR_QC_CFG.get("queue", DEFAULT_PARTITION))


rule somalier_extract:
    input:
        bam=lambda w: f"{BQSR_BAM_DIR}/{w.sample}.BQSR.bam",
        bai=lambda w: f"{BQSR_BAM_DIR}/{w.sample}.BQSR.bam.bai",
        sites=SOMALIER_SITES,
        ref=SOMALIER_REF,
    output:
        f"{OUT}/qc/somalier/{{sample}}.somalier"
    log:
        f"{OUT}/logs/qc/somalier_extract/{{sample}}.log"
    threads:
        PAIR_QC_THREADS
    resources:
        mem_mb=PAIR_QC_MEM_MB,
        runtime=PAIR_QC_RUNTIME,
        slurm_partition=PAIR_QC_QUEUE
    shell:
        r'''
        set -euo pipefail

        mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"

        "{SOMALIER_BIN}" extract \
            -d "$(dirname "{output}")" \
            --sites "{input.sites}" \
            -f "{input.ref}" \
            "{input.bam}" \
            > "{log}" 2>&1

        test -s "{output}"
        '''


rule somalier_pair_relate:
    input:
        normal=lambda w: (
            f"{OUT}/qc/somalier/{NORMAL_BY_TUMOR[w.tumor]}.somalier"
        ),
        tumor=lambda w: f"{OUT}/qc/somalier/{w.tumor}.somalier",
    output:
        pairs=f"{OUT}/qc/pair_qc/raw/{{tumor}}.pairs.tsv",
        samples=f"{OUT}/qc/pair_qc/raw/{{tumor}}.samples.tsv",
        html=f"{OUT}/qc/pair_qc/raw/{{tumor}}.html",
    log:
        f"{OUT}/logs/qc/somalier_relate/{{tumor}}.log"
    params:
        prefix=lambda w: f"{OUT}/qc/pair_qc/raw/{w.tumor}."
    threads:
        1
    resources:
        mem_mb=4000,
        runtime=30,
        slurm_partition=PAIR_QC_QUEUE
    shell:
        r'''
        set -euo pipefail

        mkdir -p "$(dirname "{output.pairs}")" "$(dirname "{log}")"

        "{SOMALIER_BIN}" relate \
            --output-prefix "{params.prefix}" \
            "{input.normal}" \
            "{input.tumor}" \
            > "{log}" 2>&1

        test -s "{output.pairs}"
        test -s "{output.samples}"
        test -s "{output.html}"
        '''


rule pair_qc_summary:
    input:
        expand(
            f"{OUT}/qc/pair_qc/raw/{{tumor}}.pairs.tsv",
            tumor=TUMORS,
        )
    output:
        f"{OUT}/qc_report/pair_qc.tsv"
    log:
        f"{OUT}/logs/qc/pair_qc_summary.log"
    run:
        import csv
        import os

        os.makedirs(os.path.dirname(str(output[0])), exist_ok=True)
        os.makedirs(os.path.dirname(str(log[0])), exist_ok=True)

        header = [
            "normal_id",
            "tumor_id",
            "relatedness",
            "ibs0",
            "ibs2",
            "concordance",
            "pair_status",
            "pair_note",
        ]

        rows = []

        for tumor, pair_file in zip(TUMORS, input):
            normal = NORMAL_BY_TUMOR[tumor]

            with open(pair_file, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rec = next(reader, None)

            if rec is None:
                raise ValueError(f"Empty Somalier pairs.tsv: {pair_file}")

            sample_a = rec.get("#sample_a", rec.get("sample_a", ""))
            sample_b = rec.get("sample_b", "")

            try:
                relatedness = float(rec["relatedness"])
                concordance = float(rec["concordance"])
                ibs0 = int(float(rec["ibs0"]))
                ibs2 = int(float(rec["ibs2"]))
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"Unexpected Somalier pairs.tsv format: {pair_file}"
                ) from exc

            notes = []

            if {sample_a, sample_b} != {normal, tumor}:
                status = "FAIL"
                notes.append(
                    f"unexpected_pair:{sample_a},{sample_b}"
                )
            else:
                if relatedness < PAIR_REL_PASS:
                    notes.append(
                        f"relatedness<{PAIR_REL_PASS:g}"
                    )
                if concordance < PAIR_CONC_PASS:
                    notes.append(
                        f"concordance<{PAIR_CONC_PASS:g}"
                    )

                status = "PASS" if not notes else "FAIL"

            rows.append(
                [
                    normal,
                    tumor,
                    f"{relatedness:.6f}",
                    str(ibs0),
                    str(ibs2),
                    f"{concordance:.6f}",
                    status,
                    "." if not notes else ";".join(notes),
                ]
            )

        with open(output[0], "w", encoding="utf-8", newline="") as out_fh:
            writer = csv.writer(
                out_fh,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writerow(header)
            writer.writerows(rows)

        n_pass = sum(1 for row in rows if row[6] == "PASS")
        n_fail = len(rows) - n_pass

        with open(log[0], "w", encoding="utf-8") as log_fh:
            log_fh.write(
                f"pairs={len(rows)} PASS={n_pass} FAIL={n_fail}\n"
            )
            log_fh.write(
                f"relatedness_pass={PAIR_REL_PASS}\n"
            )
            log_fh.write(
                f"concordance_pass={PAIR_CONC_PASS}\n"
            )
