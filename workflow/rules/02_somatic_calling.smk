# CPU somatic calling module with noPON/PON switch
# BQSR BAM -> Mutect2 -> orientation model -> normal/tumor pileups
# -> contamination -> FilterMutectCalls
# call_mode=existing reuses an existing raw VCF/index and f1r2 archive.
# In PON mode, the default PON is produced by 02_pon_build.smk.
# It can be overridden with somatic_calling.pon_file in YAML.

SOMATIC_PON_FILE = str(
    SOMATIC.get(
        "pon_file",
        config.get("pon_build", {}).get("pon_vcf", f"{OUT}/pon_db/pon.vcf.gz"),
    )
)

if SOMATIC_ENABLED:

    if CALL_MODE == "cpu":

        rule mutect2_cpu:
            input:
                tumor_bam=lambda wc: bqsr_bam(wc.tumor),
                tumor_bai=lambda wc: bqsr_bai(wc.tumor),
                normal_bam=lambda wc: bqsr_bam(NORMAL_BY_TUMOR[wc.tumor]),
                normal_bai=lambda wc: bqsr_bai(NORMAL_BY_TUMOR[wc.tumor]),
                pon=lambda wc: (
                    [SOMATIC_PON_FILE, f"{SOMATIC_PON_FILE}.tbi"]
                    if PON_MODE == "PON"
                    else []
                )
            output:
                vcf=f"{OUT}/{MUTECT_SUB}/{{tumor}}{SOMATIC_TAG}.vcf.gz",
                vcf_index=f"{OUT}/{MUTECT_SUB}/{{tumor}}{SOMATIC_TAG}.vcf.gz.tbi",
                f1r2=f"{OUT}/{F1R2_SUB}/{{tumor}}{SOMATIC_TAG}_f1r2.tar.gz"
            params:
                normal=lambda wc: NORMAL_BY_TUMOR[wc.tumor],
                pon_arg=(
                    f'-pon "{SOMATIC_PON_FILE}"'
                    if PON_MODE == "PON"
                    else ""
                )
            log:
                f"{OUT}/logs/somatic/{{tumor}}.mutect2.log"
            threads: MUTECT2_THREADS
            resources:
                mem_mb=MUTECT2_MEM_MB
            shell:
                r"""
                set -euo pipefail
                mkdir -p "$(dirname {output.vcf})" \
                         "$(dirname {output.f1r2})" \
                         "$(dirname {log})"

                [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

                "{GATK}" \
                  --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
                  Mutect2 \
                  -R "{REF_FASTA}" \
                  -L "{INTERVALS}" \
                  -germline-resource "{GNOMAD_AF}" {params.pon_arg} \
                  -I "{input.tumor_bam}" \
                  -tumor "{wildcards.tumor}" \
                  -I "{input.normal_bam}" \
                  -normal "{params.normal}" \
                  --native-pair-hmm-threads {threads} \
                  --f1r2-tar-gz "{output.f1r2}" \
                  -O "{output.vcf}" \
                  > "{log}" 2>&1

                test -s "{output.vcf}"
                test -s "{output.vcf_index}"
                test -s "{output.f1r2}"
                """


    rule learn_read_orientation_model:
        input:
            f1r2=lambda wc: raw_f1r2(wc.tumor)
        output:
            model=f"{OUT}/{F1R2_SUB}/{{tumor}}{SOMATIC_TAG}_read-orientation-model.tar.gz"
        log:
            f"{OUT}/logs/somatic/{{tumor}}.orientation_model.log"
        threads: LROM_THREADS
        resources:
            mem_mb=LROM_MEM_MB
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.model})" "$(dirname {log})"
            [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
            "{GATK}" \
              --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
              LearnReadOrientationModel \
              -I "{input.f1r2}" \
              -O "{output.model}" \
              > "{log}" 2>&1
            test -s "{output.model}"
            """


    rule normal_pileup:
        input:
            bam=lambda wc: bqsr_bam(wc.normal),
            bai=lambda wc: bqsr_bai(wc.normal)
        output:
            table=f"{OUT}/{PILEUP_SUB}/{{normal}}_normal{SOMATIC_TAG}-pileups.table"
        log:
            f"{OUT}/logs/somatic/{{normal}}.normal_pileup.log"
        threads: PILEUP_THREADS
        resources:
            mem_mb=PILEUP_MEM_MB
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.table})" "$(dirname {log})"
            [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
            "{GATK}" \
              --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
              GetPileupSummaries \
              -I "{input.bam}" \
              -L "{INTERVALS}" \
              -V "{GNOMAD_AF}" \
              -O "{output.table}" \
              > "{log}" 2>&1
            test -s "{output.table}"
            """


    rule tumor_pileup:
        input:
            bam=lambda wc: bqsr_bam(wc.tumor),
            bai=lambda wc: bqsr_bai(wc.tumor)
        output:
            table=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}_tumor-pileups.table"
        log:
            f"{OUT}/logs/somatic/{{tumor}}.tumor_pileup.log"
        threads: PILEUP_THREADS
        resources:
            mem_mb=PILEUP_MEM_MB
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.table})" "$(dirname {log})"
            [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
            "{GATK}" \
              --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
              GetPileupSummaries \
              -I "{input.bam}" \
              -L "{INTERVALS}" \
              -V "{GNOMAD_AF}" \
              -O "{output.table}" \
              > "{log}" 2>&1
            test -s "{output.table}"
            """


    rule calculate_contamination:
        input:
            tumor=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}_tumor-pileups.table",
            normal=lambda wc: (
                f"{OUT}/{PILEUP_SUB}/"
                f"{NORMAL_BY_TUMOR[wc.tumor]}_normal{SOMATIC_TAG}-pileups.table"
            )
        output:
            table=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}_contamination.table",
            segments=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}.segments.tsv"
        log:
            f"{OUT}/logs/somatic/{{tumor}}.contamination.log"
        threads: CONTAM_THREADS
        resources:
            mem_mb=CONTAM_MEM_MB
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.table})" "$(dirname {log})"
            [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
            "{GATK}" \
              --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
              CalculateContamination \
              -I "{input.tumor}" \
              -matched "{input.normal}" \
              -O "{output.table}" \
              --tumor-segmentation "{output.segments}" \
              > "{log}" 2>&1
            test -s "{output.table}"
            test -s "{output.segments}"
            """


    rule filter_mutect_calls:
        input:
            vcf=lambda wc: raw_somatic_vcf(wc.tumor),
            vcf_index=lambda wc: raw_somatic_vcf_index(wc.tumor),
            contamination=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}_contamination.table",
            segments=f"{OUT}/{PILEUP_SUB}/{{tumor}}{SOMATIC_TAG}.segments.tsv",
            model=f"{OUT}/{F1R2_SUB}/{{tumor}}{SOMATIC_TAG}_read-orientation-model.tar.gz"
        output:
            vcf=f"{OUT}/{MUTECT_SUB}/{{tumor}}{SOMATIC_TAG}_filtered.vcf.gz"
        log:
            f"{OUT}/logs/somatic/{{tumor}}.filter_mutect.log"
        threads: FILTER_THREADS
        resources:
            mem_mb=FILTER_MEM_MB
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"
            [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
            "{GATK}" \
              --java-options "-Xms{SOMATIC_JAVA_MEM_G}G -Xmx{SOMATIC_JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
              FilterMutectCalls \
              -R "{REF_FASTA}" \
              -V "{input.vcf}" \
              --contamination-table "{input.contamination}" \
              --tumor-segmentation "{input.segments}" \
              --ob-priors "{input.model}" \
              -O "{output.vcf}" \
              > "{log}" 2>&1
            test -s "{output.vcf}"
            """
