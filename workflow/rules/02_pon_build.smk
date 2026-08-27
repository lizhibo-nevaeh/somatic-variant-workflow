# 02_pon_build.smk
# Optional PON construction module.
#
# Workflow:
#   unique normal BQSR BAM
#       -> GATK Mutect2 (normal-as-tumor, --max-mnp-distance 0)
#       -> normal/{normal}.vcf.gz
#   all normal VCFs
#       -> GenomicsDBImport
#       -> CreateSomaticPanelOfNormals
#       -> pon_db/pon.vcf.gz
#
# This module only BUILDS the PON.
# Whether tumor-normal calling uses the PON is handled in 02_somatic_calling.smk.
#
# Expected globals from 00_common.smk:
#   OUT, NORMAL_BY_TUMOR, BQSR_BAM_DIR, REFS,
#   INTERVALS, GNOMAD_AF, GATK, ENV_FILE

PON_BUILD_CFG = config.get("pon_build", {})

PON_REF_FASTA = str(REFS["fasta"])
PON_INTERVALS = str(INTERVALS)
PON_GNOMAD_AF = str(GNOMAD_AF)

PON_NORMALS = sorted(set(NORMAL_BY_TUMOR.values()))

PON_STEP2A_THREADS = int(PON_BUILD_CFG.get("step2a_threads", 16))
PON_STEP2A_MEM_MB = int(PON_BUILD_CFG.get("step2a_mem_mb", 50000))
PON_STEP2A_RUNTIME = int(PON_BUILD_CFG.get("step2a_runtime", 60000))

PON_STEP2B_THREADS = int(PON_BUILD_CFG.get("step2b_threads", 24))
PON_STEP2B_MEM_MB = int(PON_BUILD_CFG.get("step2b_mem_mb", 100000))
PON_STEP2B_RUNTIME = int(PON_BUILD_CFG.get("step2b_runtime", 60000))

PON_QUEUE = str(PON_BUILD_CFG.get("queue", DEFAULT_PARTITION))
PON_BATCH_SIZE = int(PON_BUILD_CFG.get("batch_size", 50))

PON_NORMAL_VCF_DIR = str(PON_BUILD_CFG.get("normal_vcf_dir", f"{OUT}/normal"))
PON_WORKSPACE = str(PON_BUILD_CFG.get("workspace", f"{OUT}/pon_db_workspace"))
PON_FINAL_VCF = str(PON_BUILD_CFG.get("pon_vcf", f"{OUT}/pon_db/pon.vcf.gz"))

PON_NORMAL_VCFS = [
    f"{PON_NORMAL_VCF_DIR}/{normal}.vcf.gz"
    for normal in PON_NORMALS
]

PON_VCF_ARGS = " ".join(
    f'-V "{vcf}"'
    for vcf in PON_NORMAL_VCFS
)


rule pon_normal_mutect2:
    input:
        bam=lambda w: f"{BQSR_BAM_DIR}/{w.normal}.BQSR.bam",
        bai=lambda w: f"{BQSR_BAM_DIR}/{w.normal}.BQSR.bam.bai",
        ref=PON_REF_FASTA,
    output:
        vcf=f"{PON_NORMAL_VCF_DIR}/{{normal}}.vcf.gz",
        tbi=f"{PON_NORMAL_VCF_DIR}/{{normal}}.vcf.gz.tbi",
    log:
        f"{OUT}/logs/pon/step2a/{{normal}}.log"
    threads:
        PON_STEP2A_THREADS
    resources:
        mem_mb=PON_STEP2A_MEM_MB,
        runtime=PON_STEP2A_RUNTIME,
        slurm_partition=PON_QUEUE
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.vcf}")" "$(dirname "{log}")"

        [[ -s "{ENV_FILE}" ]] && source "{ENV_FILE}"

        "{GATK}" \
          --java-options "-Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -XX:ParallelGCThreads={threads}" \
          Mutect2 \
          -R "{input.ref}" \
          -I "{input.bam}" \
          --max-mnp-distance 0 \
          --native-pair-hmm-threads {threads} \
          -O "{output.vcf}" \
          > "{log}" 2>&1

        test -s "{output.vcf}"
        test -s "{output.tbi}"
        """


rule pon_genomicsdb_import:
    input:
        vcfs=PON_NORMAL_VCFS,
        ref=PON_REF_FASTA,
        intervals=PON_INTERVALS,
    output:
        workspace=directory(PON_WORKSPACE)
    log:
        f"{OUT}/logs/pon/GenomicsDBImport.log"
    threads:
        PON_STEP2B_THREADS
    resources:
        mem_mb=PON_STEP2B_MEM_MB,
        runtime=PON_STEP2B_RUNTIME,
        slurm_partition=PON_QUEUE
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{OUT}/tmp" "$(dirname "{log}")"

        [[ -s "{ENV_FILE}" ]] && source "{ENV_FILE}"

        "{GATK}" \
          --java-options "-XX:ParallelGCThreads={threads} -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={OUT}/tmp" \
          GenomicsDBImport \
          -R "{input.ref}" \
          --batch-size {PON_BATCH_SIZE} \
          --genomicsdb-shared-posixfs-optimizations true \
          --bypass-feature-reader \
          -L "{input.intervals}" \
          --genomicsdb-workspace-path "{output.workspace}" \
          {PON_VCF_ARGS} \
          > "{log}" 2>&1

        test -d "{output.workspace}"
        """


rule pon_create:
    input:
        workspace=PON_WORKSPACE,
        ref=PON_REF_FASTA,
        gnomad=PON_GNOMAD_AF,
        gnomad_tbi=f"{PON_GNOMAD_AF}.tbi",
    output:
        vcf=PON_FINAL_VCF,
        tbi=f"{PON_FINAL_VCF}.tbi",
    log:
        f"{OUT}/logs/pon/CreateSomaticPanelOfNormals.log"
    threads:
        PON_STEP2B_THREADS
    resources:
        mem_mb=PON_STEP2B_MEM_MB,
        runtime=PON_STEP2B_RUNTIME,
        slurm_partition=PON_QUEUE
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.vcf}")" "$(dirname "{log}")"

        [[ -s "{ENV_FILE}" ]] && source "{ENV_FILE}"

        cd "{OUT}"

        "{GATK}" \
          --java-options "-XX:ParallelGCThreads={threads} -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={OUT}/tmp" \
          CreateSomaticPanelOfNormals \
          -R "{input.ref}" \
          --germline-resource "{input.gnomad}" \
          -V "gendb://$(basename "{input.workspace}")" \
          -O "{output.vcf}" \
          > "{log}" 2>&1

        test -s "{output.vcf}"
        test -s "{output.tbi}"
        """
