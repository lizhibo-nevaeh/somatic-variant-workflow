# 06 MSI：MSIsensor-pro paired tumor-normal analysis

MSI_CFG = config.get("msi", {})

MSI_BIN = MSI_CFG.get("msisensor", "msisensor-pro")
MSI_SITES = MSI_CFG.get("sites", "")
MSI_BQSR_DIR = MSI_CFG.get("bqsr_dir", BQSR_BAM_DIR)
MSI_BED = str(MSI_CFG.get("bed", ""))
MSI_COV = int(MSI_CFG.get(
    "coverage_threshold",
    20 if SEQ_TYPE == "WES" else 15
))
MSI_FDR = float(MSI_CFG.get("fdr", 0.05))
MSI_THREADS = int(MSI_CFG.get("threads", 4))
MSI_MEM = int(MSI_CFG.get("mem_mb", 16000))
MSI_RUNTIME = int(MSI_CFG.get("runtime", 120))
MSI_QUEUE = str(MSI_CFG.get("queue", DEFAULT_PARTITION))
MSI_BED_ARG = f'-e "{MSI_BED}"' if MSI_BED else ""


rule msi:
    input:
        normal=lambda w: (
            f"{MSI_BQSR_DIR}/{NORMAL_BY_TUMOR[w.tumor]}.BQSR.bam"
        ),
        normal_bai=lambda w: (
            f"{MSI_BQSR_DIR}/{NORMAL_BY_TUMOR[w.tumor]}.BQSR.bam.bai"
        ),
        tumor=lambda w: f"{MSI_BQSR_DIR}/{w.tumor}.BQSR.bam",
        tumor_bai=lambda w: f"{MSI_BQSR_DIR}/{w.tumor}.BQSR.bam.bai",
        sites=MSI_SITES
    output:
        main=f"{OUT}/msi/{{tumor}}",
        all_sites=f"{OUT}/msi/{{tumor}}_all",
        distribution=f"{OUT}/msi/{{tumor}}_dis",
        unstable=f"{OUT}/msi/{{tumor}}_unstable"
    log:
        f"{OUT}/logs/msi/{{tumor}}.log"
    threads:
        MSI_THREADS
    resources:
        mem_mb=MSI_MEM,
        runtime=MSI_RUNTIME,
        slurm_partition=MSI_QUEUE
    params:
        exe=MSI_BIN,
        cov=MSI_COV,
        fdr=MSI_FDR,
        bed_arg=MSI_BED_ARG
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.main}")" "$(dirname "{log}")"

        "{params.exe}" msi \
          -d "{input.sites}" \
          -n "{input.normal}" \
          -t "{input.tumor}" \
          -o "{output.main}" \
          -c {params.cov} \
          -b {threads} \
          -f {params.fdr} \
          {params.bed_arg} \
          > "{log}" 2>&1

        test -s "{output.main}"
        test -s "{output.all_sites}"
        test -s "{output.distribution}"
        test -s "{output.unstable}"
        """
