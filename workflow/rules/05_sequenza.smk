# Sequenza purity/ploidy analysis.
# A single-pair table is generated for each tumor before calling the bundled driver.

SEQ_CFG = config.get("sequenza", {})

SEQ_PIPELINE = SEQ_CFG.get(
    "pipeline_script",
    f"{SCRIPT_DIR}/run_sequenza.sh"
)
SEQ_BQSR_DIR = SEQ_CFG.get("bqsr_dir", BQSR_BAM_DIR)
SEQ_BAM_SUFFIX = str(SEQ_CFG.get("bam_suffix", ".BQSR.bam"))
SEQ_REF_FASTA = SEQ_CFG.get("ref_fasta", "")
SEQ_GC_WIGGLE = SEQ_CFG.get("gc_wiggle", "")
SEQ_BIN_DIR = SEQ_CFG.get("sequenza_bin_dir", "/path/to/sequenza/bin")
SEQ_UTILS = SEQ_CFG.get(
    "sequenza_utils",
    f"{SEQ_BIN_DIR}/sequenza-utils"
)
SEQ_RSCRIPT = SEQ_CFG.get("rscript", "Rscript")
SEQ_DEFAULT_GENDER = str(SEQ_CFG.get("default_gender", "male"))
SEQ_BIN_SIZE = int(SEQ_CFG.get("bin_size", 50))
SEQ_THREADS = int(SEQ_CFG.get("threads", 8))
SEQ_MEM_GB = int(SEQ_CFG.get("mem_gb", 64))
SEQ_MEM_MB = int(SEQ_CFG.get("mem_mb", SEQ_MEM_GB * 1000))
SEQ_RUNTIME = int(SEQ_CFG.get("runtime", 1440))
SEQ_WALLTIME = str(SEQ_CFG.get("walltime", "24:00:00"))
SEQ_QUEUE = str(SEQ_CFG.get("queue", DEFAULT_PARTITION))


rule sequenza:
    input:
        normal=lambda w: (
            f"{SEQ_BQSR_DIR}/{NORMAL_BY_TUMOR[w.tumor]}"
            f"{SEQ_BAM_SUFFIX}"
        ),
        normal_bai=lambda w: (
            f"{SEQ_BQSR_DIR}/{NORMAL_BY_TUMOR[w.tumor]}"
            f"{SEQ_BAM_SUFFIX}.bai"
        ),
        tumor=lambda w: (
            f"{SEQ_BQSR_DIR}/{w.tumor}{SEQ_BAM_SUFFIX}"
        ),
        tumor_bai=lambda w: (
            f"{SEQ_BQSR_DIR}/{w.tumor}{SEQ_BAM_SUFFIX}.bai"
        ),
        pairs=PAIR_FILE
    output:
        purity=(
            f"{OUT}/sequenza/{{tumor}}/"
            f"{{tumor}}.purity_ploidy.txt"
        )
    log:
        f"{OUT}/logs/sequenza/{{tumor}}.log"
    threads:
        SEQ_THREADS
    resources:
        mem_mb=SEQ_MEM_MB,
        runtime=SEQ_RUNTIME,
        slurm_partition=SEQ_QUEUE
    params:
        pipeline=SEQ_PIPELINE,
        bqsr_dir=SEQ_BQSR_DIR,
        bam_suffix=SEQ_BAM_SUFFIX,
        ref_fasta=SEQ_REF_FASTA,
        gc_wiggle=SEQ_GC_WIGGLE,
        bin_dir=SEQ_BIN_DIR,
        utils=SEQ_UTILS,
        rscript=SEQ_RSCRIPT,
        default_gender=SEQ_DEFAULT_GENDER,
        bin_size=SEQ_BIN_SIZE,
        mem_gb=SEQ_MEM_GB,
        walltime=SEQ_WALLTIME,
        queue=SEQ_QUEUE
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
          "$(dirname "{log}")" \
          "{OUT}/05_sequenza/sh" \
          "{OUT}/05_sequenza/pairs"

        child="{OUT}/05_sequenza/sh/{wildcards.tumor}.step5_sequenza.sh"
        single_pair="{OUT}/05_sequenza/pairs/{wildcards.tumor}.pair.tsv"

        rm -f "$child"

        # 从总配对表提取当前 tumor 的完整记录；若有第3列 gender 会一并保留。
        tr -d '\r' < "{input.pairs}" |
        awk -v tumor="{wildcards.tumor}" '
          BEGIN {{ OFS="\t" }}
          NF >= 2 && $2 == tumor {{
              print
              found = 1
              exit
          }}
          END {{
              if (!found) exit 1
          }}
        ' > "$single_pair"

        test -s "$single_pair"

        BASE_OUT="{OUT}" \
        PAIR_LIST="$single_pair" \
        ENV_FILE="{ENV_FILE}" \
        BQSR_DIR="{params.bqsr_dir}" \
        BAM_SUFFIX="{params.bam_suffix}" \
        REF_FASTA="{params.ref_fasta}" \
        GC_WIGGLE="{params.gc_wiggle}" \
        BIN_SIZE="{params.bin_size}" \
        SEQUENZA_BIN_DIR="{params.bin_dir}" \
        SEQUENZA_UTILS="{params.utils}" \
        RSCRIPT_BIN="{params.rscript}" \
        DEFAULT_GENDER="{params.default_gender}" \
        QUEUE="{params.queue}" \
        CORE="{threads}" \
        MEM="{params.mem_gb}" \
        WALLTIME="{params.walltime}" \
        SUBMIT=0 \
        bash "{params.pipeline}" "$single_pair" \
          > "{log}" 2>&1

        test -s "$child"

        bash "$child" >> "{log}" 2>&1

        test -s "{output.purity}"
        """
