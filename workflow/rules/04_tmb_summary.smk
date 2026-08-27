# 04_tmb_summary.smk
# Per-sample TMB -> cohort-level TMB summary.
#
# This module assumes 00_common.smk already defines:
#   OUT     : workflow output root
#   TUMORS  : tumor sample IDs from the pair table
#
# Optional config overrides:
#   tools.tmb_summary
#   cohort.label
#
# Existing config values reused:
#   tools.python
#   tmb.high_cutoff

COHORT_CFG = config.get("cohort", {})
TOOLS_CFG = config.get("tools", {})

COHORT_LABEL = str(COHORT_CFG.get("label", config.get("project", "cohort")))

TMB_SUMMARY_SCRIPT = TOOLS_CFG.get(
    "tmb_summary",
    f"{SCRIPT_DIR}/summarize_tmb.py",
)
TMB_SUMMARY_PYTHON = TOOLS_CFG["python"]
COHORT_TMB_HIGH_CUTOFF = float(config.get("tmb", {}).get("high_cutoff", 10.0))

COHORT_DIR = f"{OUT}/cohort_summary"
COHORT_TMB_TABLE = f"{COHORT_DIR}/{COHORT_LABEL}_tmb_summary.tsv"
COHORT_TMB_STATS = f"{COHORT_DIR}/{COHORT_LABEL}_tmb_summary.txt"
COHORT_TMB_PLOT = f"{COHORT_DIR}/{COHORT_LABEL}_tmb_distribution.png"


def all_sample_tmbs(_wildcards):
    return expand(f"{OUT}/tmb/{{tumor}}.tmb.tsv", tumor=TUMORS)


rule summarize_tmb_cohort:
    input:
        all_sample_tmbs
    output:
        table=COHORT_TMB_TABLE,
        stats=COHORT_TMB_STATS,
        plot=COHORT_TMB_PLOT
    log:
        f"{OUT}/logs/cohort_summary/{COHORT_LABEL}.tmb_summary.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120,
        slurm_partition=DEFAULT_PARTITION
    params:
        tmb_dir=f"{OUT}/tmb",
        out_dir=COHORT_DIR,
        label=COHORT_LABEL,
        high_cutoff=COHORT_TMB_HIGH_CUTOFF,
        python=TMB_SUMMARY_PYTHON,
        script=TMB_SUMMARY_SCRIPT
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"

        "{params.python}" "{params.script}" \
          -d "{params.tmb_dir}" \
          -o "{params.out_dir}" \
          --label "{params.label}" \
          --high-cutoff "{params.high_cutoff}" \
          > "{log}" 2>&1

        test -s "{output.table}"
        test -s "{output.stats}"
        test -s "{output.plot}"
        """
