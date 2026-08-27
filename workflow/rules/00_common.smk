from pathlib import Path

OUT = config["output_root"]
SCRIPT_DIR = str(Path(workflow.basedir).resolve().parent / "scripts")
SOURCE_ROOT = config.get("source_root", "")
PAIR_FILE = config["pairs"]

MODE = config.get("mode", {})
PON_MODE = MODE.get("pon", "noPON")
SEQ_TYPE = MODE.get("seq_type", "WES")

SOMATIC = config.get("somatic_calling", {})
SOMATIC_ENABLED = bool(SOMATIC.get("enabled", False))
CALL_MODE = SOMATIC.get("call_mode", "cpu")

if SOMATIC_ENABLED:
    if PON_MODE not in {"noPON", "PON"}:
        raise ValueError(
            "mode.pon must be 'noPON' or 'PON'; "
            f"received mode.pon={PON_MODE!r}"
        )
    if CALL_MODE not in {"cpu", "existing"}:
        raise ValueError(
            "v0.3.1 supports call_mode=cpu or existing; "
            f"received somatic_calling.call_mode={CALL_MODE!r}"
        )


EXISTING_CALLING = SOMATIC.get("existing_calling", {})
EXISTING_RAW_VCF_DIR = EXISTING_CALLING.get("raw_vcf_dir", "")
EXISTING_F1R2_DIR = EXISTING_CALLING.get("f1r2_dir", "")

if SOMATIC_ENABLED and CALL_MODE == "existing":
    if not EXISTING_RAW_VCF_DIR or not EXISTING_F1R2_DIR:
        raise ValueError(
            "call_mode=existing requires somatic_calling.existing_calling."
            "raw_vcf_dir and f1r2_dir"
        )

INPUTS = config.get("inputs", {})
FILTERED_VCF_DIR = INPUTS.get("filtered_vcf_dir", "")
FILTERED_VCF_SUFFIX = INPUTS.get(
    "filtered_vcf_suffix",
    "_noPON_filtered.vcf.gz" if PON_MODE == "noPON" else "_filtered.vcf.gz",
)

REFS = config["references"]
REF_FASTA = REFS["fasta"]
CLINVAR = REFS["clinvar"]
VEP_CACHE = REFS["vep_cache"]
REF_ROOT = str(Path(REF_FASTA).parent)
KNOWN_SITES = REFS.get(
    "known_sites",
    [
        f"{REF_ROOT}/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        f"{REF_ROOT}/dbsnp_146.hg38.vcf.gz",
        f"{REF_ROOT}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        f"{REF_ROOT}/1000G_omni2.5.hg38.vcf.gz",
    ],
)
KNOWN_SITE_ARGS = " ".join(
    f'--known-sites "{path}"' for path in KNOWN_SITES
)
INTERVALS = REFS.get(
    "intervals",
    f"{REF_ROOT}/hg38_v0_wgs_calling_regions.hg38.interval_list",
)
GNOMAD_AF = REFS.get(
    "gnomad_af",
    f"{REF_ROOT}/af-only-gnomad.hg38.ref218.header218.vcf.gz",
)

TOOLS = config["tools"]
ENV_FILE = str(TOOLS.get("env_file", ""))
FASTP = str(TOOLS.get("fastp", "fastp"))
BWA = str(TOOLS.get("bwa", "bwa"))
SAMTOOLS = str(TOOLS.get("samtools", "samtools"))
SAMBAMBA = str(TOOLS.get("sambamba", "sambamba"))
VEP = str(TOOLS.get("vep", "vep"))
VCF2MAF = str(TOOLS.get("vcf2maf", "vcf2maf.pl"))
PYTHON = str(TOOLS.get("python", "python3"))
TMB_SCRIPT = str(TOOLS.get("tmb_script", f"{SCRIPT_DIR}/calculate_tmb.py"))
GATK = str(TOOLS.get("gatk", "gatk"))

CLUSTER = config.get("cluster", {})
DEFAULT_PARTITION = str(CLUSTER.get("partition", "compute"))

EXOME_SIZE_MB = config["tmb"]["exome_size_mb"]
TMB_HIGH_CUTOFF = config["tmb"]["high_cutoff"]

PREP = config.get("preprocessing", {})
KEEP_CLEAN_FASTQ = bool(PREP.get("keep_clean_fastq", False))
FASTQ_DIR = PREP.get("fastq_dir", f"{SOURCE_ROOT}/data")
FASTQ_R1_PATTERN = PREP.get(
    "fastq_r1_pattern",
    "{fastq_dir}/{sample}/{sample}_f1.fastq.gz",
)
FASTQ_R2_PATTERN = PREP.get(
    "fastq_r2_pattern",
    "{fastq_dir}/{sample}/{sample}_r2.fastq.gz",
)
RG_LB = PREP.get("rg_lb", SEQ_TYPE)
RG_PLATFORM = PREP.get("rg_platform", "Illumina")
FASTP_MIN_LENGTH = int(PREP.get("fastp_min_length", 50))
FASTP_EXTRA_ARGS = PREP.get("fastp_extra_args", "-5 -3")
BWA_MIN_SCORE = int(PREP.get("bwa_min_score", 0))
JAVA_MEM_G = int(PREP.get("java_mem_g", 32))

RESOURCE_CFG = config.get("resources", {})

def resource_value(rule_name, key, default):
    return RESOURCE_CFG.get(rule_name, {}).get(key, default)

FASTP_THREADS = int(resource_value("fastp", "threads", 8))
FASTP_MEM_MB = int(resource_value("fastp", "mem_mb", 8000))
BWA_THREADS = int(resource_value("bwa_sort", "threads", 8))
BWA_MEM_MB = int(resource_value("bwa_sort", "mem_mb", 32000))
MARKDUP_THREADS = int(resource_value("markdup", "threads", 8))
MARKDUP_MEM_MB = int(resource_value("markdup", "mem_mb", 32000))
BQSR_THREADS = int(resource_value("bqsr", "threads", 8))
BQSR_MEM_MB = int(resource_value("bqsr", "mem_mb", 40000))

MUTECT2_THREADS = int(resource_value("mutect2", "threads", 16))
MUTECT2_MEM_MB = int(resource_value("mutect2", "mem_mb", 96000))
LROM_THREADS = int(resource_value("learn_orientation", "threads", 16))
LROM_MEM_MB = int(resource_value("learn_orientation", "mem_mb", 96000))
PILEUP_THREADS = int(resource_value("pileup", "threads", 16))
PILEUP_MEM_MB = int(resource_value("pileup", "mem_mb", 96000))
CONTAM_THREADS = int(resource_value("contamination", "threads", 16))
CONTAM_MEM_MB = int(resource_value("contamination", "mem_mb", 96000))
FILTER_THREADS = int(resource_value("filter_mutect", "threads", 16))
FILTER_MEM_MB = int(resource_value("filter_mutect", "mem_mb", 96000))
SOMATIC_JAVA_MEM_G = int(SOMATIC.get("java_mem_g", 80))

BQSR_BAM_DIR = SOMATIC.get("bqsr_bam_dir", f"{OUT}/BQSR")

if PON_MODE == "noPON":
    SOMATIC_TAG = "_noPON"
    MUTECT_SUB = "mutect2_noPON"
    F1R2_SUB = "f1r2_noPON"
    PILEUP_SUB = "pileup_noPON"
else:
    SOMATIC_TAG = ""
    MUTECT_SUB = "mutect2"
    F1R2_SUB = "f1r2"
    PILEUP_SUB = "pileup"

VEP_PERL_BIN = str(TOOLS.get("vep_perl_bin", ""))
VEP_PERL5LIB = str(TOOLS.get("vep_perl5lib", ""))

pairs = []
with open(PAIR_FILE, encoding="utf-8") as handle:
    for lineno, line in enumerate(handle, 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(
                f"{PAIR_FILE}:{lineno}: expected at least 2 columns "
                "(normal tumor)"
            )
        normal, tumor = fields[0], fields[1]
        pairs.append((normal, tumor))

if not pairs:
    raise ValueError(f"No usable pairs found in {PAIR_FILE}")

NORMAL_BY_TUMOR = {tumor: normal for normal, tumor in pairs}
TUMORS = list(NORMAL_BY_TUMOR)
NORMALS = sorted({normal for normal, _ in pairs})
SAMPLES = sorted({sample for pair in pairs for sample in pair})


def filtered_vcf(wildcards):
    if SOMATIC_ENABLED:
        return (
            f"{OUT}/{MUTECT_SUB}/"
            f"{wildcards.tumor}{SOMATIC_TAG}_filtered.vcf.gz"
        )
    if not FILTERED_VCF_DIR:
        raise ValueError(
            "inputs.filtered_vcf_dir is required when "
            "somatic_calling.enabled is false"
        )
    return (
        f"{FILTERED_VCF_DIR}/"
        f"{wildcards.tumor}{FILTERED_VCF_SUFFIX}"
    )


def fastq_r1(wildcards):
    return FASTQ_R1_PATTERN.format(
        fastq_dir=FASTQ_DIR,
        sample=wildcards.sample,
    )


def fastq_r2(wildcards):
    return FASTQ_R2_PATTERN.format(
        fastq_dir=FASTQ_DIR,
        sample=wildcards.sample,
    )


def bqsr_bam(sample):
    return f"{BQSR_BAM_DIR}/{sample}.BQSR.bam"


def bqsr_bai(sample):
    return f"{BQSR_BAM_DIR}/{sample}.BQSR.bam.bai"

def raw_somatic_vcf(tumor):
    if CALL_MODE == "cpu":
        return f"{OUT}/{MUTECT_SUB}/{tumor}{SOMATIC_TAG}.vcf.gz"
    return f"{EXISTING_RAW_VCF_DIR}/{tumor}{SOMATIC_TAG}.vcf.gz"


def raw_somatic_vcf_index(tumor):
    return f"{raw_somatic_vcf(tumor)}.tbi"


def raw_f1r2(tumor):
    if CALL_MODE == "cpu":
        return f"{OUT}/{F1R2_SUB}/{tumor}{SOMATIC_TAG}_f1r2.tar.gz"
    return f"{EXISTING_F1R2_DIR}/{tumor}{SOMATIC_TAG}_f1r2.tar.gz"

