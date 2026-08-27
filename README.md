# Somatic Variant Workflow

A modular Snakemake workflow for paired tumor-normal WES/WGS analysis.

## Workflow

```text
FASTQ
  └─ fastp → BWA → duplicate marking → BQSR
                                  │
                                  ├─ BAM QC
                                  ├─ Somalier pair QC
                                  ├─ normal contamination QC
                                  ├─ Mutect2 → filtering → VEP → MAF → TMB
                                  ├─ Sequenza → CNVkit
                                  └─ MSIsensor-pro

Optional: GATK Panel of Normals (PON) construction
```

## Modules

| File | Function |
|---|---|
| `01_preprocessing.smk` | FASTQ QC, alignment, duplicate marking and BQSR |
| `01_5_qc.smk` | BAM-level `samtools flagstat` and `samtools stats` |
| `01_6_pair_qc.smk` | Tumor-normal identity QC with Somalier |
| `01_7_normal_contam_qc.smk` | Normal-sample contamination report |
| `02_pon_build.smk` | Optional GATK Panel of Normals construction |
| `02_somatic_calling.smk` | Mutect2 calling and filtering |
| `03_annotation_tmb.smk` | VEP, vcf2maf and TMB calculation |
| `04_tmb_summary.smk` | Cohort-level TMB summary |
| `05_sequenza.smk` | Tumor purity and ploidy estimation |
| `05_cnv.smk` | CNVkit copy-number analysis |
| `06_msi.smk` | MSI analysis with MSIsensor-pro |

## Requirements

Snakemake, fastp, BWA, samtools/sambamba, GATK, VEP, vcf2maf, Somalier, Sequenza, CNVkit, MSIsensor-pro, Python 3 and R.

Reference files and tool paths are configured in `config/config.example.yaml`.

## Input

Tumor-normal pairs are provided as a tab-delimited file:

```text
# normal    tumor    sex
NORMAL001   TUMOR001  male
NORMAL002   TUMOR002  female
```

The third column is optional for most modules; provide it when sex-aware Sequenza/CNVkit analysis is needed.

Copy the example configuration and edit paths/settings for your environment:

```bash
cp config/config.example.yaml config/config.yaml
```

## Usage

Preview the workflow:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --dry-run \
  --printshellcmds
```

Run with local cores:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --cores 8
```

For HPC execution, use a Snakemake executor/profile appropriate for your scheduler.

The default `rule all` produces per-tumor TMB results. Individual modules can also be requested through their rules or output files.

## Notes

Large reference files and software installations are not included. Review reference versions, parameters and compute resources before running the workflow.
