# CPU preprocessing module
# FASTQ -> fastp -> BWA/sort -> sambamba markdup -> BQSR

rule fastp:
    input:
        r1=fastq_r1,
        r2=fastq_r2
    output:
        r1=(
            f"{OUT}/tmp/{{sample}}_clean_1.fq.gz"
            if KEEP_CLEAN_FASTQ
            else temp(f"{OUT}/tmp/{{sample}}_clean_1.fq.gz")
        ),
        r2=(
            f"{OUT}/tmp/{{sample}}_clean_2.fq.gz"
            if KEEP_CLEAN_FASTQ
            else temp(f"{OUT}/tmp/{{sample}}_clean_2.fq.gz")
        ),
        json=f"{OUT}/raw_bam/{{sample}}.fastp.json",
        html=f"{OUT}/raw_bam/{{sample}}.fastp.html"
    log:
        f"{OUT}/logs/step1/{{sample}}.fastp.log"
    threads: FASTP_THREADS
    resources:
        mem_mb=FASTP_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})" \
                 "$(dirname {output.json})" \
                 "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        FASTP="{FASTP}"

        "$FASTP" \
          -l {FASTP_MIN_LENGTH} \
          -w {threads} \
          {FASTP_EXTRA_ARGS} \
          -i "{input.r1}" \
          -I "{input.r2}" \
          -o "{output.r1}" \
          -O "{output.r2}" \
          -j "{output.json}" \
          -h "{output.html}" \
          > "{log}" 2>&1

        test -s "{output.r1}"
        test -s "{output.r2}"
        test -s "{output.json}"
        test -s "{output.html}"
        """


rule bwa_sort:
    input:
        r1=f"{OUT}/tmp/{{sample}}_clean_1.fq.gz",
        r2=f"{OUT}/tmp/{{sample}}_clean_2.fq.gz"
    output:
        bam=f"{OUT}/raw_bam/{{sample}}.bam"
    log:
        f"{OUT}/logs/step1/{{sample}}.bwa_sort.log"
    threads: BWA_THREADS
    resources:
        mem_mb=BWA_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        BWA="{BWA}"
        SAMTOOLS="{SAMTOOLS}"

        "$BWA" mem \
          -t {threads} \
          -T {BWA_MIN_SCORE} \
          -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:{RG_LB}\tPL:{RG_PLATFORM}\tPU:{wildcards.sample}" \
          "{REF_FASTA}" \
          "{input.r1}" \
          "{input.r2}" \
          2> "{log}" \
        | "$SAMTOOLS" sort \
            -@ {threads} \
            -o "{output.bam}" \
            -

        test -s "{output.bam}"
        "$SAMTOOLS" quickcheck -v "{output.bam}"
        """


rule markdup:
    input:
        bam=f"{OUT}/raw_bam/{{sample}}.bam"
    output:
        bam=f"{OUT}/marked/{{sample}}.marked.bam",
        bai=f"{OUT}/marked/{{sample}}.marked.bam.bai"
    log:
        f"{OUT}/logs/step1/{{sample}}.markdup.log"
    threads: MARKDUP_THREADS
    resources:
        mem_mb=MARKDUP_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" \
                 "{OUT}/tmp/markdup" \
                 "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        SAMTOOLS="{SAMTOOLS}"
        SAMBAMBA="{SAMBAMBA}"

        cur=$(ulimit -n 2>/dev/null || echo 1024)
        hard=$(ulimit -Hn 2>/dev/null || echo 65535)
        [[ "$hard" =~ ^[0-9]+$ ]] || hard=65535
        target=$(( hard < 65535 ? hard : 65535 ))
        [[ "$cur" -lt "$target" ]] \
          && ulimit -n "$target" 2>/dev/null \
          || true

        if ! "$SAMBAMBA" markdup \
          -t {threads} \
          --tmpdir="{OUT}/tmp/markdup" \
          "{input.bam}" \
          "{output.bam}" \
          > "{log}" 2>&1
        then
          echo "[WARN] markdup failed; retrying with 4 threads" \
            >> "{log}"
          "$SAMBAMBA" markdup \
            -t 4 \
            --tmpdir="{OUT}/tmp/markdup" \
            "{input.bam}" \
            "{output.bam}" \
            >> "{log}" 2>&1
        fi

        "$SAMTOOLS" index -@ {threads} "{output.bam}"
        "$SAMTOOLS" quickcheck -v "{output.bam}"
        test -s "{output.bam}"
        test -s "{output.bai}"
        """


rule base_recalibrator:
    input:
        bam=f"{OUT}/marked/{{sample}}.marked.bam",
        bai=f"{OUT}/marked/{{sample}}.marked.bam.bai"
    output:
        table=f"{OUT}/BQSR/{{sample}}_recal_data.table1"
    log:
        f"{OUT}/logs/step1/{{sample}}.base_recalibrator.log"
    threads: BQSR_THREADS
    resources:
        mem_mb=BQSR_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.table})" "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

        "{GATK}" \
          --java-options "-Xms{JAVA_MEM_G}G -Xmx{JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
          BaseRecalibrator \
          -R "{REF_FASTA}" \
          {KNOWN_SITE_ARGS} \
          -I "{input.bam}" \
          -O "{output.table}" \
          > "{log}" 2>&1

        test -s "{output.table}"
        """


rule apply_bqsr:
    input:
        bam=f"{OUT}/marked/{{sample}}.marked.bam",
        bai=f"{OUT}/marked/{{sample}}.marked.bam.bai",
        table=f"{OUT}/BQSR/{{sample}}_recal_data.table1"
    output:
        bam=f"{OUT}/BQSR/{{sample}}.BQSR.bam",
        bai=f"{OUT}/BQSR/{{sample}}.BQSR.bam.bai"
    log:
        f"{OUT}/logs/step1/{{sample}}.apply_bqsr.log"
    threads: BQSR_THREADS
    resources:
        mem_mb=BQSR_MEM_MB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        SAMTOOLS="{SAMTOOLS}"

        "{GATK}" \
          --java-options "-Xms{JAVA_MEM_G}G -Xmx{JAVA_MEM_G}G -XX:ParallelGCThreads={threads}" \
          ApplyBQSR \
          -R "{REF_FASTA}" \
          -I "{input.bam}" \
          --bqsr-recal-file "{input.table}" \
          -O "{output.bam}" \
          > "{log}" 2>&1

        "$SAMTOOLS" index -@ {threads} "{output.bam}"
        "$SAMTOOLS" quickcheck -v "{output.bam}"
        test -s "{output.bam}"
        test -s "{output.bai}"
        """
