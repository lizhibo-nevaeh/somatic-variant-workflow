rule vep:
    input:
        vcf=filtered_vcf
    output:
        vcf=f"{OUT}/vep/{{tumor}}.vep.vcf",
        warnings=f"{OUT}/vep_warnings/{{tumor}}.vep_warnings.txt"
    log:
        f"{OUT}/logs/vep/{{tumor}}.log"
    threads: 4
    resources:
        mem_mb=40000,
        runtime=720,
        slurm_partition=DEFAULT_PARTITION
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf})" \
                 "$(dirname {output.warnings})" \
                 "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        export PATH="{VEP_PERL_BIN}:$PATH"
        export PERL5LIB="{VEP_PERL5LIB}"

        : > "{output.warnings}"

        "{VEP}" \
          -i "{input.vcf}" \
          -o "{output.vcf}" \
          --vcf --force_overwrite --cache --offline \
          --dir "{VEP_CACHE}" \
          --fasta "{REF_FASTA}" \
          --symbol --terms SO --tsl --hgvs \
          --check_existing \
          --clin_sig_allele 1 \
          --af --af_1kg --af_gnomad \
          --ccds --canonical --variant_class --uniprot \
          --protein --numbers --domains --gene_phenotype \
          --biotype --total_length --allele_number --no_escape \
          --fork {threads} \
          --custom "{CLINVAR}",clinvar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
          --failed 1 \
          --warning_file "{output.warnings}" \
          > "{log}" 2>&1

        test -s "{output.vcf}"
        '''


rule vcf2maf:
    input:
        vcf=f"{OUT}/vep/{{tumor}}.vep.vcf"
    output:
        maf=f"{OUT}/maf/{{tumor}}.vep.maf"
    params:
        normal=lambda wildcards: NORMAL_BY_TUMOR[wildcards.tumor]
    log:
        f"{OUT}/logs/vcf2maf/{{tumor}}.log"
    threads: 2
    resources:
        mem_mb=16000,
        runtime=240,
        slurm_partition=DEFAULT_PARTITION
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.maf})" "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true
        export PATH="{VEP_PERL_BIN}:$PATH"
        export PERL5LIB="{VEP_PERL5LIB}"

        "{VCF2MAF}" \
          --inhibit-vep \
          --input-vcf "{input.vcf}" \
          --output-maf "{output.maf}" \
          --ref-fasta "{REF_FASTA}" \
          --ncbi-build GRCh38 \
          --tumor-id "{wildcards.tumor}" \
          --normal-id "{params.normal}" \
          --retain-ann clinvar,clinvar_CLNSIG,clinvar_CLNREVSTAT,clinvar_CLNDN \
          > "{log}" 2>&1

        test -s "{output.maf}"
        '''


rule tmb:
    input:
        maf=f"{OUT}/maf/{{tumor}}.vep.maf"
    output:
        tsv=f"{OUT}/tmb/{{tumor}}.tmb.tsv",
        stat=f"{OUT}/tmb_filter_stats/{{tumor}}.tmb_filter.stat"
    log:
        f"{OUT}/logs/tmb/{{tumor}}.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120,
        slurm_partition=DEFAULT_PARTITION
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})" \
                 "$(dirname {output.stat})" \
                 "$(dirname {log})"

        [[ -n "{ENV_FILE}" && -s "{ENV_FILE}" ]] && source "{ENV_FILE}" || true

        "{PYTHON}" "{TMB_SCRIPT}" \
          -i "{input.maf}" \
          -o "{output.tsv}" \
          -s "{output.stat}" \
          -t_AF 0.05 \
          -t_depth 25 \
          -n_depth 10 \
          -t_alt_count 3 \
          -n_alt_count 1 \
          -gnomAD_AF 0.001 \
          --exome-size-mb "{EXOME_SIZE_MB}" \
          --high-cutoff "{TMB_HIGH_CUTOFF}" \
          > "{log}" 2>&1

        test -s "{output.tsv}"
        '''
