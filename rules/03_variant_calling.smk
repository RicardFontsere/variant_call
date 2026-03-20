# =============================================================================
# 03 - VARIANT CALLING
# 
# Workflow:
# 1. HaplotypeCaller per sample per interval -> GVCFs
# 2. GenomicsDBImport per interval -> GenomicsDB
# 3. GenotypeGVCFs per interval -> VCFs
# 4. SortVcf per interval -> sorted VCFs
# 5. MergeVcfs -> raw.vcf
# 6. VariantFiltration -> annotate filters
# 7. SelectVariants -> apply filters + biallelic SNPs only
# 8. Convert to BCF + index
# =============================================================================

rule haplotypecaller:
    input:
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam"),
        bai = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam.csi"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict",
        interval = os.path.join(INTERVALS_DIR, "{interval}.interval_list")
    output:
        gvcf = temp(os.path.join(RESULTS_DIR, "03_variants", "gvcf", "{sample}", "{interval}.g.vcf"))
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=4000,
        runtime=2400
    log:
        os.path.join(RESULTS_DIR, "logs", "03_haplotypecaller", "{sample}.{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus_per_task}
        mkdir -p $(dirname {output.gvcf})
        mkdir -p $(dirname {log})
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -L {input.interval} \
            -ERC GVCF \
            --native-pair-hmm-threads {resources.cpus_per_task} 2> {log}
        """

rule genomicsdb_import:
    input:
        gvcfs = lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "03_variants", "gvcf", "{sample}", "{interval}.g.vcf"),
            sample=SAMPLES,
            interval=wildcards.interval
        ),
        interval = os.path.join(INTERVALS_DIR, "{interval}.interval_list"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        db = directory(os.path.join(RESULTS_DIR, "03_variants", "genomicsdb", "{interval}"))
    params:
        gvcf_args = lambda wildcards, input: " ".join([f"-V {g}" for g in input.gvcfs])
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "03_genomicsdb_import", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.db})
        mkdir -p $(dirname {output.db})/tmp
        mkdir -p $(dirname {log})
        gatk GenomicsDBImport \
            {params.gvcf_args} \
            --genomicsdb-workspace-path {output.db} \
            -L {input.interval} \
            --reader-threads {resources.cpus_per_task} \
            --batch-size 50 \
            --tmp-dir $(dirname {output.db})/tmp 2> {log}
        """

rule genotype_gvcfs:
    input:
        db = os.path.join(RESULTS_DIR, "03_variants", "genomicsdb", "{interval}"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = temp(os.path.join(RESULTS_DIR, "03_variants", "genotyped", "{interval}.vcf")),
        idx = temp(os.path.join(RESULTS_DIR, "03_variants", "genotyped", "{interval}.vcf.idx"))
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=20000,
        runtime=1000
    log:
        os.path.join(RESULTS_DIR, "logs", "03_genotype_gvcfs", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            --create-output-variant-index true 2> {log}
        """

rule sort_vcf:
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "genotyped", "{interval}.vcf"),
        dict = REF_PREFIX + ".dict"
    output:
        vcf = temp(os.path.join(RESULTS_DIR, "03_variants", "sorted", "{interval}.vcf")),
        idx = temp(os.path.join(RESULTS_DIR, "03_variants", "sorted", "{interval}.vcf.idx"))
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "03_sort_vcf", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk SortVcf \
            -I {input.vcf} \
            -O {output.vcf} \
            -SD {input.dict} 2> {log}
        """

rule merge_vcfs:
    input:
        vcfs = lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "03_variants", "sorted", "{interval}.vcf"),
            interval=get_intervals()
        ),
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf")
    params:
        vcf_args = lambda wildcards, input: " ".join([f"-I {v}" for v in input.vcfs])
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000
    log:
        os.path.join(RESULTS_DIR, "logs", "03_merge_vcfs.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk MergeVcfs \
            {params.vcf_args} \
            -D {input.dict} \
            -O {output.vcf} 2> {log}
        """

# =============================================================================
# VARIANT QUALITY DIAGNOSTICS
# 
# Extract quality scores from raw variants for diagnostic plotting.
# =============================================================================

rule extract_snp_table:
    """Extract quality metrics from raw SNPs."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        table = os.path.join(RESULTS_DIR, "03_variants", "diagnostics", "raw_SNPs.table")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_diagnostics", "extract_snp_table.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.table})
        mkdir -p $(dirname {log})
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.table}.tmp.vcf 2> {log}
        gatk VariantsToTable \
            -V {output.table}.tmp.vcf \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
            --show-filtered \
            -O {output.table} 2>> {log}
        rm -f {output.table}.tmp.vcf {output.table}.tmp.vcf.idx
        """

rule extract_indel_table:
    """Extract quality metrics from raw indels."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        table = os.path.join(RESULTS_DIR, "03_variants", "diagnostics", "raw_INDELs.table")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_diagnostics", "extract_indel_table.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.table})
        mkdir -p $(dirname {log})
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {output.table}.tmp.vcf 2> {log}
        gatk VariantsToTable \
            -V {output.table}.tmp.vcf \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
            --show-filtered \
            -O {output.table} 2>> {log}
        rm -f {output.table}.tmp.vcf {output.table}.tmp.vcf.idx
        """
# =============================================================================
# 03b - VARIANT FILTRATION
#
# Workflow:
# 1. SelectVariants to split raw.vcf into SNPs and INDELs
# 2. VariantFiltration with type-specific thresholds from config
# 3. MergeVcfs to combine filtered SNPs and INDELs (site-level filtering)
# 4. Extract per-sample genotype depth (DP) from site-filtered BCF
# 5. Auto-compute global DP thresholds (median of per-sample 5th/99th pctl)
# 6. Apply genotype-level DP filter and set failed genotypes to no-call
# 7. Final output: filtered.bcf
#
# Diagnostic outputs for notebook:
#   - genotype_dp/dp_percentiles.tsv: per-sample percentile summary
#   - genotype_dp/dp_histograms.tsv: per-sample depth distributions (binned)
# =============================================================================


rule select_snps:
    """Extract SNPs from raw VCF."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw_SNPs.vcf")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "select_snps.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf} 2> {log}
        """


rule select_indels:
    """Extract INDELs from raw VCF."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw_INDELs.vcf")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "select_indels.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {output.vcf} 2> {log}
        """


rule filter_snps:
    """Apply hard filters to SNPs using thresholds from config."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw_SNPs.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_SNPs.vcf")
    params:
        QD = config["gatk_snp_filters"]["QD"],
        DP_low = config["gatk_snp_filters"]["DP_low"],
        DP_high = config["gatk_snp_filters"]["DP_high"],
        MQ = config["gatk_snp_filters"]["MQ"],
        FS = config["gatk_snp_filters"]["FS"],
        SOR = config["gatk_snp_filters"]["SOR"],
        MQRankSum_low = config["gatk_snp_filters"]["MQRankSum_low"],
        MQRankSum_high = config["gatk_snp_filters"]["MQRankSum_high"],
        ReadPosRankSum_low = config["gatk_snp_filters"]["ReadPosRankSum_low"],
        ReadPosRankSum_high = config["gatk_snp_filters"]["ReadPosRankSum_high"]
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "filter_snps.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {log})
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name "QD_filter" --filter-expression "QD < {params.QD}" \
            --filter-name "DP_low_filter" --filter-expression "DP < {params.DP_low}" \
            --filter-name "DP_high_filter" --filter-expression "DP > {params.DP_high}" \
            --filter-name "MQ_filter" --filter-expression "MQ < {params.MQ}" \
            --filter-name "FS_filter" --filter-expression "FS > {params.FS}" \
            --filter-name "SOR_filter" --filter-expression "SOR > {params.SOR}" \
            --filter-name "MQRankSum_low_filter" --filter-expression "MQRankSum < {params.MQRankSum_low}" \
            --filter-name "MQRankSum_high_filter" --filter-expression "MQRankSum > {params.MQRankSum_high}" \
            --filter-name "ReadPosRankSum_low_filter" --filter-expression "ReadPosRankSum < {params.ReadPosRankSum_low}" \
            --filter-name "ReadPosRankSum_high_filter" --filter-expression "ReadPosRankSum > {params.ReadPosRankSum_high}" 2> {log}
        """


rule filter_indels:
    """Apply hard filters to INDELs using thresholds from config."""
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw_INDELs.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_INDELs.vcf")
    params:
        QD = config["gatk_indel_filters"]["QD"],
        DP_low = config["gatk_indel_filters"]["DP_low"],
        DP_high = config["gatk_indel_filters"]["DP_high"],
        MQ = config["gatk_indel_filters"]["MQ"],
        FS = config["gatk_indel_filters"]["FS"],
        SOR = config["gatk_indel_filters"]["SOR"],
        MQRankSum_low = config["gatk_indel_filters"]["MQRankSum_low"],
        MQRankSum_high = config["gatk_indel_filters"]["MQRankSum_high"],
        ReadPosRankSum_low = config["gatk_indel_filters"]["ReadPosRankSum_low"],
        ReadPosRankSum_high = config["gatk_indel_filters"]["ReadPosRankSum_high"]
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "filter_indels.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {log})
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name "QD_filter" --filter-expression "QD < {params.QD}" \
            --filter-name "DP_low_filter" --filter-expression "DP < {params.DP_low}" \
            --filter-name "DP_high_filter" --filter-expression "DP > {params.DP_high}" \
            --filter-name "MQ_filter" --filter-expression "MQ < {params.MQ}" \
            --filter-name "FS_filter" --filter-expression "FS > {params.FS}" \
            --filter-name "SOR_filter" --filter-expression "SOR > {params.SOR}" \
            --filter-name "MQRankSum_low_filter" --filter-expression "MQRankSum < {params.MQRankSum_low}" \
            --filter-name "MQRankSum_high_filter" --filter-expression "MQRankSum > {params.MQRankSum_high}" \
            --filter-name "ReadPosRankSum_low_filter" --filter-expression "ReadPosRankSum < {params.ReadPosRankSum_low}" \
            --filter-name "ReadPosRankSum_high_filter" --filter-expression "ReadPosRankSum > {params.ReadPosRankSum_high}" 2> {log}
        """


rule merge_filtered_variants:
    """
    Merge filtered SNPs and INDELs back together.
    Keep only sites that PASS all site-level filters.
    This is an intermediate BCF before genotype-level filtering.
    """
    input:
        snps = os.path.join(RESULTS_DIR, "03_variants", "filtered_SNPs.vcf"),
        indels = os.path.join(RESULTS_DIR, "03_variants", "filtered_INDELs.vcf"),
        dict = REF_PREFIX + ".dict"
    output:
        vcf = temp(os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf")),
        idx = temp(os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf.idx"))
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "merge_filtered.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17",
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        
        gatk MergeVcfs \
            -I {input.snps} \
            -I {input.indels} \
            -D {input.dict} \
            -O {output.vcf}.tmp.vcf 2> {log}
        
        bcftools view -f PASS {output.vcf}.tmp.vcf -Ov -o {output.vcf} 2>> {log}
        gatk IndexFeatureFile -I {output.vcf} 2>> {log}
        
        rm -f {output.vcf}.tmp.vcf {output.vcf}.tmp.vcf.idx
        """


# =============================================================================
# GENOTYPE-LEVEL DEPTH FILTERING
#
# Genotypes with DP outside the 5th-99th percentile range are set to no-call.
# Thresholds are auto-computed from the data: global lower = median of
# per-sample 5th percentiles, global upper = median of per-sample 99th pctl.
# =============================================================================


rule extract_genotype_dp:
    """
    Extract per-sample genotype depth (DP) from site-filtered VCF.
    Output is a tab-delimited table: one column per sample, one row per site.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf")
    output:
        dp_table = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_table.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "extract_dp.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.dp_table})
        mkdir -p $(dirname {log})

        # Header: sample names
        bcftools query -l {input.vcf} | tr '\\n' '\\t' | \
            sed 's/\\t$/\\n/' > {output.dp_table} 2> {log}

        # Body: per-sample DP values
        bcftools query -f '[%DP\\t]\\n' {input.vcf} | \
            sed 's/\\t$//' >> {output.dp_table} 2>> {log}
        """


rule compute_dp_thresholds:
    """
    Compute per-sample DP percentiles and auto-select global thresholds.
    Global lower = median of per-sample 5th percentiles (floored).
    Global upper = median of per-sample 99th percentiles (ceiled).

    Diagnostic outputs for notebook:
      - dp_percentiles.tsv: per-sample stats
      - dp_histograms.tsv: per-sample binned distributions
    Pipeline output:
      - dp_thresholds.tsv: DP_lower and DP_upper for the filter step
    """
    input:
        dp_table = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_table.tsv")
    output:
        percentiles = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_percentiles.tsv"),
        histograms = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_histograms.tsv"),
        thresholds = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_thresholds.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "compute_dp_thresholds.log")
    shell:
        """
        mkdir -p $(dirname {log})
        python scripts/compute_dp_thresholds.py \
            {input.dp_table} \
            {output.percentiles} \
            {output.histograms} \
            {output.thresholds} &> {log}
        """


rule apply_genotype_dp_filter:
    """
    Apply genotype depth filter using auto-computed thresholds.
    Reads DP_lower and DP_upper from the thresholds file.
    Marks genotypes outside the range, then sets them to no-call.
    Converts final result to BCF and indexes.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf"),
        thresholds = os.path.join(RESULTS_DIR, "03_variants", "genotype_dp", "dp_thresholds.tsv"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        bcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "apply_genotype_dp_filter.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17",
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})

        # Read thresholds from file
        DP_LOWER=$(awk '$1=="DP_lower" {{print $2}}' {input.thresholds})
        DP_UPPER=$(awk '$1=="DP_upper" {{print $2}}' {input.thresholds})

        echo "Applying genotype DP filter: DP < $DP_LOWER || DP > $DP_UPPER" > {log}

        # Step 1: Mark genotypes outside the DP range
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.bcf}.tmp1.vcf \
            --genotype-filter-expression "DP < $DP_LOWER || DP > $DP_UPPER" \
            --genotype-filter-name "DP_${{DP_LOWER}}-${{DP_UPPER}}" 2>> {log}

        # Step 2: Set filtered genotypes to no-call
        gatk SelectVariants \
            -R {input.ref} \
            -V {output.bcf}.tmp1.vcf \
            --set-filtered-gt-to-nocall \
            -O {output.bcf}.tmp2.vcf 2>> {log}

        # Step 3: Convert to BCF and index
        bcftools view {output.bcf}.tmp2.vcf -Ob -o {output.bcf} 2>> {log}
        bcftools index {output.bcf} 2>> {log}

        rm -f {output.bcf}.tmp1.vcf {output.bcf}.tmp1.vcf.idx \
              {output.bcf}.tmp2.vcf {output.bcf}.tmp2.vcf.idx
        """
