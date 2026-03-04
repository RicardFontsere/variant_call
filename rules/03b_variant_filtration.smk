# =============================================================================
# 03b - VARIANT FILTRATION (split by type, filter independently, merge)
#
# Workflow:
# 1. SelectVariants to split raw.vcf into SNPs and INDELs
# 2. VariantFiltration with type-specific thresholds from config
# 3. MergeVcfs to combine filtered SNPs and INDELs
# 4. Convert to BCF for downstream compatibility
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
    Keep only sites that PASS all filters.
    Convert to BCF.
    """
    input:
        snps = os.path.join(RESULTS_DIR, "03_variants", "filtered_SNPs.vcf"),
        indels = os.path.join(RESULTS_DIR, "03_variants", "filtered_INDELs.vcf"),
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
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
            -O {output.vcf}.tmp.vcf.gz 2> {log}
        
        bcftools view -f PASS {output.vcf}.tmp.vcf.gz -Ob -o {output.vcf} 2>> {log}
        bcftools index {output.vcf} 2>> {log}
        
        rm -f {output.vcf}.tmp.vcf.gz {output.vcf}.tmp.vcf.gz.tbi
        """
