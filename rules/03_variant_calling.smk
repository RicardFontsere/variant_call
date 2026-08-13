# =============================================================================
# 03 - VARIANT CALLING
# 
# Workflow:
# 1. HaplotypeCaller per sample per interval -> GVCFs
# 2. GenomicsDBImport per interval -> GenomicsDB
# 3. GenotypeGVCFs per interval -> VCFs
# 4. SortVcf per interval -> sorted VCFs
# 5. MergeVcfs -> raw.vcf
# 6. SelectVariants (SNPs) + VariantFiltration -> filtered_sites.vcf
# 7. Genotype-level GQ filter -> no-call failing genotypes
# 8. Convert to BCF + index -> filtered.bcf
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
# 03b - VARIANT FILTRATION (SNPs only)
#
# Workflow:
# 1. Single site-level step straight off raw.vcf:
#    SelectVariants (SNPs only) -> VariantFiltration (config thresholds)
#    -> keep PASS  => filtered_sites.vcf
# 2. Extract per-sample genotype quality (GQ) from the site-filtered VCF
# 3. Apply the genotype-level GQ floor, set failed genotypes to no-call,
#    refresh INFO tags and drop dead sites
# 4. Final output: filtered.bcf
#
# INDELs are dropped at step 1, so there is no separate INDEL branch and no
# SNP/INDEL merge. The raw_INDELs.table diagnostic above is still produced --
# it is inspection-only and never feeds the filtering path.
#
# Diagnostic outputs for notebook:
#   - genotype_gq/gq_table.tsv: per-sample GQ per site (wide TSV)
# =============================================================================


rule filter_snps:
    """
    Select SNPs from the raw VCF and apply the hard site filters in one step.

    Emits filtered_sites.vcf so every downstream rule stays unchanged.
    Only sites passing every filter survive.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = temp(os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf")),
        idx = temp(os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf.idx"))
    params:
        QD = config["gatk_snp_filters"]["QD"],
        QUAL = config["gatk_snp_filters"].get("QUAL", 30.0),
        MQ = config["gatk_snp_filters"]["MQ"],
        FS = config["gatk_snp_filters"]["FS"],
        SOR = config["gatk_snp_filters"]["SOR"],
        MQRankSum_low = config["gatk_snp_filters"]["MQRankSum_low"],
        ReadPosRankSum_low = config["gatk_snp_filters"]["ReadPosRankSum_low"],
        # Opt-in max-depth filter. Collapsed paralogs -- two loci merged in the
        # assembly -- pile reads from both onto one site and fake heterozygosity
        # in every sample at once, which is the dominant false positive for the
        # male-het/female-hom test in 12_heterozygosity.
        #
        # Set gatk_snp_filters.DP_max to ~2x the MEDIAN INFO/DP in
        # diagnostics/raw_SNPs.table (median, not mean: the repeat tail drags
        # the mean up and would make the cut too permissive). INFO/DP is the
        # sum over all samples, so uneven per-sample coverage needs no
        # adjustment. Leave unset and no depth filter is applied.
        dp_max_arg = (
            '--filter-name "DP_filter" --filter-expression "DP > {}"'.format(
                config["gatk_snp_filters"]["DP_max"]
            )
            if config["gatk_snp_filters"].get("DP_max") is not None
            else ""
        )
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "filter_snps.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17",
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})

        # 1. SNPs only, straight off the raw joint-called VCF
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf}.snps.vcf 2> {log}

        # 2. Annotate the FILTER column with the config thresholds
        gatk VariantFiltration \
            -R {input.ref} \
            -V {output.vcf}.snps.vcf \
            -O {output.vcf}.marked.vcf \
            --filter-name "QD_filter" --filter-expression "QD < {params.QD}" \
            --filter-name "QUAL_filter" --filter-expression "QUAL < {params.QUAL}" \
            --filter-name "MQ_filter" --filter-expression "MQ < {params.MQ}" \
            --filter-name "FS_filter" --filter-expression "FS > {params.FS}" \
            --filter-name "SOR_filter" --filter-expression "SOR > {params.SOR}" \
            {params.dp_max_arg} \
            --filter-name "MQRankSum_low_filter" --filter-expression "MQRankSum < {params.MQRankSum_low}" \
            --filter-name "ReadPosRankSum_low_filter" --filter-expression "ReadPosRankSum < {params.ReadPosRankSum_low}" 2>> {log}

        # 3. Keep only sites that PASS every filter, then index
        bcftools view -f PASS {output.vcf}.marked.vcf -Ov -o {output.vcf} 2>> {log}
        gatk IndexFeatureFile -I {output.vcf} 2>> {log}

        rm -f {output.vcf}.snps.vcf {output.vcf}.snps.vcf.idx \
              {output.vcf}.marked.vcf {output.vcf}.marked.vcf.idx
        """


# =============================================================================
# GENOTYPE-LEVEL GQ EXTRACTION (for threshold inspection)
#
# Extract per-sample genotype quality (GQ) per site into a wide table for
# plotting in a notebook. Inspect the distribution, pick a GQ floor (20 lenient,
# 30 stringent), then apply it in the genotype-filter step.
# Missing genotypes (e.g. hom-ref with no GQ) appear as '.'.
# =============================================================================

rule extract_genotype_gq:
    """
    Extract per-sample GQ from the site-filtered VCF.
    Wide TSV: one column per sample, one row per site. Loadable with pandas.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf")
    output:
        gq_table = os.path.join(RESULTS_DIR, "03_variants", "genotype_gq", "gq_table.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "extract_gq.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.gq_table})
        mkdir -p $(dirname {log})

        # Header: sample names
        bcftools query -l {input.vcf} | tr '\\n' '\\t' | \
            sed 's/\\t$/\\n/' > {output.gq_table} 2> {log}

        # Body: per-sample GQ values
        bcftools query -f '[%GQ\\t]\\n' {input.vcf} | \
            sed 's/\\t$//' >> {output.gq_table} 2>> {log}
        """


rule write_batch_groups:
    output:
        os.path.join(RESULTS_DIR, "03_variants", "batch_groups.tsv")
    run:
        with open(output[0], "w") as f:
            for s, b in BATCH_MAP.items():
                f.write(f"{s}\t{b}\n")

rule apply_genotype_gq_filter:
    """
    Set genotypes with GQ below the config floor to no-call, then convert to BCF.
    GQ is depth-aware, so one floor is batch-fair (see notebook for the cut).
    After no-calling, INFO tags are refreshed and now-dead sites are dropped.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_sites.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        bcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        csi = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf.csi")
    params:
        gq_min = config["genotype_gq_min"],
        # Drop a site if more than this fraction of samples end up no-called.
        # Fraction of ALL samples, so it bites hard at small N: with 10 samples
        # 0.15 would keep only sites with at most one no-call.
        max_f_missing = config.get("max_f_missing", 0.5)
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration", "apply_genotype_gq_filter.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17",
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})

        # Step 1: Mark genotypes below the GQ floor
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.bcf}.tmp1.vcf \
            --genotype-filter-expression "GQ < {params.gq_min}" \
            --genotype-filter-name "GQ{params.gq_min}" 2> {log}

        # Step 2: Set filtered genotypes to no-call
        gatk SelectVariants \
            -R {input.ref} \
            -V {output.bcf}.tmp1.vcf \
            --set-filtered-gt-to-nocall \
            -O {output.bcf}.tmp2.vcf 2>> {log}

        # Step 3: Recompute INFO from post-no-call genotypes, drop dead sites, convert + index
        #   --trim-alt-alleles: remove ALTs no longer carried by any genotype
        #   +fill-tags: refresh AC/AN/AF/F_MISSING (they were stale after no-calling)
        #   -c1: keep only sites still carrying an ALT in some genotype. Counted
        #        from the GT column, NOT from INFO/AC: --trim-alt-alleles has
        #        already stripped the dead ALTs, so a fully no-called site has
        #        nothing left for fill-tags to count and emits no AC at all --
        #        an 'AC==0' test compares against a missing tag, evaluates false,
        #        and lets exactly the dead sites it targets through.
        #   drop over-missing sites
        bcftools view --trim-alt-alleles {output.bcf}.tmp2.vcf -Ou 2>> {log} \
          | bcftools +fill-tags -Ou -- -t AN,AC,AF,F_MISSING 2>> {log} \
          | bcftools view -c1 -e 'F_MISSING > {params.max_f_missing}' -Ob -o {output.bcf} 2>> {log}
        bcftools index {output.bcf} 2>> {log}

        rm -f {output.bcf}.tmp1.vcf {output.bcf}.tmp1.vcf.idx \
              {output.bcf}.tmp2.vcf {output.bcf}.tmp2.vcf.idx
        """

