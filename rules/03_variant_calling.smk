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
# VARIANT FILTRATION & SELECTION
# =============================================================================

rule variant_filtration:
    """
    Annotate filter tags on the raw VCF. This does NOT remove variants,
    it only marks them as PASS or with the relevant filter name.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_annotated.vcf")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "03_variant_filtration.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {log})
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name "QD_filter" --filter-expression "QD < 6.0" \
            --filter-name "DP_low_filter" --filter-expression "DP < 42" \
            --filter-name "DP_high_filter" --filter-expression "DP > 115" \
            --filter-name "FS_filter" --filter-expression "FS > 60.0" \
            --filter-name "MQ_filter" --filter-expression "MQ < 50.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
            --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5 || vc.hasAttribute('MQRankSum') == false" \
            --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0 || vc.hasAttribute('ReadPosRankSum') == false" 2> {log}
        """

rule select_variants:
    """
    Apply filters: keep only PASS biallelic SNPs.
    --exclude-filtered removes any variant not marked PASS.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_annotated.vcf"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_selected.vcf")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "03_select_variants.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {log})
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --restrict-alleles-to BIALLELIC \
            --select-type-to-include SNP \
            --exclude-filtered 2> {log}
        """

rule convert_to_bcf:
    """
    Convert filtered VCF to BCF.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered_selected.vcf")
    output:
        bcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "03_convert_to_bcf.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view {input.vcf} -Ob -o {output.bcf} 2> {log}
        """

rule index_bcf:
    """
    Index the final BCF file.
    """
    input:
        bcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    output:
        idx = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf.csi")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "03_index_bcf.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools index {input.bcf} 2> {log}
        """
