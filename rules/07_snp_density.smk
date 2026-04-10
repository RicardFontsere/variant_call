# =============================================================================
# 06 - SNP DENSITY: MALE-TO-FEMALE FOLD CHANGE
#
# Workflow:
# 1. Generate genomic windows from reference .fai
# 2. Split filtered VCF into per-sample VCFs (SNPs only)
# 3. Calculate SNP density per sample per window (bedtools coverage)
# 4. Normalize SNP counts per sample by genome-wide median
# 5. Aggregate normalized densities by sex
# 6. Calculate M:F log2 fold change per window
# =============================================================================


rule create_windows:
    """
    Generate fixed-size genomic windows from the reference .fai index.
    Used as -a input for bedtools coverage.
    """
    input:
        fai = REFERENCE + ".fai"
    output:
        bed = os.path.join(RESULTS_DIR, "07_snp_density", "windows.bed")
    params:
        window_size = config["window_size"]
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "create_windows.log")
    envmodules:
        "BEDTools/2.30.0-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        mkdir -p $(dirname {log})
        bedtools makewindows \
            -g {input.fai} \
            -w {params.window_size} > {output.bed} 2> {log}
        """


rule split_vcf_per_sample:
    """
    Extract a single sample from the joint-called filtered BCF.
    Keep only SNPs with at least one alt allele in this sample (-c1).
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    output:
        vcf = os.path.join(RESULTS_DIR, "07_snp_density", "per_sample", "{sample}.vcf.gz")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "split_{sample}.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        bcftools view \
            -c1 \
            -s {wildcards.sample} \
            --types snps \
            -Oz \
            -o {output.vcf} \
            {input.vcf} 2> {log}
        """


rule snp_density_per_sample:
    """
    Count SNPs per genomic window for a single sample using bedtools coverage.

    Output columns (7):
      1. chrom  2. start  3. end  4. SNP_count
      5. bases_covered  6. window_length  7. fraction_covered
    """
    input:
        windows = os.path.join(RESULTS_DIR, "07_snp_density", "windows.bed"),
        vcf = os.path.join(RESULTS_DIR, "07_snp_density", "per_sample", "{sample}.vcf.gz")
    output:
        density = os.path.join(RESULTS_DIR, "07_snp_density", "raw", "{sample}.txt")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=50000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "density_{sample}.log")
    envmodules:
        "BEDTools/2.30.0-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.density})
        mkdir -p $(dirname {log})
        bedtools coverage \
            -a {input.windows} \
            -b {input.vcf} > {output.density} 2> {log}
        """


rule normalize_snp_density:
    """
    Normalize per-sample SNP counts by dividing by the genome-wide median.
    Corrects for differences in sequencing depth across samples.

    Reads the 7-column bedtools coverage output.
    Replaces column 7 with: SNP_count (col4) / median(SNP_count).
    Column 7 is what extract_snp_density.py reads downstream.
    """
    input:
        density = os.path.join(RESULTS_DIR, "07_snp_density", "raw", "{sample}.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "07_snp_density", "normalized", "{sample}.txt")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    params:
        chr_prefix=config["chromosome_prefix"]
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "normalize_{sample}.log")
    envmodules:
        "R-bundle-CRAN/2025.10-foss-2025a"
    shell:
        """
        mkdir -p $(dirname {output.norm})
        mkdir -p $(dirname {log})
        Rscript --vanilla -e '
            data <- read.table("{input.density}", sep="\t", header=FALSE)
            data <- data[grepl("^{params.chr_prefix}", data$V1), ]
            data_sorted <- data[order(data$V1, data$V2), ]
            total <- sum(data_sorted$V4)
            if (total > 0) {{
                data_sorted$V7 <- (data_sorted$V4 / total) * 1e6
            }} else {{
                data_sorted$V7 <- 0
            }}
            write.table(data_sorted, "{output.norm}", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' >> {log} 2>&1
        """


rule aggregate_snp_density_males:
    """
    Aggregate normalized SNP density across all male samples.
    Calculates per-window: sum, average, and log2(average + 1).
    """
    input:
        normalized_files = expand(
            os.path.join(RESULTS_DIR, "07_snp_density", "normalized", "{sample}.txt"),
            sample=ML_SAMPLES
        )
    output:
        avg = os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_males.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "aggregate_males.log")
    shell:
        """
        mkdir -p $(dirname {output.avg})
        mkdir -p $(dirname {log})
        python {SCRIPTS_DIR}/extract_snp_density.py \
            {input.normalized_files} \
            {output.avg} >> {log} 2>&1
        """


rule aggregate_snp_density_females:
    """
    Aggregate normalized SNP density across all female samples.
    Calculates per-window: sum, average, and log2(average + 1).
    """
    input:
        normalized_files = expand(
            os.path.join(RESULTS_DIR, "07_snp_density", "normalized", "{sample}.txt"),
            sample=FL_SAMPLES
        )
    output:
        avg = os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_females.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "aggregate_females.log")
    shell:
        """
        mkdir -p $(dirname {output.avg})
        mkdir -p $(dirname {log})
        python {SCRIPTS_DIR}/extract_snp_density.py \
            {input.normalized_files} \
            {output.avg} >> {log} 2>&1
        """


rule snp_density_fold_change:
    """
    Calculate M:F log2 fold change of normalized SNP density per window.
    Output: log2((male_avg + 1) / (female_avg + 1)) per window.
    """
    input:
        females = os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_females.csv"),
        males = os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_males.csv")
    output:
        fc = os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_fc.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "07_snp_density", "fold_change.log")
    shell:
        """
        mkdir -p $(dirname {log})
        python {SCRIPTS_DIR}/snpdensity_fold_change.py \
            {input.females} \
            {input.males} \
            {output.fc} >> {log} 2>&1
        """
