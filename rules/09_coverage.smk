# =============================================================================
# 09 - COVERAGE: MALE-TO-FEMALE FOLD CHANGE
#
# Workflow:
# 1. Generate genomic windows from reference .fai index
# 2. Calculate windowed read coverage per sample (bedtools multicov)
# 3. Normalize coverage per sample (divide by genome-wide median)
# 4. Aggregate normalized coverage by sex
# 5. Calculate M:F log2 fold change per window
# =============================================================================


rule cov_create_windows:
    """
    Generate fixed-size genomic windows from the reference .fai index.
    Reuses the .fai already produced by samtools_faidx in 00_preprocess.
    """
    input:
        fai = REFERENCE + ".fai"
    output:
        bed = os.path.join(RESULTS_DIR, "09_coverage", "windows.bed")
    params:
        window_size = config["window_size"]
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "create_windows.log")
    envmodules:
        "BEDTools/2.31.0-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        mkdir -p $(dirname {log})
        bedtools makewindows \
            -g {input.fai} \
            -w {params.window_size} > {output.bed} 2> {log}
        """


rule cov_calculate_coverage:
    """
    Count reads per genomic window for a single sample using bedtools multicov.
    Uses deduplicated BAMs from the alignment step.

    Output columns (4):
      1. chrom  2. start  3. end  4. read_count
    """
    input:
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam"),
        bai = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam.csi"),
        windows = os.path.join(RESULTS_DIR, "09_coverage", "windows.bed")
    output:
        cov = os.path.join(RESULTS_DIR, "09_coverage", "raw", "{sample}_multicov.txt")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "{sample}_multicov.log")
    envmodules:
        "BEDTools/2.31.0-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.cov})
        mkdir -p $(dirname {log})
        bedtools multicov \
            -bams {input.bam} \
            -bed {input.windows} \
            > {output.cov} 2> {log}
        """


rule cov_normalize:
    """
    Normalize per-sample read counts by dividing by the genome-wide median.
    Corrects for differences in sequencing depth across samples.
    Only retains windows on chromosomes matching the configured prefix.

    Reads 4-column multicov output, filters to chromosomes, and appends
    a 5th column: read_count / median(read_count).
    """
    input:
        cov = os.path.join(RESULTS_DIR, "09_coverage", "raw", "{sample}_multicov.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "09_coverage", "normalized", "{sample}_multicov_normalized.txt")
    params:
        chr_prefix = CHROM_PREFIX
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "{sample}_normalize.log")
    envmodules:
        "R-bundle-CRAN/2025.10-foss-2025a"
    shell:
        """
        mkdir -p $(dirname {output.norm})
        mkdir -p $(dirname {log})
        Rscript --vanilla -e '
            data <- read.table("{input.cov}", sep="\t", header=FALSE)
            data <- data[grepl("^{params.chr_prefix}", data$V1), ]
            data_sorted <- data[order(data$V1, data$V2), ]
            med <- median(data_sorted$V4)
            if (med > 0) {{
                data_sorted$V5 <- data_sorted$V4 / med
            }} else {{
                data_sorted$V5 <- 0
            }}
            write.table(data_sorted, "{output.norm}", sep="\t",
                        row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' >> {log} 2>&1
        """


rule cov_aggregate_males:
    """
    Aggregate normalized coverage across all male samples.
    Calculates per-window: sum, average, and log2(average + 1).
    """
    input:
        normalized_files = expand(
            os.path.join(RESULTS_DIR, "09_coverage", "normalized", "{sample}_multicov_normalized.txt"),
            sample=ML_SAMPLES
        )
    output:
        avg = os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_males.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "aggregate_males.log")
    shell:
        """
        mkdir -p $(dirname {output.avg})
        mkdir -p $(dirname {log})
        python scripts/extract_coverage.py \
            {input.normalized_files} \
            {output.avg} >> {log} 2>&1
        """


rule cov_aggregate_females:
    """
    Aggregate normalized coverage across all female samples.
    Calculates per-window: sum, average, and log2(average + 1).
    """
    input:
        normalized_files = expand(
            os.path.join(RESULTS_DIR, "09_coverage", "normalized", "{sample}_multicov_normalized.txt"),
            sample=FL_SAMPLES
        )
    output:
        avg = os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_females.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "aggregate_females.log")
    shell:
        """
        mkdir -p $(dirname {output.avg})
        mkdir -p $(dirname {log})
        python scripts/extract_coverage.py \
            {input.normalized_files} \
            {output.avg} >> {log} 2>&1
        """


rule cov_fold_change:
    """
    Calculate M:F log2 fold change of normalized coverage per window.
    Output: log2((male_avg + 1) / (female_avg + 1)) per window.
    """
    input:
        females = os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_females.csv"),
        males = os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_males.csv")
    output:
        fc = os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_fc.csv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "09_coverage", "fold_change.log")
    shell:
        """
        mkdir -p $(dirname {log})
        python scripts/coverage_fold_change.py \
            {input.females} \
            {input.males} \
            {output.fc} >> {log} 2>&1
        """

