# =============================================================================
# 04 - PCA ANALYSIS (PER CHROMOSOME, COMBINED PLOT)
# 
# Workflow:
# 1. bgzip and index filtered VCF
# 2. Extract each chromosome from VCF
# 3. Linkage pruning with PLINK per chromosome
# 4. Run PCA per chromosome
# 5. Combine all results and plot in a SINGLE PCA plot
#    - One dot per chromosome per individual
#    - Color = individual
#    - Shape = sex
# =============================================================================


rule bgzip_vcf:
    """
    Compress filtered VCF with bgzip for indexed access.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf")
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.gz")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=2000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "bgzip_vcf.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bgzip -c -@ {resources.cpus_per_task} {input.vcf} > {output.vcf} 2> {log}
        """


rule tabix_vcf:
    """
    Index bgzipped VCF with tabix using CSI index (for large chromosomes >512Mb).
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.gz")
    output:
        idx = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.gz.csi")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "tabix_vcf.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools index -c {input.vcf} 2> {log}
        """


rule extract_chromosome:
    """
    Extract a single chromosome/interval from the filtered VCF.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.gz"),
        idx = os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.gz.csi")
    output:
        vcf = temp(os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "chr.vcf"))
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "by_chr", "{interval}", "extract_chr.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        bcftools view -r {wildcards.interval} {input.vcf} > {output.vcf} 2> {log}
        """


rule plink_prune_chr:
    """
    Perform linkage pruning per chromosome.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "chr.vcf")
    output:
        prune_in = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pruned.prune.in"),
        prune_out = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pruned.prune.out")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pruned")
    resources:
        cpus_per_task=5,
        mem_mb_per_cpu=4000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "by_chr", "{interval}", "plink_prune.log")
    envmodules:
        "PLINK/2.00a3.7-foss-2022a"
    shell:
        """
        mkdir -p $(dirname {params.out_prefix})
        mkdir -p $(dirname {log})
        plink --vcf {input.vcf} \
            --double-id \
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --indep-pairwise 50 10 0.1 \
            --out {params.out_prefix} 2> {log}
        """


rule plink_pca_chr:
    """
    Run PCA per chromosome.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "chr.vcf"),
        prune_in = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pruned.prune.in")
    output:
        eigenvec = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca.eigenvec"),
        eigenval = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca.eigenval")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=10000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "by_chr", "{interval}", "plink_pca.log")
    envmodules:
        "PLINK/2.00a3.7-foss-2022a"
    shell:
        """
        mkdir -p $(dirname {log})
        plink --vcf {input.vcf} \
            --double-id \
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --extract {input.prune_in} \
            --make-bed \
            --pca \
            --out {params.out_prefix} 2> {log}
        """


def get_pca_plots(wildcards):
    """Get PCA plot files for chromosomes only (excluding unplaced contigs)."""
    intervals = get_intervals()
    # Only include intervals that start with the chromosome prefix
    intervals = [i for i in intervals if i.startswith(CHROM_PREFIX)]
    return expand(os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca_plot.png"),
                  interval=intervals)


def get_pca_variance_plots(wildcards):
    """Get PCA variance plot files for chromosomes only (excluding unplaced contigs)."""
    intervals = get_intervals()
    # Only include intervals that start with the chromosome prefix
    intervals = [i for i in intervals if i.startswith(CHROM_PREFIX)]
    return expand(os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca_variance.png"),
                  interval=intervals)


rule pca_plot_chr:
    """
    Create PCA plot for each chromosome.
    Color = sex.
    """
    input:
        eigenvec = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca.eigenvec"),
        eigenval = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca.eigenval")
    output:
        pca_plot = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca_plot.png"),
        variance_plot = os.path.join(RESULTS_DIR, "04_pca", "by_chr", "{interval}", "pca_variance.png")
    params:
        male_pattern = config.get("male_pattern", "ML"),
        female_pattern = config.get("female_pattern", "FL")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "by_chr", "{interval}", "pca_plot.log")
    envmodules:
        "R-bundle-CRAN/2025.10-foss-2025a"
    shell:
        """
        mkdir -p $(dirname {log})
        Rscript --vanilla -e '
            library(tidyverse)

            # Read in data
            pca <- read_table("{input.eigenvec}", col_names = FALSE)
            eigenval <- scan("{input.eigenval}")

            # Process PCA data
            pca <- pca[, -1]
            names(pca)[1] <- "ind"
            names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

            # Add sex based on sample naming patterns
            pca$sex <- NA
            pca$sex[grep("{params.male_pattern}", pca$ind)] <- "male"
            pca$sex[grep("{params.female_pattern}", pca$ind)] <- "female"

            # Calculate percentage variance explained
            pve <- data.frame(PC = 1:min(6, length(eigenval)), 
                              pve = eigenval[1:min(6, length(eigenval))] / sum(eigenval) * 100)

            # Plot 1: Variance explained barplot
            a <- ggplot(pve, aes(PC, pve)) + 
                geom_bar(stat = "identity", fill = "steelblue") +
                ylab("Percentage variance explained") + 
                xlab("Principal Component") +
                theme_light() +
                ggtitle("PCA Variance Explained - {wildcards.interval}")
            ggsave("{output.variance_plot}", plot = a, width = 8, height = 6, dpi = 500)
            
            # Plot 2: PCA scatter plot (PC1 vs PC2)
            b <- ggplot(pca, aes(PC1, PC2, col = sex)) + 
                geom_point(size = 3) +
                scale_colour_manual(values = c("female" = "red", "male" = "blue"), na.value = "grey50") +
                coord_equal() + 
                theme_light() +
                xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%%)")) + 
                ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%%)")) +
                ggtitle("PCA: {wildcards.interval}")
            ggsave("{output.pca_plot}", plot = b, width = 8, height = 8, dpi = 500)
        ' >> {log} 2>&1
        """


rule aggregate_pca_plots:
    """
    Aggregate all per-chromosome PCA plots.
    This rule collects all plots after the checkpoint completes.
    """
    input:
        plots = get_pca_plots,
        variance_plots = get_pca_variance_plots
    output:
        done = touch(os.path.join(RESULTS_DIR, "04_pca", "all_plots.done"))
