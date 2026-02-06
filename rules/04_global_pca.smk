# =============================================================================
# 05 - PCA ANALYSIS
# 
# Workflow:
# 1. Linkage pruning with PLINK to identify independent SNPs
# 2. Extract pruned sites, create bed files, and run PCA
# 3. Visualize PCA results with R/ggplot2
# =============================================================================

rule plink_prune:
    """
    Perform linkage pruning to identify independent SNPs.
    Uses sliding window of 50 SNPs, step size 10, r² threshold 0.1.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    output:
        prune_in = os.path.join(RESULTS_DIR, "04_pca", "pruned.prune.in"),
        prune_out = os.path.join(RESULTS_DIR, "04_pca", "pruned.prune.out")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "04_pca", "pruned")
    resources:
        cpus_per_task=5,
        mem_mb_per_cpu=4000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "plink_prune.log")
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

rule plink_pca:
    """
    Extract pruned sites, create bed files, and perform PCA.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        prune_in = os.path.join(RESULTS_DIR, "04_pca", "pruned.prune.in")
    output:
        eigenvec = os.path.join(RESULTS_DIR, "04_pca", "pca.eigenvec"),
        eigenval = os.path.join(RESULTS_DIR, "04_pca", "pca.eigenval"),
        bed = os.path.join(RESULTS_DIR, "04_pca", "pca.bed"),
        bim = os.path.join(RESULTS_DIR, "04_pca", "pca.bim"),
        fam = os.path.join(RESULTS_DIR, "04_pca", "pca.fam")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "04_pca", "pca")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=10000,
        runtime=480
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "plink_pca.log")
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

rule pca_plot:
    """
    Visualize PCA results using R and ggplot2.
    Generates variance explained barplot and PC1 vs PC2 scatter plot.
    """
    input:
        eigenvec = os.path.join(RESULTS_DIR, "04_pca", "pca.eigenvec"),
        eigenval = os.path.join(RESULTS_DIR, "04_pca", "pca.eigenval")
    output:
        variance_plot = os.path.join(RESULTS_DIR, "04_pca", "pca_variance.png"),
        pca_plot = os.path.join(RESULTS_DIR, "04_pca", "pca_plot.png")
    params:
        male_pattern = config["male_pattern"],
        female_pattern = config["female_pattern"]
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "04_pca", "pca_plot.log")
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
            # Remove column (first column and second column show the same sampleID)
            pca <- pca[, -1]

            # Set column names
            names(pca)[1] <- "ind"
            names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1)) #Principal Component numbers

            # Add sex based on sample naming patterns
            sex <- rep(NA, length(pca$ind)) #Start vector
            sex[grep("{params.male_pattern}", pca$ind)] <- "male"
            sex[grep("{params.female_pattern}", pca$ind)] <- "female"

            pca <- as_tibble(data.frame(pca, sex)) 

            # Calculate percentage variance explained
            pve <- data.frame(PC = 1:min(6, length(eigenval)), 
                              pve = eigenval[1:min(6, length(eigenval))] / sum(eigenval) * 100)

            # Plot 1: Variance explained barplot
            a <- ggplot(pve, aes(PC, pve)) + 
                geom_bar(stat = "identity", fill = "steelblue") +
                ylab("Percentage variance explained") + 
                xlab("Principal Component") +
                theme_light() +
                ggtitle("PCA Variance Explained");
            ggsave("{output.variance_plot}", plot = a, width = 8, height = 6, dpi = 500);
            
            # Plot 2: PCA scatter plot (PC1 vs PC2)
            b <- ggplot(pca, aes(PC1, PC2, col = sex)) + 
                geom_point(size = 3) +
                scale_colour_manual(values = c("female" = "red", "male" = "blue"), na.value = "grey50") +
                coord_equal() + 
                theme_light() +
                xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%%)")) + 
                ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%%)")) +
                ggtitle("PCA: PC1 vs PC2");
            ggsave("{output.pca_plot}", plot = b, width = 8, height = 8, dpi = 500);
        ' >> {log} 2>&1
        """
