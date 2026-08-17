#!/usr/bin/env Rscript
# =============================================================================
# edgeR differential expression analysis.
#
# Original: EdgeR_script_EU_2026.R (Ezgi Unal, VUB). Kept as it was, except for
# what the header of each fix explains: the annotation file is now an argument
# instead of a hard-coded name, and the count matrix, the annotation and the
# contrast matrix are checked against the design before the model is fitted.
#
# Usage:
#   Rscript edger_de.R <analysis> <FDR> <filtering> <heatmaps> <datapath> <outpath> <annotation>
#
# <datapath> holds <analysis>_count.txt, <analysis>_design.txt and
# <analysis>_matrix.txt; results are written to <outpath>.
# =============================================================================

rm(list = ls(all.names = TRUE))
gc()

library(edgeR)
library(limma)
library(gplots)
library(dynamicTreeCut)


#General colorblind-friendly sex-based colors
color_male   <- "#009E73"
color_female <- "#E69F00"


get_group_colors <- function(group) {
  # 1) Phenotype-based patterns take priority
  # "XX_F", "XX_FB", "Am_46_F"   -> female
  # "XX_M", "XX_MB", "XYp_MB", "Am_46_M" -> male

  if (grepl("_F", group)) {
    return(color_female)
  } else if (grepl("_M", group)) {
    return(color_male)
  }

  # 2) If no explicit F/M code in the name, fall back to genotype-based rules
  if (grepl("XX", group)) {
    # XX without an explicit 'M' or 'F' tag -> default to female
    return(color_female)
  } else if (grepl("XYp|XY|XZ", group)) {
    # Y- or Z-bearing genotypes -> default to male
    return(color_male)
  }

  # 3) Anything else gets a neutral color
  return("darkgrey")
}


#####################################################################

#Setup input output
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: edger_de.R <analysis> <FDR> <filtering> <heatmaps yes|no> ",
       "<datapath> <outpath> <annotation>")
}

sub_analyse <- args[1]
FDR2use     <- as.numeric(args[2])
filtering   <- args[3]
runheat     <- args[4]
datapath    <- args[5]
outpath     <- args[6]
annot_file  <- args[7]

if (!filtering %in% c("fdr", "fdr_noss")) {
  stop("Unknown filtering '", filtering, "', expected 'fdr' or 'fdr_noss'")
}

if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)


#START!
#count file as nameofthesample_count.txt
count <- read.delim(file.path(datapath, paste(sub_analyse, '_count.txt', sep="")),
                    header=T, row.names=1, check.names=FALSE)
count <- round(count, digits=0)
annotation <- read.delim(annot_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
design <- read.delim(file.path(datapath, paste(sub_analyse, '_design.txt', sep="")), header=T)


#Fix chr name issues in the annotation if exist
unique(annotation$chr)
annotation$chr <- trimws(annotation$chr)
annotation$chr <- gsub("\\s+", "", annotation$chr)
unique(annotation$chr)


#Put the samples in the order of the design: the log-CPM table and the heatmaps
#below name the columns from design$sample, so a different order in the count
#file would silently swap samples.
design$sample <- as.character(design$sample)
missing <- setdiff(design$sample, colnames(count))
if (length(missing)) stop("Samples of the design missing from the count file: ",
                          paste(missing, collapse=", "))
count <- count[, design$sample, drop=FALSE]

#Same for the annotation: edgeR keeps it row by row next to the counts
annotation <- annotation[match(rownames(count), annotation$gid), ]
if (anyNA(annotation$gid)) stop("Transcripts of the count file are missing from the annotation")


#design check
model.formula <- as.formula("~0+group")
dmat <- model.matrix(model.formula, data=as.data.frame(design))
colnames(dmat) <- sub("^group", "", colnames(dmat))
dgl <- DGEList(counts=count, group=design$group, genes=annotation)

cat(paste("All transcripts:", nrow(dgl)), "\n",
    paste("Total samples:", ncol(dgl)), "\n",
    paste("Total reads:", sum(dgl$counts)), "\n",
    paste("Groups:", paste(unique(dgl$samples$group), collapse = ", ")))


## Filtering in EdgeR
filter_file <- file.path(outpath, paste0(sub_analyse, "_filter_file.txt"))
if (!file.exists(filter_file)) { file.create(filter_file) }

#1. Filter by average log CPM > 0
dgl <- dgl[aveLogCPM(dgl) > 0,]

#2. Filter genes expressed in at least 2 samples with CPM >= 1
dgl <- dgl[rowSums(cpm(dgl) >= 1) >= 2,]

#Log the Filtering Results
write("Filtering Summary", file = filter_file, append = TRUE)
write(paste("Number of Genes After Filtering:", nrow(dgl)), file = filter_file, append = TRUE)

#Log the summary of CPM values per gene
cpm_summary <- summary(rowSums(cpm(dgl)) / ncol(dgl))
write("CPM Summary (Normalized by Sample Size):", file = filter_file, append = TRUE)
write(capture.output(cpm_summary), file = filter_file, append = TRUE)

#Log the summary of aveLogCPM values
write("Summary of aveLogCPM After Filtering:", file = filter_file, append = TRUE)
write(capture.output(summary(aveLogCPM(dgl))), file = filter_file, append = TRUE)

#Log the summary of total counts per gene
write("Summary of Total Counts Per Gene After Filtering:", file = filter_file, append = TRUE)
write(capture.output(summary(rowSums(dgl$counts))), file = filter_file, append = TRUE)


# Estimate data normalisation factors and dispersion
xcpm <- mglmOneGroup(dgl$counts)
dgl <- calcNormFactors(dgl, method="TMM")
dgl <- estimateDisp(dgl, dmat, robust=TRUE)
write.table(dgl$samples, file=file.path(outpath, paste(sub_analyse,"_libsize.txt", sep="")), quote=F, row.names=F, sep="\t")


###Mean-Difference (MD) Plot (MA plot)

# MD plot for all
pdf(file.path(outpath,paste('MDplot_all_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
for (k in 1:ncol(dgl)) { plotMD(dgl, column=k)
  abline(h=0, col="red", lty=2, lwd=2) }
dev.off()


### MDS Multidimensional Scaling Plot

#Get unique group levels
group_levels <- unique(design$group)

#Define custom shapes with enough unique values
pchs <- c(17, 16, 15, 3, 4, 8, 10, 18, 7, 6)[1:length(group_levels)]

#Map each sample to the correct shape based on its group
sample_shapes <- pchs[match(design$group, group_levels)]

# Assign custom colors to each sample based on its group
sample_colors <- sapply(design$group, get_group_colors)

# Define output
pdf(file.path(outpath, paste('MDS_', sub_analyse, '.pdf', sep = "")), width = 8, height = 8)
par(mar = c(5, 5, 4, 3))

# Set up the dataset
y <- dgl
colnames(y) <- paste(colnames(y), sep = "\n")

# ---- 1. MDS Plot (logFC) ----
plotMDS(y, method = "logFC",pch = sample_shapes, col = sample_colors,
        cex = 2.0, lwd = 10.0, main = paste(sub_analyse, "MDS plot - logFC"),
        cex.main = 1.8, cex.lab = 1.3)

# Add legend
legend('topright', inset = 0.05, legend = group_levels,
       pch = pchs, col = sapply(group_levels, get_group_colors),
       title = "Groups", cex = 0.8)


# ---- 2. Simplified MDS (logFC) ----
plotMDS(y,
        method = "logFC", cex = 0.8,
        col = sample_colors, main = paste(sub_analyse, "MDS plot - logFC"))

dev.off()



## Fit the Quasi-Likelihood Negative Binomial Model
fitres <- glmQLFit(dgl, dmat, robust=TRUE)

## Inspect results
head(fitres$AveLogCPM)
summary(fitres$df.prior)

###Dispersion plot
## Define the threshold for extreme dispersion
dispersion_threshold <- quantile(dgl$tagwise.dispersion, 0.95)

## Identify outlier genes
outlier_genes <- rownames(dgl$counts)[dgl$tagwise.dispersion > dispersion_threshold]
## Save outlier genes to a file
write.table(outlier_genes, file=file.path(outpath, paste(sub_analyse, "_5percent_outlier_genes.txt", sep="")),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

## Save dispersion plot as PDF
pdf(file.path(outpath, paste('Dispersion_', sub_analyse, '.pdf', sep="")), width=10, height=8)
par(mar=c(5,5,4,3))

# EdgeR's Default BCV Plot
plotBCV(dgl, main="BCV Plot - EdgeR")

# EdgeR's Quasi-Likelihood Dispersion Plot
plotQLDisp(fitres, main="Quasi-Likelihood Dispersion Plot - EdgeR")

# Custom Dispersion Plot with Outliers Highlighted
ylim_range <- c(0, max(dgl$tagwise.dispersion, na.rm = TRUE))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.7, col="black",
     xlab="log2CPM", ylab="Dispersion",
     main=paste(sub_analyse, "Custom Dispersion Plot with Outliers"),
     ylim=ylim_range)

points(xcpm[dgl$tagwise.dispersion > dispersion_threshold],
       dgl$tagwise.dispersion[dgl$tagwise.dispersion > dispersion_threshold],
       col="red", pch=16, cex=0.7)

# Overlay trended dispersion if available
if (!is.null(dgl$trended.dispersion)) {
  points(xcpm, dgl$trended.dispersion,
         pch=16, cex=0.7, col="blue")
}

# Add horizontal line for common dispersion
abline(h=dgl$common.dispersion, col="green", lwd=2)

legend("topright", c("Common Dispersion", "Trended Dispersion", "Tagwise Dispersion", "Outliers 5%"),
       pch=16, col=c("green", "blue", "black", "red"),
       title="Dispersion Types")

dev.off()


#########################################################################################################


##ANALYSIS AND CONTRAST

# Read contrast matrix
x <- read.delim(file.path(datapath, paste(sub_analyse, "_matrix.txt", sep="")), sep="\t", header=T)
sortedX <- data.frame(x[order(x$model_coefficients, decreasing=F),])
cmat <- as.matrix(sortedX[,-1, drop=FALSE])
rownames(cmat) <- as.character(sortedX[,1])

# Sorting the rows of the contrast matrix by group name puts them in the order
# of the design matrix columns (model.matrix sorts the factor levels the same
# way). Anything else means the contrasts would be applied to the wrong groups.
if (!identical(rownames(cmat), colnames(dmat))) {
  stop("Groups of the contrast matrix (", paste(rownames(cmat), collapse=", "),
       ") do not match the groups of the design (", paste(colnames(dmat), collapse=", "), ")")
}

# Initialize result list
lrtres <- list()

# Choose appropriate DE test method
for(k in 1:ncol(cmat)) lrtres[[k]] <- glmQLFTest(fitres, contrast=cmat[,k])

# Store results
logFC <- NULL
PV <- NULL
FDR <- NULL

for(k in 1:ncol(cmat)) {
  PV <- cbind(PV, lrtres[[k]]$table[,"PValue"])
  FDR <- cbind(FDR, p.adjust(PV[,k], method="BH"))  # Benjamini-Hochberg
  logFC <- cbind(logFC, lrtres[[k]]$table[,"logFC"])}

# Extract logCPM values
xcpm <- lrtres[[1]]$table[,"logCPM"]

# Filter genes based on expression levels
allzeros <- which(rowSums(cpm(dgl)) < 3,)
allused <- which(rowSums(cpm(dgl)) >= 3,)

# Number of genes filtered out
length(allzeros)
# Number of genes retained
length(allused)

# Rename columns for clarity
cname <- colnames(cmat)
colnames(logFC) <- paste("logFC", cname, sep=".")
colnames(PV) <- paste("PV", cname, sep=".")
colnames(FDR) <- paste("FDR", cname, sep=".")
rownames(logFC) <- rownames(PV) <- rownames(FDR) <- rownames(fitres$coefficients)

# Save results
write.table(logFC, file=file.path(outpath, paste(sub_analyse, "_logFC.txt", sep="")), sep="\t", quote=F)
write.table(PV, file=file.path(outpath, paste(sub_analyse, "_PValue.txt", sep="")), sep="\t", quote=F)
write.table(FDR, file=file.path(outpath, paste(sub_analyse, "_FDR.txt", sep="")), sep="\t", quote=F)

nrow(logFC)
nrow(PV)
nrow(FDR)

###Comprehensive results table
# Identify genes with low counts
idxzeros <- allzeros
length(idxzeros)  # Number of low-expression genes

#Create Results Table
#Columns taken by name from the annotation: gene ID, name, product, chromosome,
#chromosome ID, start, end, full length and effective length.
annot_cols <- c("gid", "gname", "gproduct", "chr", "chr_id", "start", "end", "length", "eff_length")
missing_cols <- setdiff(annot_cols, colnames(dgl$genes))
if (length(missing_cols)) stop("Annotation columns missing: ", paste(missing_cols, collapse=", "))

restab <- data.frame(
  dgl$genes[, annot_cols],
  logCPM = xcpm,                         # Average log2 counts per million
  logFC,                                 # Log-fold change
  PV,                                    # Raw p-values
  FDR )                                  # Adjusted p-values (FDR)

sapply(restab, function(x) sum(is.na(x)))

# Remove Genes with Missing Values
restab <- restab[complete.cases(restab), ]
head(restab)
tail(restab)
nrow(restab)


## Save full DE results table
write.table(restab, file=file.path(outpath, paste(sub_analyse, "_DE_results.txt", sep="")),
            sep="\t", quote=F, row.names=F)

# Save id and name of analyzed genes
write.table(restab[, c("gid", "gname")],
            file = file.path(outpath, paste(sub_analyse, "_used_gene_names.txt", sep = "")),
            quote = FALSE, row.names = FALSE, sep = "\t")


##########################
## P-value Histogram Plot

# Define output file
pdf(file.path(outpath, paste('p_hist_', FDR2use, '_', sub_analyse, '.pdf', sep = "")), width = 8, height = 8)

# Determine the number of panels
npanel <- ncol(logFC)
np <- ceiling(sqrt(npanel))
if (np * (np - 1) >= npanel) {
  mfcol <- c(np - 1, np)
} else if ((np - 1) * (np + 1) >= npanel) {
  mfcol <- c(np + 1, np - 1)
} else {
  mfcol <- c(np, np)}

# Set plotting layout margins
par(mfcol = mfcol, mar = c(4, 4, 2, 1))

# Loop to generate histograms for each contrast
for (k in 1:npanel) {
  hist(PV[, k],
       breaks = 50,
       col = "steelblue",
       border = "white",
       xlab = "P-value",
       main = colnames(logFC)[k],
       xlim = c(0, 1),
       ylim = c(0, max(table(cut(PV[, k], breaks = 50))) * 1.1))}

dev.off()


## Within group pairwise scatter plot
# Extract raw counts and group names
wx <- dgl$counts
wg <- as.character(dgl$samples$group)
wug <- unique(wg)

# Assign colors consistently based on group endings
group_colors <- sapply(wg, get_group_colors)

# Generate pairwise scatter plots for each group
for (k in 1:length(wug)) {
  ix <- wg %in% wug[k]

  # Plot only if at least two samples are available
  if (sum(ix) > 1) {
    xmat <- log2(wx[, ix] + 1)
    pdf(file.path(outpath, paste('pairwise_raw_count_', wug[k], '_', sub_analyse, '.pdf', sep = "")), width = 8, height = 8)
    par(mfrow = c(1, 1))
    pairs(xmat, pch = 16, col = group_colors[ix], cex = 0.4,
          main = paste("Pairwise Scatter Plot:", wug[k]))

    dev.off()}}


##################################
## Build a summary table of DE counts per contrast and logFC threshold
colnames(logFC) <- paste("logFC", cname, sep=".")
contrast_names <- gsub("^logFC\\.", "", colnames(logFC))
summary_rows <- list()

for (logFC_use in c(3, 2, 1, 0)) {
  for (contrast in contrast_names) {

    fdr_col   <- paste0("FDR.", contrast)
    logfc_col <- paste0("logFC.", contrast)

    # Define DE logic
    fdr_vals   <- restab[[fdr_col]]
    logfc_vals <- restab[[logfc_col]]

    if (logFC_use == 0) {
      de_idx <- !is.na(fdr_vals) & (fdr_vals < FDR2use) & (logfc_vals != 0)
    } else {
      de_idx <- !is.na(fdr_vals) & (fdr_vals < FDR2use) & (abs(logfc_vals) >= logFC_use)}

    # --- write the per-contrast gene table
    de_genes <- restab[de_idx, ]
    write.table(
      de_genes,
      file = file.path(outpath, paste0("de_", contrast, "_FDR", FDR2use,
                                       "_logFC", logFC_use, "_", sub_analyse, ".txt")),
      sep = "\t", quote = FALSE, row.names = FALSE)

    # --- accumulate summary row
    n_total <- sum(!is.na(fdr_vals))
    n_de    <- sum(de_idx)
    n_up    <- sum(de_idx & (logfc_vals > 0))
    n_down  <- sum(de_idx & (logfc_vals < 0))

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      contrast  = contrast,
      logFC_cut = logFC_use,
      n_total   = n_total,
      n_DE      = n_de,
      n_up      = n_up,
      n_down    = n_down,
      prop_DE   = ifelse(n_total > 0, n_de / n_total, NA_real_))}}

de_summary_table <- do.call(rbind, summary_rows)
de_summary_table <- de_summary_table[order(de_summary_table$contrast, -de_summary_table$logFC_cut), ]
write.table(
  de_summary_table,
  file = file.path(outpath, paste0("Number_de_FDR", FDR2use, "_", sub_analyse, ".txt")),
  sep = "\t", quote = FALSE, row.names = FALSE)


###########################

# logFC plot using EdgeR's `plotMD()`


for (contrast in contrast_names) {

  # Extract groups from contrast name
  contrast_split <- unlist(strsplit(contrast, '\\.'))
  group1 <- contrast_split[1]
  group2 <- contrast_split[2]

  # Define column names in restab
  fdr_col <- paste0("FDR.", contrast)
  logfc_col <- paste0("logFC.", contrast)

  # Loop through logFC thresholds
  for (logFC_use in c(3, 2, 1, 0)) {

    key_name <- paste(contrast, logFC_use, sep = "_")

    # Define DE status using your logic
    fdr_vals <- restab[[fdr_col]]
    logfc_vals <- restab[[logfc_col]]

    if (logFC_use == 0) {
      is_de <- ifelse(!is.na(fdr_vals) & fdr_vals < FDR2use & logfc_vals != 0,
                      ifelse(logfc_vals > 0, paste("up", group1), paste("down", group2)),
                      "NS")
    } else {
      is_de <- ifelse(!is.na(fdr_vals) & fdr_vals < FDR2use & abs(logfc_vals) >= logFC_use,
                      ifelse(logfc_vals > 0, paste("up", group1), paste("down", group2)),
                      "NS")}

    # Make it a factor with consistent levels
    is_de <- factor(is_de, levels = c(paste("up", group1), paste("down", group2), "NS"))

    # Assign colors
    color_up <- get_group_colors(group1)
    color_down <- get_group_colors(group2)
    plot_colors <- c(color_up, color_down, "grey")

    # Output PDF
    pdf(file.path(outpath, paste0("FC-CPMplot_", sub_analyse, "_", contrast, "_logFC", logFC_use, ".pdf")),
        width = 8, height = 8)

    # Define y = logFC; x = CPM or proxy
    y_vals <- logfc_vals
    x_vals <- restab$logCPM

    # Plot manually
    plot(x_vals, y_vals,
         col = plot_colors[is_de],
         pch = 20,
         cex = 0.5,
         xlab = expression(log[2]~"CPM"),
         ylab = paste(group1, "-", group2),
         main = paste("LogFC plot:", group1, "vs", group2, "| logFC >", logFC_use))

    grid()
    abline(h = 0, col = "black", lty = 2)
    legend("bottomright", legend = levels(is_de), col = plot_colors, pch = 20)
    dev.off()}}

############################

## === Initialize DE lists ===

list_de <- list()
list_nonde <- list()

# Loop over logFC thresholds and contrasts
for (logFC_use in c(3, 2, 1, 0)) {
  for (contrast in contrast_names) {

    fdr_col <- paste0("FDR.", contrast)
    logfc_col <- paste0("logFC.", contrast)
    key_name <- paste(contrast, logFC_use, sep = "_")

    if (logFC_use == 0) {
      list_de[[key_name]] <- subset(restab,
                                    restab[[fdr_col]] < FDR2use &
                                      restab[[logfc_col]] != 0)
      list_nonde[[key_name]] <- subset(restab,
                                       restab[[fdr_col]] >= FDR2use |
                                         restab[[logfc_col]] == 0)
    } else {
      list_de[[key_name]] <- subset(restab,
                                    restab[[fdr_col]] < FDR2use &
                                      abs(restab[[logfc_col]]) >= logFC_use)
      list_nonde[[key_name]] <- subset(restab,
                                       restab[[fdr_col]] >= FDR2use |
                                         abs(restab[[logfc_col]]) < logFC_use)}}}


### Output subsets FDR05 for GO annotation
### 5% FDR, GO input files per contrast
for (k in 1:ncol(cmat)) {

  #Current contrast name
  contrast_name <- colnames(cmat)[k]

  #Column names in 'restab' for this contrast
  fdr_col   <- paste0("FDR.",   contrast_name)
  logfc_col <- paste0("logFC.", contrast_name)

  #Split contrast into the two groups
  group1 <- strsplit(contrast_name, "\\.")[[1]][1]
  group2 <- strsplit(contrast_name, "\\.")[[1]][2]

  #Genes significant at FDR < FDR2use
  FDR05 <- subset(restab, restab[[fdr_col]] < FDR2use)

  #Genes *not* significant (FDR > FDR2use)
  FDR05_unbiased <- subset(restab, restab[[fdr_col]] >= FDR2use)

  #Among significant genes, those up in group1 (logFC > 1)
  UP <- subset(FDR05, FDR05[[logfc_col]] > 1)

  #Among significant genes, those up in group2 (logFC < -1)
  DOWN <- subset(FDR05, FDR05[[logfc_col]] < -1)

  #Write FDR-significant gene IDs (no logFC threshold) – main DE gene list
  write.table(FDR05$gid,file = file.path(outpath, paste0("gene_", FDR2use, "_", contrast_name, "_", sub_analyse, ".txt")),
    quote = FALSE, row.names = FALSE, sep = "\t")

  #Write *non-significant* genes with logFC and FDR (background/unbiased set)
  write.table(data.frame(
      gene  = FDR05_unbiased$gid,
      logFC = FDR05_unbiased[[logfc_col]],
      FDR   = FDR05_unbiased[[fdr_col]]),
    file = file.path(outpath, paste0("unbiased_", FDR2use, "_", contrast_name, "_", sub_analyse, ".txt")),
    quote = FALSE, sep = "\t")

  # Write significant UP genes (group1-biased) with logFC and FDR
  write.table(data.frame(
      gene  = UP$gid,
      logFC = UP[[logfc_col]],
      FDR   = UP[[fdr_col]]),
    file = file.path(outpath, paste0("UP_", FDR2use, "_", contrast_name, "_", sub_analyse, ".txt")),
    quote = FALSE, sep = "\t")

  #Write significant DOWN genes (group2-biased) with logFC and FDR
  write.table(
    data.frame(
      gene  = DOWN$gid,
      logFC = DOWN[[logfc_col]],
      FDR   = DOWN[[fdr_col]]),
    file = file.path(outpath,paste0("DOWN_", FDR2use, "_", contrast_name, "_", sub_analyse, ".txt")),
    quote = FALSE, sep = "\t")

  #Create GO output directory for this contrast
  gopath <- file.path(outpath, paste0("GO_", FDR2use, "_", contrast_name))
  dir.create(gopath, showWarnings = FALSE)

  # GO file: list of FDR-significant genes (no direction)
  write.table(data.frame(gene = FDR05$gid), file = file.path(gopath, "FDR_genes.txt"),
    quote = FALSE, row.names = FALSE, sep = "\t")

  #GO file: group1-biased (UP) genes with logFC
  write.table(data.frame(
      gene  = UP$gid,
      logFC = UP[[logfc_col]]),
    file = file.path(gopath, paste0(group1, "_biased.txt")),
    quote = FALSE, row.names = FALSE, sep = "\t")

  # GO file: group2-biased (DOWN) genes with logFC
  write.table(data.frame(
      gene  = DOWN$gid,
      logFC = DOWN[[logfc_col]]),
    file = file.path(gopath, paste0(group2, "_biased.txt")),
    quote = FALSE, row.names = FALSE, sep = "\t")

  # GO file: full universe with FDR and logFC for this contrast
  GO_pvalues <- data.frame(
    gene  = as.character(restab$gid),
    FDR   = restab[[fdr_col]],
    logFC = restab[[logfc_col]])
  write.table(
    GO_pvalues,file = file.path(gopath, "GO_pvalues.txt"),
    quote = FALSE, row.names = FALSE, sep = "\t")


  # GO file: only group1-biased genes (FDR < threshold and logFC > 1) keep their FDR;
  # everything else is set to 1
  GO_pvalues_UP <- data.frame(gene = as.character(restab$gid), FDR  = restab[[fdr_col]])
  GO_pvalues_UP$FDR[!(restab[[fdr_col]] < FDR2use & restab[[logfc_col]] > 1)] <- 1
  write.table(GO_pvalues_UP, file = file.path(gopath, "GO_pvalues_UP.txt"),
              quote = FALSE, row.names = FALSE, sep = "\t")

  # GO file: only group2-biased genes (FDR < threshold and logFC < -1) keep their FDR;
  # everything else is set to 1
  GO_pvalues_DOWN <- data.frame(gene = as.character(restab$gid), FDR  = restab[[fdr_col]])
  GO_pvalues_DOWN$FDR[!(restab[[fdr_col]] < FDR2use & restab[[logfc_col]] < -1)] <- 1
  write.table(GO_pvalues_DOWN, file = file.path(gopath, "GO_pvalues_DOWN.txt"),
    quote = FALSE, row.names = FALSE, sep = "\t")}

## Export ALL results for EVERY gene analyzed and moderated log-counts-per-million
nc <- cpm(dgl, prior.count=1, log=T)
nc2 <- data.frame(row.names(nc), nc)
colnames(nc2) <- c('ID', as.character(design$sample))
restab_logCPM = merge(restab, nc2, by.x="gid", by.y="ID", all=F )
write.table(restab_logCPM, file=file.path(outpath, paste('LogCPM_',FDR2use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

nrow(nc)
nrow(nc2)
nrow(restab_logCPM)



##HEATMAPS

# Define colors
cbred <- "red"
cbblue <- "blue"

# Heatmaps
if (runheat == 'yes') {
  colnames(restab_logCPM)
  tail(restab_logCPM)
  cmat <- data.frame(cmat)

  for(k in 1:ncol(cmat)) {
    if (ncol(cmat) == 1) {
      cmat_subset <- cmat
    } else {
      cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]] != 0)
    }

    design_subset <- design[design$group %in% row.names(cmat_subset),]
    length(as.character(design_subset$sample))
    rownames(restab_logCPM) <- restab_logCPM$gid
    DE_counts <- subset(restab_logCPM, restab_logCPM[[paste('FDR.', colnames(cmat)[k], sep="")]] < FDR2use)
    DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))

    # row.names(DE_counts_relevant) <- DE_counts$gname
    colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group)

    if (nrow(DE_counts_relevant) > 2) {

      pdf(file.path(outpath, paste('Heatmap_', FDR2use, '_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)

      heatmap.2(as.matrix(DE_counts_relevant),
                col=colorpanel(100, cbred, 'white', cbblue),
                scale="row",
                key=TRUE,
                keysize=1,
                density.info="density",
                trace="none",
                cexCol=0.5,
                cexRow=0.6,
                main=paste(sub_analyse, strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs',
                           strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use),
                srtCol=0,
                key.title=NA)

      cmethods <- c('average', 'ward.D')

      for(m in 1:length(cmethods)) {
        d <- as.matrix(DE_counts_relevant)
        d[is.na(d)] <- 0

        myheatcol <- colorpanel(100, cbred, 'white', cbblue)
        distmatrix <- as.dist(1 - cor(t(d), method="pearson", use="pairwise.complete.obs"))
        hr <- hclust(distmatrix, method=cmethods[m])
        mycl <- cutreeDynamic(hr, distM=as.matrix(distmatrix))
        clusterCols <- rainbow(length(unique(mycl)))

        # Fix: Assign correct cluster colors
        myClusterSideBar <- clusterCols[as.factor(mycl)]

        heatmap.2(d, col=myheatcol,
                  Rowv=reorder(as.dendrogram(hr), wts=mycl),
                  keysize=1,
                  scale="row",
                  density.info="density",
                  trace="none",
                  cexCol=0.5,
                  cexRow=0.6,
                  RowSideColors=myClusterSideBar,
                  main=paste(sub_analyse, cmethods[m], strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs',
                             strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use),
                  srtCol=45,
                  key.title=NA)

        # Fix: Use correct color mapping for legend
        legend(x=0.15, y=1.12, legend=unique(mycl), col=clusterCols, lty=1, lwd=3, cex=.5, title="clusters", xpd=TRUE)

        DE_counts[[cmethods[m]]] = as.factor(mycl)

        gopath_method <- file.path(outpath, paste('GO_', FDR2use, '_', colnames(cmat)[k], '/', cmethods[m], sep=""))
        dir.create(gopath_method, showWarnings = FALSE, recursive = TRUE)

        for(l in 1:length(levels(DE_counts[[cmethods[m]]]))) {
          FDR05 <- data.frame(gene=DE_counts$gid,
                              cluster=DE_counts[[cmethods[m]]],
                              FDR=DE_counts[[paste('FDR.', colnames(cmat)[k], sep="")]])

          # Fix: Prevent incorrect FDR updates
          FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1

          GO_pvalues <- data.frame(gene=as.character(restab$gid),
                                   FDR=restab[[paste('FDR.', colnames(cmat)[k], sep="")]],
                                   logFC=restab[[paste('logFC.', colnames(cmat)[k], sep="")]])

          merged1 <- merge(GO_pvalues, FDR05, by.x="gene", by.y="gene", all=TRUE)
          merged1$FDR.y[is.na(merged1$FDR.y)] <- as.character(merged1$FDR.x[is.na(merged1$FDR.y)])

          write.table(data.frame(gene=merged1[,1], FDR=merged1[,5]),
                      file=file.path(gopath_method, paste('cluster_', l, '_', cmethods[m], '.txt', sep="")),
                      quote=FALSE, row.names=FALSE, sep='\t')
        }
      }

      dev.off()

      write.table(DE_counts,
                  file=file.path(outpath, paste('clusters', FDR2use, '_', sub_analyse, '_', colnames(cmat)[k], '.txt', sep="")),
                  quote=FALSE, row.names=FALSE, sep='\t')
    }
  }
}
