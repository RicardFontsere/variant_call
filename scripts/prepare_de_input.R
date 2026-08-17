#!/usr/bin/env Rscript
# =============================================================================
# Build the edgeR input files for one differential-expression analysis.
#
# Adapted from Making_counts_file_EU_2026.R (Ezgi Unal, VUB) for use inside the
# Snakemake pipeline: paths come from the command line and the samples of the
# analysis are taken from the design file instead of being hard-coded.
#
# Counts, TPM and annotation are all reduced to the transcripts they have in
# common and written in annotation order, so that edgeR can attach the
# annotation to the count matrix row by row.
#
# Usage:
#   Rscript prepare_de_input.R counts.txt tpm.txt eff_length.txt annotation.tsv \
#                              design.txt out_count.txt out_annotation.txt out_tpm.txt
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: prepare_de_input.R counts tpm eff_length annotation design ",
       "out_count out_annotation out_tpm")
}

count_file  <- args[1]
tpm_file    <- args[2]
efflen_file <- args[3]
annot_file  <- args[4]
design_file <- args[5]
out_count   <- args[6]
out_annot   <- args[7]
out_tpm     <- args[8]

# check.names = FALSE: sample names are matched against the design file, they
# must not be mangled by make.names().
count  <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
tpm    <- read.delim(tpm_file, check.names = FALSE, stringsAsFactors = FALSE)
efflen <- read.delim(efflen_file, check.names = FALSE, stringsAsFactors = FALSE)
annot  <- read.delim(annot_file, check.names = FALSE, stringsAsFactors = FALSE)
design <- read.delim(design_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!"gid" %in% colnames(annot)) stop("The annotation needs a 'gid' column")
if (!"sample" %in% colnames(design)) stop("The design needs a 'sample' column")

samples <- as.character(design$sample)
missing <- setdiff(samples, colnames(count))
if (length(missing)) {
  stop("Samples of the design that were not quantified: ",
       paste(missing, collapse = ", "))
}

# Transcripts shared by the counts and the annotation, kept in annotation order
ids <- annot$gid[annot$gid %in% count$id]
if (!length(ids)) stop("Counts and annotation share no transcript ID")

annot  <- annot[match(ids, annot$gid), ]
count  <- count[match(ids, count$id), ]
tpm    <- tpm[match(ids, tpm$id), ]
efflen <- efflen[match(ids, efflen$id), ]

# The effective length is sample dependent, so the annotation carries the mean
# over the samples of this analysis. edgeR reports it next to the gene length.
annot$eff_length <- rowMeans(efflen[, samples, drop = FALSE])

cat("Transcripts in the annotation:", nrow(read.delim(annot_file)), "\n",
    "Transcripts kept (also quantified):", length(ids), "\n",
    "Samples in the analysis:", length(samples), "\n")

write.table(data.frame(id = ids, count[, samples, drop = FALSE],
                       check.names = FALSE),
            out_count, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(annot, out_annot, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(annot, tpm[, samples, drop = FALSE], check.names = FALSE),
            out_tpm, sep = "\t", quote = FALSE, row.names = FALSE)
