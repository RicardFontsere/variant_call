#!/usr/bin/env bash
# =============================================================================
# ONE-OFF: re-run all downstream metrics with one sample dropped
#
# Not part of the pipeline. Nothing in rules/ or Snakefile is aware of this
# script -- it only stages a second results tree that the *unmodified*
# workflow can then be pointed at.
#
# What it does
# ------------
# 1. Builds a pruned reads dir (symlinks to every sample dir except the
#    excluded one). SAMPLES / ML_SAMPLES / FL_SAMPLES in the Snakefile are
#    derived from `os.listdir(READS_DIR)`, so dropping the symlink is all it
#    takes for the sample to disappear from every rule at once.
# 2. Builds a new results dir and symlinks in the per-sample artefacts that
#    do NOT change when a *different* sample is removed:
#       02_aligned/*.dedup.bam(.csi)      alignment
#       02_coverage/*.coverage*           genomecov histograms
#       06_kmer/<sample>/*                KMC counts + kmers_with_strand
#       09_coverage/raw|normalized/*      per-sample windowed read counts
#    These are the expensive steps and they are all pure functions of a
#    single sample, so reusing them is exact, not an approximation.
# 3. Writes results/03_variants/filtered.bcf by subsetting the existing
#    joint-called BCF, reproducing the cohort-level parts of
#    apply_genotype_gq_filter for N-1 samples (see the comment on that step
#    below for what is and is not recomputed).
# 4. Writes a config.yaml pointing at the two new directories and prints the
#    snakemake command to run.
#
# Everything cohort-dependent (PCA, FST, k-mer table + association, SNP
# density aggregation, GEMMA, coverage fold change, heterozygosity) is left
# for snakemake to rebuild from the new BCF and the reduced sample list.
#
# Usage
# -----
#   oneoff/exclude_sample.sh \
#       --sample        ML_07 \
#       --config        master/config/config_RT.yaml \
#       --out-reads     /scratch/.../reads_noML07 \
#       --out-results   /scratch/.../variantcalling_noML07 \
#       [--out-config   master/config/config_noML07.yaml] \
#       [--force]
#
# reads_dir and results_dir of the *original* run are read from --config.
# =============================================================================

set -euo pipefail

SAMPLE=""
CONFIG=""
NEW_READS=""
NEW_RESULTS=""
NEW_CONFIG=""
FORCE=0

usage() { sed -n '2,45p' "$0"; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)      SAMPLE="$2";      shift 2 ;;
        --config)      CONFIG="$2";      shift 2 ;;
        --out-reads)   NEW_READS="$2";   shift 2 ;;
        --out-results) NEW_RESULTS="$2"; shift 2 ;;
        --out-config)  NEW_CONFIG="$2";  shift 2 ;;
        --force)       FORCE=1;          shift ;;
        -h|--help)     usage ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$SAMPLE"      ]] || { echo "--sample is required" >&2; exit 1; }
[[ -n "$CONFIG"      ]] || { echo "--config is required" >&2; exit 1; }
[[ -n "$NEW_READS"   ]] || { echo "--out-reads is required" >&2; exit 1; }
[[ -n "$NEW_RESULTS" ]] || { echo "--out-results is required" >&2; exit 1; }
[[ -f "$CONFIG"      ]] || { echo "config not found: $CONFIG" >&2; exit 1; }

NEW_CONFIG="${NEW_CONFIG:-$(dirname "$CONFIG")/$(basename "${CONFIG%.yaml}")_no${SAMPLE}.yaml}"

# ---------------------------------------------------------------------------
# Read the original paths out of the config. Plain sed rather than a YAML
# parser so this has no dependency beyond coreutils + bcftools.
# ---------------------------------------------------------------------------
yaml_get() {
    sed -nE "s/^$1:[[:space:]]*\"?([^\"#]*[^\"# ])\"?[[:space:]]*(#.*)?$/\1/p" "$CONFIG" | head -1
}

OLD_READS="$(yaml_get reads_dir)"
OLD_RESULTS="$(yaml_get results_dir)"
MAX_F_MISSING="$(yaml_get max_f_missing)"
MAX_F_MISSING="${MAX_F_MISSING:-0.5}"   # matches config.get("max_f_missing", 0.5)

[[ -d "$OLD_READS"   ]] || { echo "reads_dir from config not found: '$OLD_READS'" >&2; exit 1; }
[[ -d "$OLD_RESULTS" ]] || { echo "results_dir from config not found: '$OLD_RESULTS'" >&2; exit 1; }
[[ -e "$OLD_READS/$SAMPLE" ]] || { echo "sample '$SAMPLE' not found in $OLD_READS" >&2; exit 1; }

OLD_BCF="$OLD_RESULTS/03_variants/filtered.bcf"
NEW_BCF="$NEW_RESULTS/03_variants/filtered.bcf"
[[ -f "$OLD_BCF" ]] || { echo "joint-called BCF not found: $OLD_BCF" >&2; exit 1; }

if [[ -e "$NEW_RESULTS" || -e "$NEW_READS" ]] && [[ "$FORCE" -eq 0 ]]; then
    echo "output dirs already exist; pass --force to add to them" >&2
    exit 1
fi

command -v bcftools >/dev/null || {
    echo "bcftools not on PATH -- module load BCFtools/1.18-GCC-12.3.0 first" >&2
    exit 1
}

echo "excluding    : $SAMPLE"
echo "reads   from : $OLD_READS  ->  $NEW_READS"
echo "results from : $OLD_RESULTS  ->  $NEW_RESULTS"
echo "max_f_missing: $MAX_F_MISSING"
echo

# ---------------------------------------------------------------------------
# 1. Pruned reads dir
#
# Link to the fully resolved target so get_batch_map() still recovers the
# sequencing run from the parent of the realpath, exactly as in the original
# run.
# ---------------------------------------------------------------------------
mkdir -p "$NEW_READS"
KEPT=()
for d in "$OLD_READS"/*; do
    s="$(basename "$d")"
    [[ -d "$d" ]] || continue
    [[ "$s" == "$SAMPLE" ]] && continue
    ln -sfn "$(readlink -f "$d")" "$NEW_READS/$s"
    KEPT+=("$s")
done
echo "kept ${#KEPT[@]} samples in $NEW_READS"

# ---------------------------------------------------------------------------
# 2. Symlink the per-sample artefacts worth reusing
# ---------------------------------------------------------------------------
link_if_exists() {   # link_if_exists <src> <dest>
    [[ -e "$1" ]] || return 0
    mkdir -p "$(dirname "$2")"
    ln -sfn "$(readlink -f "$1")" "$2"
}

for s in "${KEPT[@]}"; do
    # alignment
    link_if_exists "$OLD_RESULTS/02_aligned/$s.dedup.bam"     "$NEW_RESULTS/02_aligned/$s.dedup.bam"
    link_if_exists "$OLD_RESULTS/02_aligned/$s.dedup.bam.csi" "$NEW_RESULTS/02_aligned/$s.dedup.bam.csi"

    # genomecov histogram + per-sample summary (02b)
    link_if_exists "$OLD_RESULTS/02_coverage/$s.coverage.txt"         "$NEW_RESULTS/02_coverage/$s.coverage.txt"
    link_if_exists "$OLD_RESULTS/02_coverage/$s.coverage_summary.tsv" "$NEW_RESULTS/02_coverage/$s.coverage_summary.tsv"

    # k-mer counting: the step this whole script exists to avoid repeating.
    # kmers_with_strand is the only file rule combine_kmers consumes; the KMC
    # databases are linked too so the branch stays self-consistent if anything
    # upstream is ever asked for explicitly.
    for f in "$OLD_RESULTS/06_kmer/$s"/output_kmc_* "$OLD_RESULTS/06_kmer/$s/kmers_with_strand"; do
        [[ -e "$f" ]] || continue
        link_if_exists "$f" "$NEW_RESULTS/06_kmer/$s/$(basename "$f")"
    done

    # windowed read counts (09): pure function of one BAM, and its
    # normalisation is per-sample too (divide by that sample's own median)
    link_if_exists "$OLD_RESULTS/09_coverage/raw/${s}_multicov.txt" \
                   "$NEW_RESULTS/09_coverage/raw/${s}_multicov.txt"
    link_if_exists "$OLD_RESULTS/09_coverage/normalized/${s}_multicov_normalized.txt" \
                   "$NEW_RESULTS/09_coverage/normalized/${s}_multicov_normalized.txt"
done
echo "linked reusable per-sample artefacts into $NEW_RESULTS"

# The window BEDs are derived from the reference .fai alone.
link_if_exists "$OLD_RESULTS/09_coverage/windows.bed"    "$NEW_RESULTS/09_coverage/windows.bed"
link_if_exists "$OLD_RESULTS/07_snp_density/windows.bed" "$NEW_RESULTS/07_snp_density/windows.bed"

# Deliberately NOT linked, because they are cohort-level and must be rebuilt:
#   02_coverage/all_samples_coverage_summary.tsv
#   03_variants/{all_samples,males,females}.txt
#   06_kmer/combined/**            (k-mer table, --mac/--maf depend on N)
#   07_snp_density/{per_sample,raw,normalized,results}
#   09_coverage/results/**
#   04_pca, 05_FST, 08_GWAS, 12_hetero

# ---------------------------------------------------------------------------
# 3. Subset the joint-called BCF
#
# Mirrors step 3 of rule apply_genotype_gq_filter (03_variant_calling.smk) so
# the cohort-level bookkeeping matches what the pipeline would have produced
# for N-1 samples:
#   --trim-alt-alleles  drop ALTs no longer carried by any retained genotype
#   +fill-tags          recompute AN/AC/AF/F_MISSING for the smaller cohort
#   -c1                 drop sites whose only ALT carrier was the excluded
#                       sample (otherwise they linger as monomorphic records
#                       and inflate SNP-density denominators / distort FST)
#   -e F_MISSING > x    re-apply the missingness ceiling against the new N
#
# NOT recomputed: the GATK site annotations (QD, QUAL, MQ, FS, SOR, the rank
# sums) still summarise reads from all N samples, because they were baked in
# at joint genotyping. Recomputing them would mean redoing
# GenomicsDBImport + GenotypeGVCFs, and the per-sample GVCFs are temp() so
# that path starts back at HaplotypeCaller. Dropping one sample out of N from
# an already-filtered callset is the normal post-hoc move -- just say so in
# the methods. The genotype-level GQ filter is per-genotype and so is
# unaffected.
# ---------------------------------------------------------------------------
mkdir -p "$NEW_RESULTS/03_variants"

bcftools query -l "$OLD_BCF" | grep -qx "$SAMPLE" || {
    echo "sample '$SAMPLE' is not in $OLD_BCF" >&2; exit 1; }

echo "subsetting $OLD_BCF ..."
bcftools view -s "^$SAMPLE" --trim-alt-alleles "$OLD_BCF" -Ou \
  | bcftools +fill-tags -Ou -- -t AN,AC,AF,F_MISSING \
  | bcftools view -c1 -e "F_MISSING > $MAX_F_MISSING" -Ob -o "$NEW_BCF"
bcftools index -f "$NEW_BCF"

echo "  samples: $(bcftools query -l "$OLD_BCF" | wc -l) -> $(bcftools query -l "$NEW_BCF" | wc -l)"
echo "  sites  : $(bcftools index -n "$OLD_BCF") -> $(bcftools index -n "$NEW_BCF")"

# ---------------------------------------------------------------------------
# 4. Config for the re-run
# ---------------------------------------------------------------------------
sed -E "s|^reads_dir:.*|reads_dir: \"$NEW_READS\"|; s|^results_dir:.*|results_dir: \"$NEW_RESULTS\"|" \
    "$CONFIG" > "$NEW_CONFIG"
grep -qE "^reads_dir: \"$NEW_READS\"$"     "$NEW_CONFIG" || { echo "failed to rewrite reads_dir" >&2;   exit 1; }
grep -qE "^results_dir: \"$NEW_RESULTS\"$" "$NEW_CONFIG" || { echo "failed to rewrite results_dir" >&2; exit 1; }
echo "wrote $NEW_CONFIG"

# ---------------------------------------------------------------------------
# 5. What to run next
# ---------------------------------------------------------------------------
cat <<EOF

Done staging. Now dry-run the workflow against the new config and check that
the only jobs listed are the cohort-level ones -- no fastp, align,
mark_duplicates, kmer_analysis or add_strand_information should appear:

  snakemake --profile master/profile --configfile $NEW_CONFIG -n \\
    $NEW_RESULTS/04_pca/pca_plot.png \\
    $NEW_RESULTS/04_pca/pca_variance.png \\
    $NEW_RESULTS/05_FST/males_vs_females.snp \\
    $NEW_RESULTS/06_kmer/combined/blast/male_blast.out \\
    $NEW_RESULTS/06_kmer/combined/blast/female_blast.out \\
    $NEW_RESULTS/07_snp_density/results/snpdensity_fc.csv \\
    $NEW_RESULTS/08_GWAS/gwas_gemma.assoc.txt \\
    $NEW_RESULTS/09_coverage/results/coverage_fc.csv \\
    $NEW_RESULTS/12_hetero/heterozygosity/heterozygosity.done \\
    $NEW_RESULTS/02_coverage/all_samples_coverage_summary.tsv

Targets are listed explicitly rather than using 'rule all' so the re-run does
not ask for raw.vcf and the 03_variants diagnostics, which would drag the
whole calling stage back in. Drop the -n to execute.

If the dry run wants to rebuild something that was linked, snakemake decided
it was stale rather than missing; --touch on that specific file fixes it
without recomputing.
EOF
