import os
import glob

# =============================================================================
# CONFIGURATION
# =============================================================================
READS_DIR = config["reads_dir"]
RESULTS_DIR = config["results_dir"]
REFERENCE = config["reference"]
REF_PREFIX = os.path.splitext(REFERENCE)[0]
CHROM_PREFIX = config["chromosome_prefix"]
CONTIG_PREFIX = config["contig_prefix"]
INTERVALS_DIR = os.path.join(os.path.dirname(REFERENCE), "intervals")

# Discover samples: subdirectories in READS_DIR
SAMPLES = sorted([d for d in os.listdir(READS_DIR) if os.path.isdir(os.path.join(READS_DIR, d))])

ML_SAMPLES = [s for s in SAMPLES if config["male_pattern"] in s]
FL_SAMPLES = [s for s in SAMPLES if config["female_pattern"] in s]

# Discover intervals from pre-generated interval files
def get_intervals(wildcards=None):
    """Get intervals after checkpoint completes."""
    checkpoint_output = checkpoints.generate_intervals.get().output[0]
    intervals = [os.path.basename(f).replace(".interval_list", "") 
                 for f in glob.glob(os.path.join(checkpoint_output, "*.interval_list"))]
    return intervals

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def get_read_file(sample, read_num):
    """
    Find the reads for each sample.
    
    Args:
        sample: Sample directory name
        read_num: Either "1" or "2" for R1/R2
        
    Returns:
        Path to the read file
    """
    patterns = [
        f"*_R{read_num}.fastq.gz",
        f"*_R{read_num}.fq.gz",
        f"*_{read_num}.fastq.gz",
        f"*_{read_num}.fq.gz"
    ]
    
    for pattern in patterns:
        files = glob.glob(os.path.join(READS_DIR, sample, pattern))
        if files:
            return files[0]
    
    raise FileNotFoundError(
        f"No read file found for sample '{sample}' with read number {read_num}. "
        f"Searched in {os.path.join(READS_DIR, sample)} for patterns: {patterns}"
    )

# =============================================================================
# INCLUDE RULE MODULES
# =============================================================================
include: "rules/00_preprocess.smk"
include: "rules/01_trimming.smk"
include: "rules/02_alignment.smk"
include: "rules/02b_coverage.smk"
include: "rules/03_variant_calling.smk"
include: "rules/04_global_pca.smk"
include: "rules/05_FST.smk"
include: "rules/06_Kmer_GWAS.smk"
include: "rules/07_snp_density.smk"
include: "rules/08_GWAS.smk"
include: "rules/09_coverage.smk"
include: "rules/10_kmer_bb.smk"

# =============================================================================
# TARGET RULE
# =============================================================================
rule all:
    input:
        # Reference indices
        multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        REFERENCE + ".fai",
        REF_PREFIX + ".dict",
        # Intervals (checkpoint)
        INTERVALS_DIR,
        # Trimmed reads
        expand(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1.fq.gz"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2.fq.gz"), sample=SAMPLES),
        # Aligned and processed BAMs
        expand(os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam.csi"), sample=SAMPLES),
        # Coverage stats
        os.path.join(RESULTS_DIR, "02_coverage", "all_samples_coverage_summary.tsv"),
        #VCF stats
        os.path.join(RESULTS_DIR, "03_variants", "diagnostics", "raw_SNPs.table"),
        os.path.join(RESULTS_DIR, "03_variants", "diagnostics", "raw_INDELs.table"),
        # Final joint-called variants
        os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        # PCA results
        #os.path.join(RESULTS_DIR, "04_pca", "pca_plot.done"),
        os.path.join(RESULTS_DIR, "04_pca", "pca_plot.png"),
        os.path.join(RESULTS_DIR, "04_pca", "pca_variance.png"),
        # FST results
        os.path.join(RESULTS_DIR, "05_FST", "males_vs_females.snp"),
        # KMER_GWAS
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "male_blast.out"),
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "female_blast.out"),
        #SNP density
        os.path.join(RESULTS_DIR, "07_snp_density", "results", "snpdensity_fc.csv"),
        #GWAS
        os.path.join(RESULTS_DIR, "08_GWAS", "gwas_gemma.assoc.txt"),
        #M:F coverage
        os.path.join(RESULTS_DIR, "09_coverage", "results", "coverage_fc.csv"),
        #Kmer_BB
        expand(os.path.join(RESULTS_DIR, "10_kmer_bb", "blast", "{sex}_contigs_blast.out"), sex=["male", "female"])