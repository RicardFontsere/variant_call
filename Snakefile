import os
import glob

configfile: "master/config/config.yaml"

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
    if read_num not in ("1", "2"):
        raise ValueError(f"read_num must be '1' or '2', got {read_num}")
    
    # Try common patterns for read files
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
include: "rules/03_variant_calling.smk"

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
        # Final joint-called variants
        os.path.join(RESULTS_DIR, "03_variants", "raw.vcf"),
        os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf"),
        os.path.join(RESULTS_DIR, "03_variants", "filtered.vcf.idx")
