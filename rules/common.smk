# Configuration
READS_DIR = config["INPUT"]["reads"]
BASE_NAME = os.path.basename(READS_DIR.rstrip('/'))
RESULTS_DIR = os.path.join(config["OUTPUT"]["results"], BASE_NAME)

# Reference genome paths
REFERENCE_GENOME = config["INPUT"]["reference_genome"]
REF_BASE = os.path.splitext(REFERENCE_GENOME)[0]  # Remove extension
REF_INDEX_FILES = expand(f"{REF_BASE}.{{ext}}", ext=["amb", "ann", "bwt", "pac", "sa"])

# Genome length and windows files (for coverage analysis)
# Get the directory where the reference genome is located
REF_DIR = os.path.dirname(REFERENCE_GENOME)
GENOME_LENGTH = os.path.join(REF_DIR, f"{BASE_NAME}_length.txt")
WINDOWS_FILE = os.path.join(REF_DIR, f"{BASE_NAME}_{config['coverage_params']['window_file_suffix']}.txt")

# Sample discovery
SPECIES = os.path.basename(config["INPUT"]["reads"])
SAMPLES = sorted(
    d for d in os.listdir(READS_DIR)
    if (os.path.isdir(os.path.join(READS_DIR, d)) and 
        d.startswith(BASE_NAME))
)

ML_SAMPLES = [s for s in SAMPLES if "ML" in s]
FL_SAMPLES = [s for s in SAMPLES if "FL" in s]


def get_resource_path(resource_type, sample, read_num=None): 
    """Get path to either read files or reference genome."""
    if resource_type == "reads":
        if read_num not in ("1", "2"):
            raise ValueError("read_num must be '1' or '2'")
        files = glob.glob(os.path.join(READS_DIR, sample, f"*_{read_num}.fq.gz"))
        if not files:
            raise ValueError(f"No read files found for sample {sample}, read {read_num}")
        return files[0]
    
    if resource_type == "reference":
        return REFERENCE_GENOME
    
    raise ValueError("resource_type must be 'reads' or 'reference'")
