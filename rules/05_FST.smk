# =============================================================================
# 05 - SEX-BASED VARIANT SPLITTING
# 
# Workflow:
# 1. Extract sample names and create male/female sample lists
# 2. Split VCF by sex and filter for biallelic SNPs only
# =============================================================================

rule create_sex_sample_lists:
    """
    Extract sample names from VCF and create separate lists for males and females
    based on naming patterns.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    output:
        samples = os.path.join(RESULTS_DIR, "03_variants", "all_samples.txt"),
        males = os.path.join(RESULTS_DIR, "03_variants", "males.txt"),
        females = os.path.join(RESULTS_DIR, "03_variants", "females.txt")
    params:
        male_pattern = config.get("male_pattern", "ML"),
        female_pattern = config.get("female_pattern", "FL")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "05_sex_split", "create_sample_lists.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        
        # Extract all sample names
        bcftools query -l {input.vcf} > {output.samples} 2> {log}
        
        # Create male and female sample lists
        grep "{params.male_pattern}" {output.samples} > {output.males} 2>> {log} || true
        grep "{params.female_pattern}" {output.samples} > {output.females} 2>> {log} || true
        
        # Log counts
        echo "Total samples: $(wc -l < {output.samples})" >> {log}
        echo "Males: $(wc -l < {output.males})" >> {log}
        echo "Females: $(wc -l < {output.females})" >> {log}
        """

rule split_vcf_males:
    """
    Extract male samples and filter for SNPs only.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        samples = os.path.join(RESULTS_DIR, "03_variants", "males.txt")
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "males_biallelic.bcf")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "05_sex_split", "split_males.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view \
            -S {input.samples} \
            -v snps \
            -o {output.vcf} \
            {input.vcf} 2> {log}
        """

rule split_vcf_females:
    """
    Extract female samples and filter for SNPs only.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        samples = os.path.join(RESULTS_DIR, "03_variants", "females.txt")
    output:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "females_biallelic.bcf")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "05_sex_split", "split_females.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view \
            -S {input.samples} \
            -v snps \
            -o {output.vcf} \
            {input.vcf} 2> {log}
        """

rule fst_analysis:
    """
    Run FST analysis comparing males vs females using FSTest.
    Uses Nei's estimator with sliding window approach.
    
    Parameters:
        --m 2: Nei's FST estimator
        --zt 1: Z-transformation
    """
    input:
        males = os.path.join(RESULTS_DIR, "03_variants", "males_biallelic.bcf"),
        females = os.path.join(RESULTS_DIR, "03_variants", "females_biallelic.bcf")
    output:
        fst = os.path.join(RESULTS_DIR, "05_FST", "males_vs_females.snp")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "05_fst", "males_vs_females"),
        window = 50,
        step = 25
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=400000,
        runtime=6000
    log:
        os.path.join(RESULTS_DIR, "logs", "05_fst", "fst_analysis.log")
    shell:
        """
        mkdir -p $(dirname {output.fst})
        mkdir -p $(dirname {log})
        /user/brussel/109/vsc10945/home/scratch/Software/FSTest/FSTest1.3 \
            --pop1 {input.males} \
            --pop2 {input.females} \
            --m 2 \
            --zt 1 \
            --o {params.out_prefix} 2> {log}
        """
