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
INTERVALS_DIR = config["intervals_dir"]
CHROM_PREFIX = config["chromosome_prefix"]
CONTIG_PREFIX = config["contig_prefix"]

# Discover samples: subdirectories in READS_DIR
SAMPLES = sorted([d for d in os.listdir(READS_DIR) if os.path.isdir(os.path.join(READS_DIR, d))])

# Discover intervals from pre-generated interval files
# Run scripts/generate_intervals.py manually before running the pipeline
INTERVALS = [os.path.basename(f).replace(".interval_list", "") 
             for f in glob.glob(os.path.join(INTERVALS_DIR, "*.interval_list"))]

# =============================================================================
# TARGET RULE
# =============================================================================
rule all:
    input:
        # Reference indices
        multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        REFERENCE + ".fai",
        REF_PREFIX + ".dict",
        # Trimmed reads
        expand(os.path.join(RESULTS_DIR, "trimmed", "{sample}_1.fq.gz"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "trimmed", "{sample}_2.fq.gz"), sample=SAMPLES),
        # Mapped and processed BAMs
        expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam.csi"), sample=SAMPLES),
        # Final joint-called variants
        os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz"),
        os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz.tbi"),
        os.path.join(RESULTS_DIR, "variants", "filtered.vcf.gz"),
        os.path.join(RESULTS_DIR, "variants", "filtered.vcf.gz.tbi")

# =============================================================================
# PREPROCESSING
# =============================================================================
rule trimmomatic:
    input:
        r1 = os.path.join(READS_DIR, "{sample}", "{sample}_R1.fastq.gz"),
        r2 = os.path.join(READS_DIR, "{sample}", "{sample}_R2.fastq.gz")
    output:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2.fq.gz"),
        r1_unpaired = os.path.join(RESULTS_DIR, "trimmed", "unpaired", "{sample}_1_unpaired.fq.gz"),
        r2_unpaired = os.path.join(RESULTS_DIR, "trimmed", "unpaired", "{sample}_2_unpaired.fq.gz")
    params:
        slidingwindow = config["trimmomatic"]["slidingwindow"],
        minlen = config["trimmomatic"]["minlen"],
        leading = config["trimmomatic"]["leading"],
        trailing = config["trimmomatic"]["trailing"]
    resources:
        cpus_per_task=16,
        mem_mb_per_cpu=4000
    log:
        os.path.join(RESULTS_DIR, "logs", "trimmomatic", "{sample}.log")
    envmodules:
        "Trimmomatic/0.39-Java-11"
    shell:
        """
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.r1_unpaired})
        mkdir -p $(dirname {log})
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads {resources.cpus_per_task} -phred33 \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            SLIDINGWINDOW:{params.slidingwindow} \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            MINLEN:{params.minlen} &> {log}
        """

# =============================================================================
# REFERENCE PREPARATION
# =============================================================================
rule bwa_mem2_index:
    input:
        REFERENCE
    output:
        multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        os.path.join(RESULTS_DIR, "logs", "bwa_mem2_index.log")
    envmodules:
        "bwa-mem2/2.2.1-intel-compilers-2023.1.0"
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=60000
    shell:
        """
        mkdir -p $(dirname {log})
        bwa-mem2 index -p {REF_PREFIX} {input} 2> {log}
        """

rule samtools_faidx:
    input:
        REFERENCE
    output:
        fai = REFERENCE + ".fai"
    log:
        os.path.join(RESULTS_DIR, "logs", "samtools_faidx.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        samtools faidx {input} 2> {log}
        """

rule create_sequence_dictionary:
    input:
        REFERENCE
    output:
        dict = REF_PREFIX + ".dict"
    log:
        os.path.join(RESULTS_DIR, "logs", "create_dict.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=10000
    shell:
        """
        mkdir -p $(dirname {log})
        gatk CreateSequenceDictionary \
            -R {input} \
            -O {output.dict} 2> {log}
        """

# =============================================================================
# MAPPING WITH GATK REQUIREMENTS
# =============================================================================
rule bwa_mem2_mem:
    input:
        reads = [
            os.path.join(RESULTS_DIR, "trimmed", "{sample}_1.fq.gz"),
            os.path.join(RESULTS_DIR, "trimmed", "{sample}_2.fq.gz")
        ],
        idx = multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        ref = REFERENCE
    output:
        bam = temp(os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam"))
    params:
        rg = lambda wildcards: fr"@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}"
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=3000,
        runtime=2000,
        partition="zen4"
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_bwa_mem2.log")
    envmodules:
        "bwa-mem2/2.2.1-intel-compilers-2023.1.0",
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        bwa-mem2 mem -t {resources.cpus_per_task} -R '{params.rg}' \
            {REF_PREFIX} {input.reads[0]} {input.reads[1]} 2> {log} | \
        samtools sort -@ {resources.cpus_per_task} -o {output.bam} -
        """

rule mark_duplicates:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam")
    output:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam"),
        metrics = os.path.join(RESULTS_DIR, "mapped", "metrics", "{sample}.dedup_metrics.txt")
    resources:
        cpus_per_task=8,
        mem_mb_per_cpu=2000
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_markdup.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX false 2> {log}
        """

rule index_bam:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam")
    output:
        bai = os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam.csi")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_index.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        samtools index -c -@ {resources.cpus_per_task} {input.bam} 2> {log}
        """

# =============================================================================
# VARIANT CALLING WITH GATK - PARALLELIZED BY INTERVAL
# 
# Workflow:
# 1. HaplotypeCaller per sample per interval -> per-sample per-interval GVCFs
# 2. GenomicsDBImport per interval (combines all samples) -> per-interval GenomicsDB
# 3. GenotypeGVCFs per interval -> per-interval VCFs
# 4. MergeVcfs -> single raw.vcf.gz
# 5. VariantFiltration -> filtered.vcf.gz
# =============================================================================

# Step 1: HaplotypeCaller per sample per interval
rule haplotypecaller_per_interval:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam"),
        bai = os.path.join(RESULTS_DIR, "mapped", "{sample}.dedup.bam.csi"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict",
        interval = os.path.join(INTERVALS_DIR, "{interval}.interval_list")
    output:
        gvcf = os.path.join(RESULTS_DIR, "variants", "gvcf", "per_interval", "{sample}", "{interval}.g.vcf.gz")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=4000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "haplotypecaller", "{sample}", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus_per_task}
        mkdir -p $(dirname {output.gvcf})
        mkdir -p $(dirname {log})
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -L {input.interval} \
            -ERC GVCF \
            --create-output-variant-index false \
            --native-pair-hmm-threads {resources.cpus_per_task} 2> {log}
        """

# Step 1.1: Index GVCFs with CSI format
rule index_gvcf_csi:
    input:
        gvcf = os.path.join(RESULTS_DIR, "variants", "gvcf", "per_interval", "{sample}", "{interval}.g.vcf.gz")
    output:
        csi = os.path.join(RESULTS_DIR, "variants", "gvcf", "per_interval", "{sample}", "{interval}.g.vcf.gz.csi")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "index_gvcf", "{sample}", "{interval}.log")
    envmodules:
        "BCFtools/1.21-GCC-13.3.0"  # or whatever version you have
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools index --csi {input.gvcf} 2> {log}
        echo "Index created: {output.csi}" >> {log}
        """
        
# Step 2: GenomicsDBImport per interval (combines all samples for that interval)
rule genomicsdb_import_per_interval:
    input:
        gvcfs = expand(
            os.path.join(RESULTS_DIR, "variants", "gvcf", "per_interval", "{sample}", "{{interval}}.g.vcf.gz"),
            sample=SAMPLES
        ),
        csis = expand(
        os.path.join(RESULTS_DIR, "variants", "gvcf", "per_interval", "{sample}", "{{interval}}.g.vcf.gz.csi"),
        sample=SAMPLES
        ),
        interval = os.path.join(INTERVALS_DIR, "{interval}.interval_list"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        db = directory(os.path.join(RESULTS_DIR, "variants", "genomicsdb", "{interval}"))
    params:
        gvcf_args = lambda wildcards, input: " ".join([f"-V {g}" for g in input.gvcfs])
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "genomicsdb", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.db})
        mkdir -p $(dirname {output.db})/tmp
        mkdir -p $(dirname {log})
        gatk GenomicsDBImport \
            {params.gvcf_args} \
            --genomicsdb-workspace-path {output.db} \
            -L {input.interval} \
            --reader-threads {resources.cpus_per_task} \
            --batch-size 50 \
            --tmp-dir $(dirname {output.db})/tmp 2> {log}
        """

# Step 3: GenotypeGVCFs per interval
rule genotype_gvcfs_per_interval:
    input:
        db = os.path.join(RESULTS_DIR, "variants", "genomicsdb", "{interval}"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "variants", "per_interval", "{interval}.vcf.gz"),
        tbi = os.path.join(RESULTS_DIR, "variants", "per_interval", "{interval}.vcf.gz.tbi")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "genotype", "{interval}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output.vcf} 2> {log}
        """

# Step 4: Merge all per-interval VCFs into one raw VCF
rule merge_vcfs:
    input:
        vcfs = expand(
            os.path.join(RESULTS_DIR, "variants", "per_interval", "{interval}.vcf.gz"),
            interval=INTERVALS
        ),
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz")
    params:
        vcf_args = lambda wildcards, input: " ".join([f"-I {v}" for v in input.vcfs])
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "merge_vcfs.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        mkdir -p $(dirname {log})
        gatk MergeVcfs \
            {params.vcf_args} \
            -D {input.dict} \
            -O {output.vcf} 2> {log}
        """

rule index_raw_vcf:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz")
    output:
        tbi = os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz.tbi")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "index_raw.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf} 2> {log}
        """

# Step 5: Hard filtering
rule variant_filtration:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz"),
        tbi = os.path.join(RESULTS_DIR, "variants", "raw.vcf.gz.tbi"),
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        vcf = os.path.join(RESULTS_DIR, "variants", "filtered.vcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "filtration.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-name "QD_filter" --filter-expression "QD < 2.0" \
            --filter-name "FS_filter" --filter-expression "FS > 60.0" \
            --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
            --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" 2> {log}
        """

rule index_filtered_vcf:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "filtered.vcf.gz")
    output:
        tbi = os.path.join(RESULTS_DIR, "variants", "filtered.vcf.gz.tbi")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "index_filtered.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf} 2> {log}
        """
