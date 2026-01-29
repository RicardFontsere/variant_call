# =============================================================================
# 00 - REFERENCE PREPARATION
# =============================================================================

rule bwa_index:
    input:
        REFERENCE
    output:
        multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        os.path.join(RESULTS_DIR, "logs", "00_bwa_index.log")
    envmodules:
        "bwa-mem2/2.2.1-intel-compilers-2023.1.0"
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=40000,
    shell:
        """
        mkdir -p $(dirname {log})
        bwa-mem2.avx2 index -p {REF_PREFIX} {input} 2> {log}
        """

rule samtools_faidx:
    input:
        REFERENCE
    output:
        fai = REFERENCE + ".fai"
    log:
        os.path.join(RESULTS_DIR, "logs", "00_samtools_faidx.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        samtools faidx {input} 2> {log}
        """

rule sequence_dictionary:
    input:
        REFERENCE
    output:
        dict = REF_PREFIX + ".dict"
    log:
        os.path.join(RESULTS_DIR, "logs", "00_sequence_dictionary.log")
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

checkpoint generate_intervals:
    input:
        fai = REFERENCE + ".fai",
        dict = REF_PREFIX + ".dict"
    output:
        directory(INTERVALS_DIR)
    params:
        chrom_prefix = CHROM_PREFIX,
        contig_prefix = CONTIG_PREFIX
    log:
        os.path.join(RESULTS_DIR, "logs", "00_generate_intervals.log")
    shell:
        """
        mkdir -p $(dirname {log})
        python scripts/generate_intervals.py \
            {input.fai} {input.dict} {output} \
            {params.chrom_prefix} {params.contig_prefix} 2> {log}
        """