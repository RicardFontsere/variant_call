# =============================================================================
# 02 - ALIGNMENT
# =============================================================================

rule align:
    input:
        reads = [
            os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1.fq.gz"),
            os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2.fq.gz")
        ],
        idx = multiext(REF_PREFIX, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        ref = REFERENCE
    output:
        bam = temp(os.path.join(RESULTS_DIR, "02_aligned", "{sample}.sorted.bam"))
    params:
        rg = lambda wildcards: fr"@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}"
    resources:
        cpus_per_task=10,
        mem_mb_per_cpu=10000,
        runtime=1000,
        partition="zen4"
    log:
        os.path.join(RESULTS_DIR, "logs", "02_align", "{sample}.log")
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
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.sorted.bam")
    output:
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam"),
        metrics = os.path.join(RESULTS_DIR, "00_qc", "{sample}.dedup_metrics.txt")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=10000
    log:
        os.path.join(RESULTS_DIR, "logs", "02_mark_duplicates", "{sample}.log")
    envmodules:
        "GATK/4.5.0.0-GCCcore-12.3.0-Java-17"
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        mkdir -p $(dirname {log})
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX false 2> {log}
        """

rule index_bam:
    input:
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam")
    output:
        bai = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam.csi")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=10000
    log:
        os.path.join(RESULTS_DIR, "logs", "02_index_bam", "{sample}.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        samtools index -c -@ {resources.cpus_per_task} {input.bam} 2> {log}
        """
