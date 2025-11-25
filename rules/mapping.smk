rule bwa_index:
    input:
        REFERENCE_GENOME
    output:
        REF_INDEX_FILES
    params:
        index_prefix = REF_BASE
    log:
        os.path.join(RESULTS_DIR, "logs", "bwa_index.log")
    envmodules:
        "BWA/0.7.17-GCCcore-11.3.0"
    shell:
        """
        mkdir -p {RESULTS_DIR}/logs/bwa
        bwa index -p {params.index_prefix} {input} 2> {log}
        """

rule bwa_mem:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        ref = REFERENCE_GENOME,
        idx = REF_INDEX_FILES
    output:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.unsorted.bam")
    resources:
        cpus_per_task = 64,
        mem_mb_per_cpu = 4000,  
        partition = "zen5_mpi"
    params:
        min_mapq = 20,
        ref_prefix = REF_BASE
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_mem.log")
    envmodules:
        "BWA/0.7.17-GCCcore-12.3.0",
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p {RESULTS_DIR}/mapped
        mkdir -p {RESULTS_DIR}/logs/bwa
        bwa mem -t {resources.cpus_per_task} {params.ref_prefix} {input.r1} {input.r2} | \
        /user/brussel/109/vsc10945/home/scratch/Software/samblaster/samblaster -r | \
        samtools view -@ {resources.cpus_per_task} -bh -q {params.min_mapq} -o {output.bam} # add "-F 0x100 -F 0x800" for removing secondary and supplementary alignments
        """

rule sort_bam:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.unsorted.bam")
    output:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam")
    resources:
        cpus_per_task = 20,
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_sort.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        samtools sort -@ {resources.cpus_per_task} -o {output.bam} {input.bam} 
        rm {input.bam} # remove unsorted bam to save space, have to check if this works in snakemake
        """

rule index_bam:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam")
    output:
        bai = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam.csi")
    resources:
        cpus_per_task = 20,
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}_index.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0"
    shell:
        """
        samtools index -@ {resources.cpus_per_task} -c {input.bam}
        """

