# =============================================================================
# 01 - TRIMMING
# =============================================================================

rule fastp:
    input:
        # ancient(): the raw FASTQs are >100 kB, so Snakemake never checksums
        # them (io.py is_checksum_eligible) and falls back to mtime alone. Any
        # re-sync or re-symlink of the source data then bumps their mtime and
        # re-triggers fastp -> align -> the entire downstream pipeline.
        # ancient() makes is_newer() return False, pinning that trigger off.
        r1 = lambda wildcards: ancient(get_read_file(wildcards.sample, "1")),
        r2 = lambda wildcards: ancient(get_read_file(wildcards.sample, "2"))
    output:
        # temp(): deleted once rule align has consumed them. The QC html/json
        # below are small and kept as the permanent QC record.
        r1 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1.fq.gz")),
        r2 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2.fq.gz")),
        html = os.path.join(RESULTS_DIR, "00_qc", "{sample}.html"),
        json = os.path.join(RESULTS_DIR, "00_qc", "{sample}.json")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=2000,
        runtime=720
    log:
        os.path.join(RESULTS_DIR, "logs", "01_fastp", "{sample}.log")
    envmodules:
        "fastp/1.0.1-GCC-13.3.0" 
    shell:
        """
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.html})
        mkdir -p $(dirname {log})
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --thread {resources.cpus_per_task} \
            --html {output.html} \
            --json {output.json} \
            &> {log}
        """
