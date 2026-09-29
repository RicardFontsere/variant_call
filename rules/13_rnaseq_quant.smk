# =============================================================================
# 13 - RNA-seq QUANTIFICATION (fastp + Salmon)
#
# The reads are already on disk, one pair per sample:
#   {rna_reads_dir}/{sample}_1.fastq.gz and {rna_reads_dir}/{sample}_2.fastq.gz
#
# fastp replaces the Trimmomatic + FastQC steps of the original pipeline: it
# trims and reports the before/after metrics in one pass, and MultiQC reads its
# json directly.
# =============================================================================

rule rna_fastp:
    """
    Trim the RNA-seq reads. The flags reproduce the Trimmomatic settings that
    were used before: HEADCROP:12 -> --trim_front, SLIDINGWINDOW:4:15 ->
    --cut_right, MINLEN:36 -> --length_required. Adapters are detected from the
    reads instead of being read from adapters.fasta.
    """
    input:
        # ancient(): keeps a re-sync of the raw reads from re-running everything,
        # see the comment on rule fastp in 01_trimming.smk
        r1 = ancient(os.path.join(RNA_READS_DIR, "{sample}_1.fastq.gz")),
        r2 = ancient(os.path.join(RNA_READS_DIR, "{sample}_2.fastq.gz"))
    output:
        # temp(): deleted once salmon_quant has consumed them, the json/html QC
        # reports are the permanent record
        r1 = temp(os.path.join(RESULTS_DIR, "13_rnaseq", "trimmed", "{sample}_1.fq.gz")),
        r2 = temp(os.path.join(RESULTS_DIR, "13_rnaseq", "trimmed", "{sample}_2.fq.gz")),
        html = os.path.join(RESULTS_DIR, "13_rnaseq", "qc", "{sample}.html"),
        json = os.path.join(RESULTS_DIR, "13_rnaseq", "qc", "{sample}.json")
    params:
        trim_front = RNA_TRIM_FRONT
    resources:
        cpus_per_task=10,
        mem_mb_per_cpu=3200,
        runtime=720
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "fastp_{sample}.log")
    envmodules:
        "fastp/1.0.1-GCC-13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {output.html}) $(dirname {log})
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --detect_adapter_for_pe \
            --trim_front1 {params.trim_front} \
            --trim_front2 {params.trim_front} \
            --cut_right \
            --cut_right_window_size 4 \
            --cut_right_mean_quality 15 \
            --length_required 36 \
            --thread {resources.cpus_per_task} \
            --html {output.html} \
            --json {output.json} \
            &> {log}
        """

rule rna_multiqc:
    """
    One MultiQC report over the fastp json files of all samples.
    """
    input:
        expand(os.path.join(RESULTS_DIR, "13_rnaseq", "qc", "{sample}.json"), sample=RNA_SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "13_rnaseq", "qc", "multiqc_report.html")
    params:
        qc_dir = os.path.join(RESULTS_DIR, "13_rnaseq", "qc")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "multiqc.log")
    envmodules:
        "MultiQC/1.28-foss-2024a"
    shell:
        """
        mkdir -p $(dirname {log})
        multiqc {params.qc_dir} -o {params.qc_dir} -n $(basename {output}) -f &> {log}
        """

rule rna_transcriptome:
    """
    Extract the mRNA sequences from the genome with gffread (-E cleans up the
    features, -w writes the spliced exons).
    """
    input:
        genome = RNA_GENOME,
        gff = RNA_GFF
    output:
        fasta = os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "transcriptome.fa")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=10000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "transcriptome.log")
    envmodules:
        "SAMtools/1.18-GCC-12.3.0",
        "gffread"
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        # gffread is much faster with a genome index next to the fasta. Only
        # built if absent, writing it would restamp the index the DNA rules use.
        [ -f {input.genome}.fai ] || samtools faidx {input.genome}
        gffread -E -w {output.fasta} -g {input.genome} {input.gff} &> {log}
        """

rule salmon_index:
    """
    Decoy-aware index: the genome is appended to the transcriptome so that reads
    coming from unannotated regions are not forced onto a transcript.
    """
    input:
        transcriptome = os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "transcriptome.fa"),
        genome = RNA_GENOME
    output:
        index_dir = directory(os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "salmon_index"))
    params:
        decoys = os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "decoys.txt"),
        gentrome = os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "gentrome.fa")
    resources:
        cpus_per_task=8,
        mem_mb_per_cpu=10000,
        runtime=720
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "salmon_index.log")
    envmodules:
        "Salmon/1.10.3-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        grep '^>' {input.genome} | cut -d ' ' -f1 | sed 's/^>//' > {params.decoys}
        cat {input.transcriptome} {input.genome} > {params.gentrome}
        salmon index \
            -t {params.gentrome} \
            -d {params.decoys} \
            -p {resources.cpus_per_task} \
            -i {output.index_dir} &> {log}
        # the gentrome is a copy of the genome, the index does not need it
        rm -f {params.gentrome}
        """

rule salmon_quant:
    """
    Quantify one sample. -l A lets salmon infer the library type.
    """
    input:
        index_dir = os.path.join(RESULTS_DIR, "13_rnaseq", "reference", "salmon_index"),
        r1 = os.path.join(RESULTS_DIR, "13_rnaseq", "trimmed", "{sample}_1.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "13_rnaseq", "trimmed", "{sample}_2.fq.gz")
    output:
        sf = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "{sample}", "quant.sf")
    params:
        outdir = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "{sample}")
    resources:
        cpus_per_task=8,
        mem_mb_per_cpu=4000,
        runtime=720
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "salmon_quant_{sample}.log")
    envmodules:
        "Salmon/1.10.3-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        salmon quant \
            -i {input.index_dir} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {resources.cpus_per_task} \
            -o {params.outdir} \
            --gcBias &> {log}
        """

rule salmon_compile:
    """
    One matrix per quantity, transcripts in rows and samples in columns.
    quant.sf columns: 1 Name, 2 Length, 3 EffectiveLength, 4 TPM, 5 NumReads.
    All samples share the index, so the transcript order is the same in every
    file and the columns can simply be pasted together.
    """
    input:
        sf = expand(os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "{sample}", "quant.sf"), sample=RNA_SAMPLES)
    output:
        counts = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonCountNumReads.txt"),
        tpm = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonTPM.txt"),
        efflen = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonEffLength.txt")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=10000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "13_rnaseq", "salmon_compile.log")
    shell:
        """
        mkdir -p $(dirname {log})
        tmp=$(mktemp -d)
        {{ echo id; tail -n +2 {input.sf[0]} | cut -f1; }} > $tmp/id

        for sf in {input.sf}; do
            sample=$(basename $(dirname $sf))
            for col in 3 4 5; do
                {{ echo "$sample"; tail -n +2 "$sf" | cut -f$col; }} > "$tmp/$sample.$col"
            done
        done

        paste $tmp/id $tmp/*.3 > {output.efflen}
        paste $tmp/id $tmp/*.4 > {output.tpm}
        paste $tmp/id $tmp/*.5 > {output.counts}
        rm -rf $tmp
        """
