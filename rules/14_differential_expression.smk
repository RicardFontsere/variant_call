# =============================================================================
# 14 - DIFFERENTIAL EXPRESSION (edgeR)
#
# One analysis per entry of `de_analyses` in the config. For an analysis named
# <A>, <de_input_dir> must hold:
#   <A>_design.txt  sample <tab> ... <tab> group   (one line per sample, the
#                   sample names are the RNA-seq sample names)
#   <A>_matrix.txt  model_coefficients <tab> one column per contrast, with the
#                   group weights (1 / -1 / 0) below
#
# Run with:  snakemake de_all --profile master/profile
# =============================================================================

DE_DIR = os.path.join(RESULTS_DIR, "14_DE")
DE_STAGED = os.path.join(DE_DIR, "input")

rule de_all:
    input:
        os.path.join(RESULTS_DIR, "13_rnaseq", "qc", "multiqc_report.html"),
        expand(os.path.join(DE_DIR, "{analysis}", "{analysis}_DE_results.txt"), analysis=DE_ANALYSES)

rule de_input:
    """
    Reduce the Salmon matrices to the transcripts of the annotation and to the
    samples of the design. The design and the contrast matrix are staged next to
    the count file, so that everything edgeR reads sits in one directory.
    """
    input:
        salmon_counts = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonCountNumReads.txt"),
        salmon_tpm = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonTPM.txt"),
        salmon_efflen = os.path.join(RESULTS_DIR, "13_rnaseq", "salmon", "salmonEffLength.txt"),
        annotation = RNA_ANNOTATION,
        design = os.path.join(DE_INPUT_DIR, "{analysis}_design.txt"),
        matrix = os.path.join(DE_INPUT_DIR, "{analysis}_matrix.txt")
    output:
        counts = os.path.join(DE_STAGED, "{analysis}_count.txt"),
        annotation = os.path.join(DE_STAGED, "{analysis}_annotation.txt"),
        tpm = os.path.join(DE_STAGED, "{analysis}_TPM_annotated.txt"),
        design = os.path.join(DE_STAGED, "{analysis}_design.txt"),
        matrix = os.path.join(DE_STAGED, "{analysis}_matrix.txt")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "14_DE", "de_input_{analysis}.log")
    envmodules:
        "R-bundle-CRAN/2025.10-foss-2025a"
    shell:
        """
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        Rscript --vanilla {SCRIPTS_DIR}/prepare_de_input.R \
            {input.salmon_counts} \
            {input.salmon_tpm} \
            {input.salmon_efflen} \
            {input.annotation} \
            {input.design} \
            {output.counts} \
            {output.annotation} \
            {output.tpm} &> {log}
        cp {input.design} {output.design}
        cp {input.matrix} {output.matrix}
        """

rule edger:
    """
    edgeR quasi-likelihood test per contrast. Writes the DE tables, the GO input
    files, the diagnostic plots and (if de_heatmaps is 'yes') the heatmaps.
    """
    input:
        counts = os.path.join(DE_STAGED, "{analysis}_count.txt"),
        annotation = os.path.join(DE_STAGED, "{analysis}_annotation.txt"),
        design = os.path.join(DE_STAGED, "{analysis}_design.txt"),
        matrix = os.path.join(DE_STAGED, "{analysis}_matrix.txt")
    output:
        results = os.path.join(DE_DIR, "{analysis}", "{analysis}_DE_results.txt")
    params:
        datapath = DE_STAGED,
        outpath = os.path.join(DE_DIR, "{analysis}"),
        fdr = DE_FDR,
        filtering = DE_FILTERING,
        heatmaps = DE_HEATMAPS
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=32000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "14_DE", "edger_{analysis}.log")
    envmodules:
        BIOCONDUCTOR_MODULE
    shell:
        """
        mkdir -p {params.outpath} $(dirname {log})
        Rscript --vanilla {SCRIPTS_DIR}/edger_de.R \
            {wildcards.analysis} \
            {params.fdr} \
            {params.filtering} \
            {params.heatmaps} \
            {params.datapath} \
            {params.outpath} \
            {input.annotation} &> {log}
        """
