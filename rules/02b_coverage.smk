# =============================================================================
# 02b - COVERAGE STATISTICS
# =============================================================================

rule genome_coverage:
    """
    Compute coverage histogram per chromosome using bedtools genomecov.
    Only includes chromosomes matching the configured prefix.
    """
    input:
        bam = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam"),
        bai = os.path.join(RESULTS_DIR, "02_aligned", "{sample}.dedup.bam.csi"),
        fai = REFERENCE + ".fai"
    output:
        cov = os.path.join(RESULTS_DIR, "02_coverage", "{sample}.coverage.txt")
    params:
        chrom_prefix = CHROM_PREFIX
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=8000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "02_coverage", "{sample}.coverage.log")
    envmodules:
        "BEDTools/2.31.0-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.cov})
        mkdir -p $(dirname {log})
        
        # Create genome file with only chromosomes matching prefix
        awk -v prefix="{params.chrom_prefix}" '$1 ~ "^"prefix {{print $1"\\t"$2}}' \
            {input.fai} > {output.cov}.genome.txt 2> {log}
        
        genomeCoverageBed -ibam {input.bam} -g {output.cov}.genome.txt \
            > {output.cov} 2>> {log}
        
        rm -f {output.cov}.genome.txt
        """


rule coverage_summary:
    """
    Summarize bedtools genomecov histogram into per-chromosome and genome-wide
    mean coverage. Only uses chromosome lines (matching prefix) and the 'genome' line.
    """
    input:
        cov = os.path.join(RESULTS_DIR, "02_coverage", "{sample}.coverage.txt")
    output:
        stats = os.path.join(RESULTS_DIR, "02_coverage", "{sample}.coverage_summary.tsv")
    params:
        chrom_prefix = CHROM_PREFIX
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=1000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "02_coverage", "{sample}.coverage_summary.log")
    shell:
        """
        mkdir -p $(dirname {log})
        awk -v prefix="{params.chrom_prefix}" -v sample="{wildcards.sample}" '
        BEGIN {{
            OFS="\\t"
            print "sample\\tchromosome\\tmean_depth\\tbreadth_pct"
        }}
        ($1 ~ "^"prefix || $1 == "genome") {{
            # col1=chrom, col2=depth, col3=bases_at_depth, col4=chrom_size, col5=fraction
            weighted[($1)] += ($2) * ($3)
            size[($1)] = ($4)
            if ($2 > 0) covered[($1)] += ($3)
        }}
        END {{
            for (chr in size) {{
                mean = weighted[chr] / size[chr]
                breadth = (covered[chr]+0) / size[chr] * 100
                print sample "\\t" chr "\\t" sprintf("%.2f", mean) "\\t" sprintf("%.2f", breadth)
            }}
        }}' {input.cov} > {output.stats} 2> {log}
        """


rule aggregate_coverage:
    """
    Combine all per-sample coverage summaries into one table.
    """
    input:
        stats = expand(os.path.join(RESULTS_DIR, "02_coverage", "{sample}.coverage_summary.tsv"), sample=SAMPLES)
    output:
        summary = os.path.join(RESULTS_DIR, "02_coverage", "all_samples_coverage_summary.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "02_coverage", "aggregate_coverage.log")
    shell:
        """
        mkdir -p $(dirname {log})
        head -1 {input.stats[0]} > {output.summary} 2> {log}
        for f in {input.stats}; do
            tail -n +2 "$f" >> {output.summary}
        done 2>> {log}
        echo "Aggregated $(echo {input.stats} | wc -w) samples" >> {log}
        """