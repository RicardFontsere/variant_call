# =============================================================================
# 09 - KMER BBDUK: READ-LEVEL ASSEMBLY OF SEX-SPECIFIC K-MERS
#
# Uses the strictly sex-specific k-mers already identified in step 06:
#   - male_specific_kmers.assoc   (F_A=1, F_U=0: present in all males, absent in all females)
#   - female_specific_kmers.assoc (F_A=0, F_U=1: absent in all males, present in all females)
#
# These are the same criteria a BBDuk+R validation loop would arrive at,
# since KMC's ci2 threshold already removes single-read noise.
#
# Workflow:
# 1. Convert sex-specific k-mer sequences to FASTA
# 2. BBDuk: pull reads matching sex-specific k-mers (per sample)
# 3. MEGAHIT: assemble matched reads per sample
# 4. CD-HIT: reduce redundant contigs across samples (Note changed to 95% identity and added -aS 0.8 minimum 80% needs to match. Different from Sophie)
# 5. BLAST: map reduced contigs to reference genome
#
# Compared to the ABySS assembly in step 06 (which assembles the raw 31bp
# k-mers directly), this approach goes back to the full-length reads (~150bp)
# with paired-end information. This yields longer, higher-quality contigs
# that are more informative for downstream BLAST and annotation.
# =============================================================================


rule make_sexspecific_fasta:
    """
    Convert sex-specific k-mer association results from step 06 into FASTA.
    Column 2 (SNP) in the .assoc file contains the k-mer sequence.
    {sex} wildcard: male or female.
    """
    input:
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "{sex}_specific_kmers.assoc")
    output:
        fasta = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_kmers.fasta")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "make_{sex}_fasta.log")
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})
        awk '{{print ">" NR "_{wildcards.sex}_kmer\\n" $2}}' {input.assoc} \
            > {output.fasta} 2> {log}
        echo "K-mers in FASTA: $(grep -c '^>' {output.fasta})" >> {log}
        """


rule bbduk_sexspecific:
    """
    Use BBDuk to extract reads containing sex-specific k-mers from trimmed FASTQs.
    Runs per sample per sex. outm = matched reads for assembly.
    """
    input:
        r1 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_kmers.fasta")
    output:
        discard1 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_discard_1.fq.gz"),
        discard2 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_discard_2.fq.gz"),
        matched1 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_1.fq"),
        matched2 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_2.fq"),
        stats = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_stats")
    params:
        k = 31
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=40000,
        runtime=720
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "{sex}_bbduk_{sample}.log")
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.matched1})
        mkdir -p $(dirname {log})
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.discard1} \
            out2={output.discard2} \
            outm={output.matched1} \
            outm2={output.matched2} \
            ref={input.ref} \
            k={params.k} \
            stats={output.stats} \
            usejni=t &> {log}
        """


rule megahit_sexspecific:
    """
    Assemble sex-specific matched reads per sample using MEGAHIT.
    Separated from BBDuk so failed assemblies don't re-run matching.
    """
    input:
        r1 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_1.fq"),
        r2 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_2.fq")
    output:
        contigs = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_megahit", "{sample}", "{sample}.contigs.fa")
    params:
        outdir = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_megahit", "{sample}")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=2000,
        runtime=2040
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "{sex}_megahit_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        /user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -f \
            -t {resources.cpus_per_task} \
            -m 1 \
            --out-dir {params.outdir} \
            --out-prefix {wildcards.sample} 2> {log}
        """


rule reduce_contigs:
    """
    Concatenate per-sample contigs and cluster with CD-HIT-EST at 98% identity
    to remove redundancy across samples.
    """
    input:
        contigs = lambda wildcards: expand(
            os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_megahit", "{sample}", "{sample}.contigs.fa"),
            sex=wildcards.sex,
            sample=ML_SAMPLES if wildcards.sex == "male" else FL_SAMPLES
        )
    output:
        all_contigs = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_all_contigs.fasta"),
        reduced = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_reduced_contigs.fasta")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=1000,
        runtime=1200
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "reduce_{sex}_contigs.log")
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        cat {input.contigs} > {output.all_contigs} 2> {log}
        cd-hit-est \
            -i {output.all_contigs} \
            -o {output.reduced} \
            -c 0.95 -n 10 -d 0 \
            -aS 0.8 \
            -M 16000 \
            -T {resources.cpus_per_task} >> {log} 2>&1
        echo "Input contigs: $(grep -c '^>' {output.all_contigs})" >> {log}
        echo "Reduced contigs: $(grep -c '^>' {output.reduced})" >> {log}
        """


rule blast_contigs:
    """
    BLAST reduced sex-specific contigs against reference genome.
    Reuses the BLAST database created in step 06.
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_reduced_contigs.fasta"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "10_kmer_bb", "blast", "{sex}_contigs_blast.out")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=36000,
        runtime=4000
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "blast_{sex}_contigs.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        blastn \
            -query {input.contigs} \
            -db {input.ref} \
            -task megablast \
            -word_size 16 \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
            -evalue 1e-5 \
            -perc_identity 90 \
            -num_threads 20 \
            -out {output.blast} 2> {log}
        """
