# =============================================================================
# 10 - KMER BBDUK: READ-LEVEL ASSEMBLY OF SEX-SPECIFIC K-MERS
#
# Uses the strictly sex-specific k-mers produced by step 06:
#   - male_specific_kmers.txt    (present in every male, absent in every female)
#   - female_specific_kmers.txt  (present in every female, absent in every male)
#
# These are KMC dumps (kmer<TAB>count), the k-mer in COLUMN 1. The old
# PLINK association chain that produced {sex}_specific_kmers.assoc (k-mer in
# column 2) was removed from step 06, so this module reads the .txt dumps.
# The counts in column 2 are set-operation minima across the present group
# and are not used here - only the sequences matter.
#
# Reads searched: the STANDARD trimmed reads from 01_trimming
# (01_trimmed/{sample}_{1,2}.fq.gz), NOT the k-mer branch reads from step 06
# (01_trimmed/{sample}_{1,2}_kmer.fq.gz). The k-mer branch hard-crops the
# first 15 bp of every mate to suppress BGI/MGI priming bias, which is right
# for counting k-mers but throws away sequence we want back when assembling
# full-length reads. A k-mer called on the cropped reads is still a genuine
# 31-mer of the uncropped read, so matching against the standard trimmed
# reads recovers the full ~150 bp read plus its mate.
#
# NOTE: the standard trimmed FASTQs are temp() in 01_trimming.smk. If they
# have already been consumed by the alignment branch and deleted, requesting
# this module re-runs fastp for the samples it needs.
#
# Workflow:
# 1. Convert sex-specific k-mer sequences to FASTA
# 2. BBDuk: pull reads matching sex-specific k-mers (per sample)
# 3. MEGAHIT: assemble matched reads per sample
# 4. CD-HIT: reduce redundant contigs across samples
# 5. BLAST: map reduced contigs to reference genome
#
# Compared to the ABySS assembly in step 06 (which assembles the raw 31 bp
# k-mers directly), this approach goes back to the full-length reads (~150 bp)
# with paired-end information. This yields longer, higher-quality contigs
# that are more informative for downstream BLAST and annotation.
# =============================================================================

# MEGAHIT is not one of the kmersGWAS helper binaries in binaries/, so it may
# live anywhere. Override with `megahit: /path/to/megahit` in the config.
MEGAHIT = config.get(
    "megahit",
    "/user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"
)


rule make_sexspecific_fasta:
    """
    Convert the sex-specific k-mer set from step 06 into FASTA.
    Input is the KMC dump {sex}_specific_kmers.txt: k-mer in column 1,
    minimum count across the present group in column 2.
    {sex} wildcard: male or female.
    """
    input:
        kmers = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "{sex}_specific_kmers.txt")
    output:
        fasta = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_kmers.fasta")
    wildcard_constraints:
        sex = "male|female"
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
        awk '{{print ">" NR "_{wildcards.sex}_kmer\\n" $1}}' {input.kmers} \
            > {output.fasta} 2> {log}
        echo "K-mers in FASTA: $(grep -c '^>' {output.fasta})" >> {log}
        """


rule bbduk_sexspecific:
    """
    Use BBDuk to extract reads containing sex-specific k-mers from the
    STANDARD trimmed FASTQs of 01_trimming (see the module header: not the
    _kmer.fq.gz reads, which are cropped by 15 bp at the 5' end).
    Runs per sample per sex. outm = matched reads for assembly.

    mm=f: BBDuk masks the middle base of every reference k-mer by default,
    which accepts a read carrying a different base at that position. The
    k-mers here were called by exact set membership, so the match must be
    exact too - a middle-base mismatch is a different k-mer, and one that
    step 06 may well have assigned to the other sex.

    rcomp=t is BBDuk's default and is stated explicitly: the KMC databases
    behind these k-mers are both-strand, so a read may carry either
    orientation.
    """
    input:
        r1 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_kmers.fasta")
    output:
        # temp(): the non-matching reads are a full copy of the input FASTQs
        # and nothing downstream reads them. BBDuk still writes them, they
        # are just not kept.
        discard1 = temp(os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_discard_1.fq.gz")),
        discard2 = temp(os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_discard_2.fq.gz")),
        matched1 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_1.fq"),
        matched2 = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_matched_2.fq"),
        stats = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_bbduk", "{sample}_stats")
    params:
        k = 31
    wildcard_constraints:
        sex = "male|female"
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
            mm=f \
            rcomp=t \
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
    wildcard_constraints:
        sex = "male|female"
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=2000,
        runtime=2040
    log:
        os.path.join(RESULTS_DIR, "logs", "10_kmer_bb", "{sex}_megahit_{sample}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {MEGAHIT} \
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
    Concatenate per-sample contigs and cluster with CD-HIT-EST at 95%
    identity to remove redundancy across samples. Male contigs are pooled
    over the males only, female contigs over the females only.
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
    wildcard_constraints:
        sex = "male|female"
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=1000,
        runtime=1000
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
            -M 16000 \
            -T {resources.cpus_per_task} >> {log} 2>&1
        echo "Input contigs: $(grep -c '^>' {output.all_contigs})" >> {log}
        echo "Reduced contigs: $(grep -c '^>' {output.reduced})" >> {log}
        """


rule blast_contigs:
    """
    BLAST reduced sex-specific contigs against reference genome.
    Reuses the BLAST database created in step 06 (rule blast_db).
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "10_kmer_bb", "{sex}_reduced_contigs.fasta"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "10_kmer_bb", "blast", "{sex}_contigs_blast.out")
    wildcard_constraints:
        sex = "male|female"
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
            -word_size 30 \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
            -evalue 1e-50 \
            -perc_identity 90 \
            -num_threads {resources.cpus_per_task} \
            -out {output.blast} 2> {log}
        """
