# =============================================================================
# 06 - SEX-SPECIFIC K-MERS (KMC set operations)
#
# Workflow:
# 1. fastp preprocessing of the k-mer branch reads
# 2. KMC3 k-mer counting (canonical + all) per sample
# 3. kmc_tools set operations -> strictly sex-specific k-mers
# 4. ABySS assembly of sex-specific k-mers
# 5. BLAST assembled contigs against reference
#
# The kmersGWAS/PLINK association chain (add_strand_information ->
# combine_kmers -> build_kmers_table -> kmers_table_to_bed -> PLINK --assoc)
# was removed: with n=8 there is no association power, and "present in all
# males, absent in all females" is a deterministic set-membership question
# that kmc_tools answers directly and far more cheaply.
#
# NOTE for the currently disabled modules 10_kmer_bb.smk and
# 11_kmer_coverage.smk: they read {sex}_specific_kmers.assoc and take the
# k-mer from column 2. The new outputs are {sex}_specific_kmers.txt with the
# k-mer in column 1 (KMC dump format: kmer<TAB>count). Either adjust their
# awk to `$1`, or point them at the ready-made FASTA {sex}_abyss.input.
# =============================================================================

rule fastp_kmer:
    """
    Preprocess BGI/MGI reads for the k-mer branch.
    BGI/MGI libraries show non-random base composition in the first
    ~10-15 bp (random-hexamer priming bias), which produces spurious
    k-mers. fastp performs adapter + quality trimming and hard-crops the
    first 15 bp of both mates.
    """
    input:
        # ancient(): see rule fastp in 01_trimming.smk - raw FASTQs are too
        # large to checksum, so a bumped mtime alone would re-trigger the
        # whole k-mer branch.
        r1 = lambda wildcards: ancient(get_read_file(wildcards.sample, "1")),
        r2 = lambda wildcards: ancient(get_read_file(wildcards.sample, "2"))
    output:
        # temp(): deleted once rule kmer_analysis has consumed them. The QC
        # html/json below are small and kept as the permanent QC record.
        r1 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1_kmer.fq.gz")),
        r2 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2_kmer.fq.gz")),
        html = os.path.join(RESULTS_DIR, "00_qc", "{sample}_kmer.html"),
        json = os.path.join(RESULTS_DIR, "00_qc", "{sample}_kmer.json")
    params:
        crop_front = 15
    resources:
        cpus_per_task=10,
        mem_mb_per_cpu=2000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "fastp_kmer_{sample}.log")
    envmodules:
        "fastp/1.0.1-GCC-13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.json})
        mkdir -p $(dirname {log})
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --trim_front1 {params.crop_front} \
            --trim_front2 {params.crop_front} \
            --thread {resources.cpus_per_task} \
            -j {output.json} -h {output.html} &> {log}
        """

rule kmer_analysis:
    """
    Count k-mers (k=31) per sample using KMC3.
    Produces canonical (ci2) and all (ci0) k-mer databases.
    """
    input:
        r1 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1_kmer.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2_kmer.fq.gz")
    output:
        canon = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        all_kmers = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all.kmc_pre")
    params:
        outdir = os.path.join(RESULTS_DIR, "06_kmer", "{sample}")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=4000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "{sample}_kmer_analysis.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        cd {params.outdir}
        echo {input.r1} > input_files.txt
        echo {input.r2} >> input_files.txt
        {SOFTWARE_DIR}/kmc_v3 \
            -t{resources.cpus_per_task} -k31 -ci2 \
            @input_files.txt output_kmc_canon ./ \
            1> kmc_canon.1 2> kmc_canon.2
        {SOFTWARE_DIR}/kmc_v3 \
            -t{resources.cpus_per_task} -k31 -ci0 -b \
            @input_files.txt output_kmc_all ./ \
            1> kmc_all.1 2> kmc_all.2
        cat kmc_canon.2 kmc_all.2 > {log} 2>&1
        """


rule sex_specific_kmers:
    """
    Strictly sex-specific k-mers by set operations on the KMC databases.
    No association test: with n=8 there is no power, and "present in all
    males, absent in all females" is a set-membership question.

        male-specific   = intersect(males,   -ci5) - union(females, -ci1)
        female-specific = intersect(females, -ci5) - union(males,   -ci1)
        both then capped at -cx30 on the per-k-mer MINIMUM count.

    The thresholds are asymmetric on purpose: a k-mer must be solidly
    present (>=kmer_min_present) in every sample of its own sex, but a
    SINGLE read in one sample of the other sex disqualifies it.

    kmer_max_present drops repeat-derived k-mers. It is applied to the
    minimum count across the present group - never per sample, never to
    the subtraction side:
      - per sample, one outlier vetoes the cohort (6/7/8/36 across four
        males is dropped by the 36 alone, and on a Y that copy-number
        variation is normal);
      - on the subtraction side, a high-copy k-mer would go invisible in
        the other sex, fail to veto, and be reported as specific.
    Set it to 0 to disable.

    Evaluated as pairwise `kmc_tools simple` calls rather than one n-way
    `complex`, so only two databases are open at a time regardless of
    cohort size. The n-way form OOM-killed on 8 samples at 64 GB.

    Uses output_kmc_all (-ci0 -b) on both sides; output_kmc_canon is
    canonicalised and not comparable with it. Being a both-strands
    database, each k-mer and its reverse complement both appear: counts
    are per oriented k-mer, roughly half the locus depth.
    """
    input:
        male_dbs = expand(
            os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all.kmc_pre"),
            sample=ML_SAMPLES),
        female_dbs = expand(
            os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all.kmc_pre"),
            sample=FL_SAMPLES)
    output:
        male_kmers   = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "male_specific_kmers.txt"),
        female_kmers = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "female_specific_kmers.txt"),
        male_abyss   = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.input"),
        female_abyss = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.input")
    params:
        workdir     = os.path.join(RESULTS_DIR, "06_kmer", "combined", "kmc_sets"),
        db_dir      = os.path.join(RESULTS_DIR, "06_kmer"),
        males       = ML_SAMPLES,
        females     = FL_SAMPLES,
        min_present = config.get("kmer_min_present", 4),
        max_present = config.get("kmer_max_present", 30),
        min_absent  = config.get("kmer_min_absent", 1)
    resources:
        # Headroom, not a per-sample requirement: the pairwise form below
        # holds two databases open whatever the cohort size.
        cpus_per_task=6,
        mem_mb_per_cpu=10000,
        runtime=220
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "sex_specific_kmers.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir} $(dirname {output.male_kmers}) \
                 $(dirname {output.male_abyss}) $(dirname {log})

        MALES="{params.males}"
        FEMALES="{params.females}"

        # Upper bound is applied once, at the dump, to the minimum count
        # across the present group - see the docstring. 0 means no cap.
        CAP=""
        if [ "{params.max_present}" -gt 0 ]; then
            if [ "{params.max_present}" -le "{params.min_present}" ]; then
                echo "ERROR: kmer_max_present ({params.max_present}) must exceed" \
                     "kmer_min_present ({params.min_present}); the window is empty." >&2
                exit 1
            fi
            CAP="-cx{params.max_present}"
        fi

        cd {params.workdir}
        rm -f ./*_step*.kmc_pre ./*_step*.kmc_suf

        echo "=== group assignment ===" > {log}
        echo "males   ($MALES): $(echo $MALES | wc -w)"     >> {log}
        echo "females ($FEMALES): $(echo $FEMALES | wc -w)" >> {log}
        echo "=== thresholds ===" >> {log}
        echo "present, per sample:  -ci{params.min_present}" >> {log}
        echo "absent,  per sample:  -ci{params.min_absent}"  >> {log}
        echo "cap on min of present group: ${{CAP:-none}}"   >> {log}

        # Never removes anything but our own intermediates.
        drop_tmp () {{
            case "${{1:-}}" in
                *_step*) rm -f "$1.kmc_pre" "$1.kmc_suf" ;;
            esac
        }}

        # $1 present samples, $2 absent samples, $3 output db, $4 tag
        build_set () {{
            local present="$1" absent="$2" outdb="$3" tag="$4"
            local step=0 cur="" cur_params="" nxt="" prev=""

            # Fold the intersection two databases at a time. -ocmin keeps the
            # minimum counter, which is what the cap is later applied to.
            for s in $present; do
                if [ -z "$cur" ]; then
                    cur="{params.db_dir}/$s/output_kmc_all"
                    cur_params="-ci{params.min_present}"
                    continue
                fi
                step=$((step+1)); nxt="${{tag}}_step$step"
                echo "+ intersect [$cur $cur_params] x [$s -ci{params.min_present}] -> $nxt" >> {log}
                {KMC_TOOLS} -t{resources.cpus_per_task} simple \
                    "$cur" $cur_params \
                    "{params.db_dir}/$s/output_kmc_all" -ci{params.min_present} \
                    intersect "$nxt" -ci1 -ocmin >> {log} 2>&1
                drop_tmp "$prev"
                prev="$nxt"; cur="$nxt"; cur_params="-ci1"
            done

            # A - (B1+B2+B3+B4) == (((A-B1)-B2)-B3)-B4, so subtract one at a
            # time and never build the union, which would be larger than any
            # single database.
            for s in $absent; do
                step=$((step+1)); nxt="${{tag}}_step$step"
                echo "+ subtract [$s -ci{params.min_absent}] from [$cur] -> $nxt" >> {log}
                {KMC_TOOLS} -t{resources.cpus_per_task} simple \
                    "$cur" $cur_params \
                    "{params.db_dir}/$s/output_kmc_all" -ci{params.min_absent} \
                    kmers_subtract "$nxt" -ci1 >> {log} 2>&1
                drop_tmp "$prev"
                prev="$nxt"; cur="$nxt"; cur_params="-ci1"
            done

            # Refuse to move a sample database if no operation ever ran.
            case "$cur" in
                *_step*) ;;
                *) echo "ERROR: no set operation ran for $tag." >&2; exit 1 ;;
            esac
            mv "$cur.kmc_pre" "$outdb.kmc_pre"
            mv "$cur.kmc_suf" "$outdb.kmc_suf"
        }}

        echo "=== male set ===" >> {log}
        build_set "$MALES"   "$FEMALES" male_set   male
        echo "=== female set ===" >> {log}
        build_set "$FEMALES" "$MALES"   female_set female

        # -ci1 stops kmc_tools applying the database's own cutoff; $CAP is the
        # upper bound, hitting the minimum-across-present-group counter.
        {KMC_TOOLS} -t{resources.cpus_per_task} transform male_set   -ci1 $CAP dump {output.male_kmers}   >> {log} 2>&1
        {KMC_TOOLS} -t{resources.cpus_per_task} transform female_set -ci1 $CAP dump {output.female_kmers} >> {log} 2>&1

        awk '{{printf ">m%d\n%s\n", NR, $1}}' {output.male_kmers}   > {output.male_abyss}
        awk '{{printf ">f%d\n%s\n", NR, $1}}' {output.female_kmers} > {output.female_abyss}

        echo "=== counts ===" >> {log}
        echo "male-specific k-mers:   $(wc -l < {output.male_kmers})"   >> {log}
        echo "female-specific k-mers: $(wc -l < {output.female_kmers})" >> {log}
        overlap=$(LC_ALL=C comm -12 \
            <(cut -f1 {output.male_kmers}   | LC_ALL=C sort) \
            <(cut -f1 {output.female_kmers} | LC_ALL=C sort) | wc -l)
        echo "overlap (MUST be 0): $overlap" >> {log}
        if [ "$overlap" -ne 0 ]; then
            echo "ERROR: $overlap k-mers are in both sex-specific sets." >&2
            exit 1
        fi
        """


rule abyss_sex:
    """
    Assemble sex-specific k-mers with ABySS. {sex} is male or female.

    The input is the k-mer set itself, one 31-mer per FASTA record, so what
    limits contig formation here is how densely those k-mers tile the
    underlying region - not the assembler. Every k-mer lost upstream (to
    kmer_min_present, to the cap, or to a single stray read in the other
    sex) punches a hole in the tiling, and a hole wider than the assembly
    k splits one contig into two. If contigs come back short, the k-mer
    thresholds are usually the thing to revisit, not these parameters.

    Overridable from the config:
      abyss_k        assembly k, must be < 31. Lower bridges wider holes in
                     the tiling, giving longer and fewer-broken contigs, at
                     the cost of joining k-mers that share only abyss_k
                     bases, which can be chimeric. 25 is the historical
                     value; 21 is the first thing to try if the assembly is
                     fragmented.
      abyss_coverage -c, drop contigs below this mean k-mer coverage. The
                     counts on these k-mers are set-operation minima, not
                     assembly coverage, so filtering on them is not
                     meaningful - leave at 0.
      abyss_erode    -e, erode contig ends below this coverage. 0 = off.
      abyss_trim     -t, longest dangling edge to trim. 0 keeps tips as
                     contigs of their own rather than discarding them,
                     which is what you want when every input k-mer has
                     already passed the presence/absence filter.

    abyss-fac output is appended to the log so parameter sweeps can be
    compared directly (contig count, N50, max, total bp).
    """
    input:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "{sex}_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "{sex}_abyss.output")
    params:
        k        = config.get("abyss_k", 21),
        coverage = config.get("abyss_coverage", 0),
        erode    = config.get("abyss_erode", 0),
        trim     = config.get("abyss_trim", 0)
    wildcard_constraints:
        sex = "male|female"
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "abyss_{sex}.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        if [ -s {input} ]; then
            ABYSS -k{params.k} -c{params.coverage} -e{params.erode} -t{params.trim} \
                {input} -o {output} &> {log}
            echo "=== abyss-fac (-k{params.k} -c{params.coverage} -e{params.erode} -t{params.trim}) ===" >> {log}
            abyss-fac {output} >> {log} 2>&1 || true
        else
            echo "Input {input} is empty - no k-mers to assemble; writing empty output." > {log}
            touch {output}
        fi
        """


rule blast_db:
    """
    Create BLAST nucleotide database from reference genome.
    """
    input:
        ref = REFERENCE
    output:
        db = REFERENCE + ".nin"
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_db.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {log})
        makeblastdb -in {input.ref} -dbtype nucl 2> {log}
        """


rule blast_male_kmers:
    """
    BLAST male-assembled contigs against reference genome.
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.output"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "male_blast.out")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=15000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_male.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        if [ -s {input.contigs} ]; then
            blastn \
                -query {input.contigs} \
                -db {input.ref} \
                -task blastn \
                -dust yes \
                -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
                -evalue 1e-10 \
                -num_threads {resources.cpus_per_task} \
                > {output.blast} 2> {log}
        else
            echo "Input {input.contigs} is empty - no contigs to BLAST; writing empty output." > {log}
            touch {output.blast}
        fi
        """


rule blast_female_kmers:
    """
    BLAST female-assembled contigs against reference genome.
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.output"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "female_blast.out")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_female.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        if [ -s {input.contigs} ]; then
            blastn \
                -query {input.contigs} \
                -db {input.ref} \
                -task blastn \
                -dust yes \
                -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
                -evalue 1e-10 \
                -num_threads {resources.cpus_per_task} \
                > {output.blast} 2> {log}
        else
            echo "Input {input.contigs} is empty - no contigs to BLAST; writing empty output." > {log}
            touch {output.blast}
        fi
        """