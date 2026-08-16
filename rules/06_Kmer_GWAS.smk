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
    Extract strictly sex-specific k-mers by set operations on the KMC
    databases (no association test - with n=8 there is no power anyway,
    and this is a deterministic set membership question).

        male-specific   = (intersection of all males, -ci5 -cx30)
                          - (union of all females, -ci1)
        female-specific = (intersection of all females, -ci5 -cx30)
                          - (union of all males, -ci1)

    Asymmetric -ci thresholds are deliberate: a k-mer must be solidly
    present (>=kmer_min_present copies) to count as present, but a SINGLE
    read is enough to count as present on the subtraction side. This is
    what keeps low-coverage samples from producing false "specific" k-mers.

    -cx (kmer_max_present) caps the count on the intersection side to drop
    satellite/repeat-derived k-mers, which sail through the intersection
    and pollute the assembly. It is applied to the PRESENT side ONLY, and
    deliberately never to the subtraction side: an -cx there would make a
    high-copy k-mer invisible in the other sex, so it would fail to veto
    and would be reported as sex-specific - the exact opposite of the
    intent. Set kmer_max_present to 0 to disable the cap.

    Both bounds are counts of a single ORIENTED k-mer in a -b database
    (see below), i.e. roughly half the locus depth. Budget accordingly if
    you ever switch these rules to the canonical database.

    Uses output_kmc_all (built -ci0 -b) on BOTH sides. Do not mix with
    output_kmc_canon: that database is canonicalised and the k-mer
    representations are not comparable.

    Because output_kmc_all is a both-strands (-b) database, a sex-specific
    k-mer and its reverse complement both appear in the output. That is
    self-consistent and harmless for ABySS, but it means the reported
    counts are ~2x the number of distinct canonical k-mers.
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
        min_present = config.get("kmer_min_present", 5),
        max_present = config.get("kmer_max_present", 30),
        min_absent  = config.get("kmer_min_absent", 1)
    resources:
        cpus_per_task=8,
        mem_mb_per_cpu=8000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "sex_specific_kmers.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir} $(dirname {output.male_kmers}) \
                 $(dirname {output.male_abyss}) $(dirname {log})

        # Resolve everything to absolute paths BEFORE cd: {{log}} and the
        # outputs are relative to the workflow root, not to params.workdir.
        DB=$(realpath -m {params.db_dir})
        LOG=$(realpath -m {log})
        MALE_TXT=$(realpath -m {output.male_kmers})
        FEMALE_TXT=$(realpath -m {output.female_kmers})
        MALE_FA=$(realpath -m {output.male_abyss})
        FEMALE_FA=$(realpath -m {output.female_abyss})
        MALE_DB=$(realpath -m {params.workdir})/male_set
        FEMALE_DB=$(realpath -m {params.workdir})/female_set

        MALES="{params.males}"
        FEMALES="{params.females}"
        if [ -z "$MALES" ] || [ -z "$FEMALES" ]; then
            echo "ERROR: one of the sex groups is empty (males='$MALES' females='$FEMALES')." \
                 "Check male_pattern/female_pattern in the config." >&2
            exit 1
        fi

        # Count window for the intersection (present) side. -cx0 would mean
        # "no k-mer may exceed 0 copies", so 0 is treated as "no cap" and the
        # flag is omitted entirely rather than passed as -cx0.
        PRESENT_PARAMS="-ci{params.min_present}"
        if [ "{params.max_present}" -gt 0 ]; then
            if [ "{params.max_present}" -le "{params.min_present}" ]; then
                echo "ERROR: kmer_max_present ({params.max_present}) must be greater than" \
                     "kmer_min_present ({params.min_present}); the window is empty." >&2
                exit 1
            fi
            PRESENT_PARAMS="$PRESENT_PARAMS -cx{params.max_present}"
        fi

        cd {params.workdir}

        echo "=== group assignment ===" > "$LOG"
        echo "males   ($MALES): $(echo $MALES | wc -w)"     >> "$LOG"
        echo "females ($FEMALES): $(echo $FEMALES | wc -w)" >> "$LOG"
        echo "=== thresholds ===" >> "$LOG"
        echo "present side (intersection): $PRESENT_PARAMS"  >> "$LOG"
        echo "absent side  (subtraction):  -ci{params.min_absent}" >> "$LOG"

        # ---- build ops files -------------------------------------------
        build_ops () {{
            # $1 = present-group samples, $2 = absent-group samples
            # $3 = output database path,  $4 = ops file to write
            local present="$1" absent="$2" outdb="$3" ops="$4"
            local i=0 p_terms="" a_terms=""
            echo "INPUT:" > "$ops"
            for s in $present; do
                i=$((i+1))
                echo "p$i = $DB/$s/output_kmc_all $PRESENT_PARAMS" >> "$ops"
                p_terms="${{p_terms:+$p_terms*}}p$i"
            done
            i=0
            for s in $absent; do
                i=$((i+1))
                echo "a$i = $DB/$s/output_kmc_all -ci{params.min_absent}" >> "$ops"
                a_terms="${{a_terms:+$a_terms+}}a$i"
            done
            echo "OUTPUT:" >> "$ops"
            # The name left of '=' is the output database path, so each run
            # writes straight to its own db - no shared result.kmc_* to rename.
            echo "$outdb = ($p_terms)-($a_terms)" >> "$ops"
            echo "OUTPUT_PARAMS:" >> "$ops"
            echo "-ci1" >> "$ops"
        }}

        build_ops "$MALES"   "$FEMALES" "$MALE_DB"   ops_male.txt
        build_ops "$FEMALES" "$MALES"   "$FEMALE_DB" ops_female.txt

        echo "=== ops_male.txt ==="   >> "$LOG"; cat ops_male.txt   >> "$LOG"
        echo "=== ops_female.txt ===" >> "$LOG"; cat ops_female.txt >> "$LOG"

        # ---- run set operations ----------------------------------------
        {KMC_TOOLS} -t{resources.cpus_per_task} complex ops_male.txt   >> "$LOG" 2>&1
        {KMC_TOOLS} -t{resources.cpus_per_task} complex ops_female.txt >> "$LOG" 2>&1

        # -ci1 on the dump side too: kmc_tools would otherwise apply the
        # database's own default cutoff and silently drop singleton k-mers.
        {KMC_TOOLS} -t{resources.cpus_per_task} transform "$MALE_DB"   -ci1 dump "$MALE_TXT"   >> "$LOG" 2>&1
        {KMC_TOOLS} -t{resources.cpus_per_task} transform "$FEMALE_DB" -ci1 dump "$FEMALE_TXT" >> "$LOG" 2>&1

        # ---- FASTA for ABySS -------------------------------------------
        awk '{{printf ">m%d\n%s\n", NR, $1}}' "$MALE_TXT"   > "$MALE_FA"
        awk '{{printf ">f%d\n%s\n", NR, $1}}' "$FEMALE_TXT" > "$FEMALE_FA"

        # ---- sanity ----------------------------------------------------
        echo "=== counts ===" >> "$LOG"
        echo "male-specific k-mers:   $(wc -l < "$MALE_TXT")"   >> "$LOG"
        echo "female-specific k-mers: $(wc -l < "$FEMALE_TXT")" >> "$LOG"
        overlap=$(LC_ALL=C comm -12 \
            <(cut -f1 "$MALE_TXT"   | LC_ALL=C sort) \
            <(cut -f1 "$FEMALE_TXT" | LC_ALL=C sort) | wc -l)
        echo "overlap (MUST be 0): $overlap" >> "$LOG"
        if [ "$overlap" -ne 0 ]; then
            echo "ERROR: $overlap k-mers are in both sex-specific sets." >&2
            exit 1
        fi
        """


rule abyss_male:
    """
    Assemble male-specific k-mers with ABySS.
    """
    input:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.output")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "abyss_male.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        if [ -s {input} ]; then
            ABYSS -k25 -c0 -e0 {input} -o {output} &> {log}
        else
            echo "Input {input} is empty - ...; writing empty output." > {log}
            touch {output}
        fi
        """


rule abyss_female:
    """
    Assemble female-specific k-mers with ABySS.
    """
    input:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.output")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "abyss_female.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        if [ -s {input} ]; then
            ABYSS -k25 -c0 -e0 {input} -o {output} &> {log}
        else
            echo "Input {input} is empty - ...; writing empty output." > {log}
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
                -perc_identity 80 \
                -qcov_hsp_perc 80 \
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
                -perc_identity 80 \
                -qcov_hsp_perc 80 \
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