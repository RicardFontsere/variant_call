# =============================================================================
# 11 - KMER COVERAGE: EXTRACT COUNTS FOR SEX-SPECIFIC & AUTOSOMAL KMERS
#
# Uses the PLINK association output from step 06 for sex-specific k-mers:
#   - Male-specific:   F_A=1, F_U=0 (present in all males, absent in all females)
#   - Female-specific: F_A=0, F_U=1 (absent in all males, present in all females)
#
# For autosomal/shared k-mers, we bypass the .assoc file entirely.
# K-mers present in ALL samples are monomorphic and never appear in PLINK
# output (filtered by --maf and --mac). Instead, we intersect all sample
# KMC canonical databases directly to find k-mers present in every sample.
# These are genuinely shared/autosomal and provide the diploid coverage
# baseline for comparison.
#
# Workflow:
# 1a. Sex-specific: extract k-mer sequences from .assoc → FASTA
# 1b. Autosomal: sequential KMC intersection of all samples → dump → subsample → FASTA
# 2.  Build KMC reference database from each FASTA
# 3.  Intersect ref db with each sample's canonical db (-ocleft for counts)
# 4.  Dump counts, combine into tables for plotting
# =============================================================================

AUTOSOMAL_SUBSAMPLE = config.get("autosomal_kmer_subsample", 10000000)


# =============================================================================
# STEP 1a: Sex-specific k-mers from PLINK association
# =============================================================================

rule male_kmers_to_fasta:
    """
    Extract male-specific k-mer sequences from PLINK association output.
    Uses the same filtered .assoc file produced by rule 06 (kmer_association).
    Explicit tab delimiter to match the sed-cleaned format.
    """
    input:
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "male_specific_kmers.assoc")
    output:
        fasta = os.path.join(RESULTS_DIR, "11_kmer_coverage", "fasta", "male_specific_kmers.fasta")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "male_kmers_to_fasta.log")
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})
        awk -F'\\t' '{{print ">male_kmer_" NR "\\n" $2}}' {input.assoc} > {output.fasta} 2> {log}
        echo "Male-specific k-mers: $(grep -c '^>' {output.fasta})" >> {log}
        """


rule female_kmers_to_fasta:
    """
    Extract female-specific k-mer sequences from PLINK association output.
    Uses the same filtered .assoc file produced by rule 06 (kmer_association).
    Explicit tab delimiter to match the sed-cleaned format.
    """
    input:
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "female_specific_kmers.assoc")
    output:
        fasta = os.path.join(RESULTS_DIR, "11_kmer_coverage", "fasta", "female_specific_kmers.fasta")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "female_kmers_to_fasta.log")
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})
        awk -F'\\t' '{{print ">female_kmer_" NR "\\n" $2}}' {input.assoc} > {output.fasta} 2> {log}
        echo "Female-specific k-mers: $(grep -c '^>' {output.fasta})" >> {log}
        """


# =============================================================================
# STEP 1b: Autosomal k-mers via KMC intersection of all samples
# =============================================================================

rule autosomal_kmer_intersection:
    """
    Find k-mers present in ALL samples by sequentially intersecting
    every sample's canonical KMC database.

    kmc_tools intersect operates on two databases at a time, so we chain:
        intersect(sample1, sample2) → tmp_2
        intersect(tmp_2, sample3)   → tmp_3
        ...
        intersect(tmp_N-1, sampleN) → final

    The result contains k-mers found in every single sample — truly shared
    (autosomal) k-mers that were invisible to PLINK due to being monomorphic.

    We then dump the full set, subsample, and convert to FASTA.
    """
    input:
        kmc_dbs = expand(
            os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_pre"),
            sample=SAMPLES
        )
    output:
        fasta = os.path.join(RESULTS_DIR, "11_kmer_coverage", "fasta", "autosomal_kmers.fasta")
    params:
        sample_db_prefixes = [
            os.path.join(RESULTS_DIR, "06_kmer", s, "output_kmc_canon") for s in SAMPLES
        ],
        workdir = os.path.join(RESULTS_DIR, "11_kmer_coverage", "autosomal_tmp"),
        n_subsample = AUTOSOMAL_SUBSAMPLE
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=16000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "autosomal_kmer_intersection.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.workdir}
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})

        SAMPLES=({params.sample_db_prefixes})
        N=${{#SAMPLES[@]}}
        echo "Intersecting $N sample databases" > {log}

        # Start with the first sample
        CURRENT="${{SAMPLES[0]}}"
        echo "  Starting with: $CURRENT" >> {log}

        # Sequentially intersect with each remaining sample
        for (( i=1; i<N; i++ )); do
            NEXT="${{SAMPLES[$i]}}"
            OUT="{params.workdir}/tmp_$i"
            echo "  Intersecting with: $NEXT -> $OUT" >> {log}

            /user/brussel/109/vsc10945/home/scratch/Software/kmc/bin/kmc_tools simple \
                "$CURRENT" \
                "$NEXT" \
                intersect \
                "$OUT" >> {log} 2>&1

            # Clean up previous temp (not the original sample dbs!)
            if [ $i -gt 1 ]; then
                PREV="{params.workdir}/tmp_$((i-1))"
                rm -f "${{PREV}}.kmc_pre" "${{PREV}}.kmc_suf"
            fi

            CURRENT="$OUT"
        done

        # Dump the final intersection
        DUMP_FILE="{params.workdir}/all_shared_kmers.txt"
        echo "Dumping shared k-mers..." >> {log}
        /user/brussel/109/vsc10945/home/scratch/Software/kmc/bin/kmc_tools transform \
            "$CURRENT" \
            dump "$DUMP_FILE" >> {log} 2>&1

        TOTAL=$(wc -l < "$DUMP_FILE")
        echo "Total k-mers present in ALL samples: $TOTAL" >> {log}

        # Subsample and convert to FASTA
        if [ "$TOTAL" -le "{params.n_subsample}" ]; then
            echo "Using all $TOTAL shared k-mers (below subsample limit)" >> {log}
            awk '{{print ">auto_kmer_" NR "\\n" $1}}' "$DUMP_FILE" > {output.fasta}
        else
            echo "Subsampling {params.n_subsample} from $TOTAL shared k-mers" >> {log}
            shuf -n {params.n_subsample} "$DUMP_FILE" | \
                awk '{{print ">auto_kmer_" NR "\\n" $1}}' > {output.fasta}
        fi

        echo "Autosomal k-mers in FASTA: $(grep -c '^>' {output.fasta})" >> {log}

        # Clean up
        rm -f "${{CURRENT}}.kmc_pre" "${{CURRENT}}.kmc_suf"
        rm -f "$DUMP_FILE"
        rmdir {params.workdir} 2>/dev/null || true
        """


# =============================================================================
# STEP 2: Build KMC reference databases from FASTA
# =============================================================================

rule build_kmer_ref_db:
    """
    Build a KMC database from a k-mer FASTA file.
    Serves as the reference set for intersection.
    -ci1: keep all (each k-mer appears once in the FASTA)
    -fm:  input is FASTA format
    -k31: same k-mer size as the sample databases

    Wildcard {category}: male_specific, female_specific, or autosomal
    """
    input:
        fasta = os.path.join(RESULTS_DIR, "11_kmer_coverage", "fasta", "{category}_kmers.fasta")
    output:
        kmc_pre = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}.kmc_pre"),
        kmc_suf = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}.kmc_suf")
    params:
        db_prefix = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}"),
        tmpdir = os.path.join(RESULTS_DIR, "11_kmer_coverage", "tmp", "{category}")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "build_db_{category}.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.tmpdir}
        mkdir -p $(dirname {log})
        {KMERSGWAS_DIR}/external_programs/kmc_v3 \
            -t{resources.cpus_per_task} \
            -k31 \
            -fm \
            -ci1 \
            {input.fasta} \
            {params.db_prefix} \
            {params.tmpdir} &> {log}
        rm -rf {params.tmpdir}
        """


# =============================================================================
# STEP 3: Intersect with sample KMC databases and dump counts
# =============================================================================

rule intersect_kmer_counts:
    """
    Intersect a sample's canonical KMC database with a k-mer reference database.
    -ocleft: keep counts from the sample (left) database.

    Without -ocleft, kmc_tools takes min(sample_count, ref_count).
    Since the ref db was built with -ci1, every ref k-mer has count=1,
    and min(anything, 1) = 1 — giving wrong results.

    Wildcards:
        {category}: male_specific, female_specific, or autosomal
        {sample}:   sample name
    """
    input:
        sample_pre = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        sample_suf = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_suf"),
        ref_pre = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}.kmc_pre"),
        ref_suf = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}.kmc_suf")
    output:
        dump = os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "{category}", "{sample}.txt")
    params:
        sample_db = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon"),
        ref_db = os.path.join(RESULTS_DIR, "11_kmer_coverage", "ref_db", "{category}"),
        intersect_db = os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "{category}", "{sample}_intersect")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "intersect_{category}_{sample}.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.dump})
        mkdir -p $(dirname {log})

        # Intersect: -ocleft keeps counts from the sample (left) database
        /user/brussel/109/vsc10945/home/scratch/Software/kmc/bin/kmc_tools simple \
            {params.sample_db} \
            {params.ref_db} \
            intersect \
            {params.intersect_db} \
            -ocleft &> {log}

        # Dump to text: kmer<tab>count
        /user/brussel/109/vsc10945/home/scratch/Software/kmc/bin/kmc_tools transform \
            {params.intersect_db} \
            dump {output.dump} >> {log} 2>&1

        # Clean up intermediate db
        rm -f {params.intersect_db}.kmc_pre {params.intersect_db}.kmc_suf

        echo "K-mers with counts: $(wc -l < {output.dump})" >> {log}
        """


# =============================================================================
# STEP 4: Combine per-sample counts into summary tables
# =============================================================================

rule combine_male_kmer_counts:
    """
    Combine male-specific and autosomal k-mer counts across all male samples.
    Output: sample, category, count — for histogram plotting.
    """
    input:
        male_dumps = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "male_specific", "{sample}.txt"),
            sample=ML_SAMPLES
        ),
        auto_dumps = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "autosomal", "{sample}.txt"),
            sample=ML_SAMPLES
        )
    output:
        table = os.path.join(RESULTS_DIR, "11_kmer_coverage", "results", "male_samples_counts.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "combine_male_counts.log")
    shell:
        """
        mkdir -p $(dirname {output.table})
        mkdir -p $(dirname {log})

        echo -e "sample\\tcategory\\tcount" > {output.table}

        # Male-specific k-mer counts
        for f in {input.male_dumps}; do
            sample=$(basename "$f" .txt)
            awk -v s="$sample" -v c="male_specific" 'BEGIN{{OFS="\\t"}} {{print s, c, $2}}' "$f" \
                >> {output.table}
        done

        # Autosomal k-mer counts in male samples
        for f in {input.auto_dumps}; do
            sample=$(basename "$f" .txt)
            awk -v s="$sample" -v c="autosomal" 'BEGIN{{OFS="\\t"}} {{print s, c, $2}}' "$f" \
                >> {output.table}
        done

        echo "Male-specific rows: $(grep -c 'male_specific' {output.table})" > {log}
        echo "Autosomal rows:     $(grep -c 'autosomal' {output.table})" >> {log}
        """


rule combine_female_kmer_counts:
    """
    Combine female-specific and autosomal k-mer counts across all female samples.
    Output: sample, category, count — for histogram plotting.
    """
    input:
        female_dumps = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "female_specific", "{sample}.txt"),
            sample=FL_SAMPLES
        ),
        auto_dumps = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "autosomal", "{sample}.txt"),
            sample=FL_SAMPLES
        )
    output:
        table = os.path.join(RESULTS_DIR, "11_kmer_coverage", "results", "female_samples_counts.tsv")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "combine_female_counts.log")
    shell:
        """
        mkdir -p $(dirname {output.table})
        mkdir -p $(dirname {log})

        echo -e "sample\\tcategory\\tcount" > {output.table}

        # Female-specific k-mer counts
        for f in {input.female_dumps}; do
            sample=$(basename "$f" .txt)
            awk -v s="$sample" -v c="female_specific" 'BEGIN{{OFS="\\t"}} {{print s, c, $2}}' "$f" \
                >> {output.table}
        done

        # Autosomal k-mer counts in female samples
        for f in {input.auto_dumps}; do
            sample=$(basename "$f" .txt)
            awk -v s="$sample" -v c="autosomal" 'BEGIN{{OFS="\\t"}} {{print s, c, $2}}' "$f" \
                >> {output.table}
        done

        echo "Female-specific rows: $(grep -c 'female_specific' {output.table})" > {log}
        echo "Autosomal rows:       $(grep -c 'autosomal' {output.table})" >> {log}
        """


# =============================================================================
# STEP 5: Coverage-based filtering of sex-specific k-mers
#
# Male-specific k-mers from PLINK (F_A=1, F_U=0) include both genuine
# Y-linked k-mers (hemizygous) and autosomal artifacts (diploid).
# We filter by requiring the count to fall in the hemizygous range [3,11]
# in ALL samples of the relevant sex.
#
# Y-mers: male-specific k-mers with count in [3,11] in ALL males
# W-mers: female-specific k-mers with count in [3,11] in ALL females
# =============================================================================

KMER_COUNT_LOWER = config.get("kmer_count_lower", 3)
KMER_COUNT_UPPER = config.get("kmer_count_upper", 11)


rule filter_ymers:
    """
    Filter male-specific k-mers to retain only those with hemizygous
    coverage (count in [LOWER, UPPER]) in ALL male samples.

    Reads per-sample count dumps, joins by k-mer, keeps only k-mers
    where every sample's count is in range. Then filters the original
    .assoc file to keep matching rows.
    """
    input:
        count_files = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "male_specific", "{sample}.txt"),
            sample=ML_SAMPLES
        ),
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "male_specific_kmers.assoc")
    output:
        assoc = os.path.join(RESULTS_DIR, "11_kmer_coverage", "filtered", "ymers.assoc"),
        discarded = os.path.join(RESULTS_DIR, "11_kmer_coverage", "filtered", "ymers_discarded.assoc")
    params:
        lower = KMER_COUNT_LOWER,
        upper = KMER_COUNT_UPPER,
        n_samples = len(ML_SAMPLES)
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "filter_ymers.log")
    shell:
        """
        mkdir -p $(dirname {output.assoc})
        mkdir -p $(dirname {log})

        # Build list of k-mers that pass: count in [lower,upper] in ALL males
        # Each count file is: kmer<tab>count
        # Step 1: cat all files, tag each line, then use awk to group by kmer
        cat {input.count_files} | \
        awk -v lo={params.lower} -v hi={params.upper} -v n={params.n_samples} '
        {{
            if ($2 >= lo && $2 <= hi)
                pass[$1]++
        }}
        END {{
            for (k in pass)
                if (pass[k] == n)
                    print k
        }}' | sort > {output.assoc}.pass_kmers 2> {log}

        TOTAL_IN=$(wc -l < {input.assoc})
        PASS=$(wc -l < {output.assoc}.pass_kmers)
        echo "Input k-mers: $TOTAL_IN" >> {log}
        echo "Passing hemizygous filter [{params.lower},{params.upper}] in all {params.n_samples} males: $PASS" >> {log}

        # Filter the .assoc file: keep rows where column 2 (kmer) is in pass set
        awk -F'\\t' 'NR==FNR {{keep[$1]; next}} ($2 in keep)' \
            {output.assoc}.pass_kmers {input.assoc} > {output.assoc} 2>> {log}

        # Write discarded for inspection
        awk -F'\\t' 'NR==FNR {{keep[$1]; next}} !($2 in keep)' \
            {output.assoc}.pass_kmers {input.assoc} > {output.discarded} 2>> {log}

        echo "Y-mers (kept): $(wc -l < {output.assoc})" >> {log}
        echo "Discarded: $(wc -l < {output.discarded})" >> {log}

        rm -f {output.assoc}.pass_kmers
        """


rule filter_wmers:
    """
    Filter female-specific k-mers to retain only those with hemizygous
    coverage (count in [LOWER, UPPER]) in ALL female samples.

    In an XY system, W-mers are not expected. Any that pass this filter
    are worth investigating — they could indicate X-linked heterozygous
    variants or structural variation.
    """
    input:
        count_files = expand(
            os.path.join(RESULTS_DIR, "11_kmer_coverage", "counts", "female_specific", "{sample}.txt"),
            sample=FL_SAMPLES
        ),
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "female_specific_kmers.assoc")
    output:
        assoc = os.path.join(RESULTS_DIR, "11_kmer_coverage", "filtered", "wmers.assoc"),
        discarded = os.path.join(RESULTS_DIR, "11_kmer_coverage", "filtered", "wmers_discarded.assoc")
    params:
        lower = KMER_COUNT_LOWER,
        upper = KMER_COUNT_UPPER,
        n_samples = len(FL_SAMPLES)
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=16000,
        runtime=30
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "filter_wmers.log")
    shell:
        """
        mkdir -p $(dirname {output.assoc})
        mkdir -p $(dirname {log})

        # Build list of k-mers that pass: count in [lower,upper] in ALL females
        cat {input.count_files} | \
        awk -v lo={params.lower} -v hi={params.upper} -v n={params.n_samples} '
        {{
            if ($2 >= lo && $2 <= hi)
                pass[$1]++
        }}
        END {{
            for (k in pass)
                if (pass[k] == n)
                    print k
        }}' | sort > {output.assoc}.pass_kmers 2> {log}

        TOTAL_IN=$(wc -l < {input.assoc})
        PASS=$(wc -l < {output.assoc}.pass_kmers)
        echo "Input k-mers: $TOTAL_IN" >> {log}
        echo "Passing hemizygous filter [{params.lower},{params.upper}] in all {params.n_samples} females: $PASS" >> {log}

        # Filter the .assoc file
        awk -F'\\t' 'NR==FNR {{keep[$1]; next}} ($2 in keep)' \
            {output.assoc}.pass_kmers {input.assoc} > {output.assoc} 2>> {log}

        # Write discarded for inspection
        awk -F'\\t' 'NR==FNR {{keep[$1]; next}} !($2 in keep)' \
            {output.assoc}.pass_kmers {input.assoc} > {output.discarded} 2>> {log}

        echo "W-mers (kept): $(wc -l < {output.assoc})" >> {log}
        echo "Discarded: $(wc -l < {output.discarded})" >> {log}

        rm -f {output.assoc}.pass_kmers
        """


# =============================================================================
# STEP 6: Assemble and BLAST filtered k-mers
#
# Same pipeline as rule 06 (ABySS → BLAST) but using coverage-filtered
# Y-mers and W-mers instead of the raw sex-specific k-mers.
# =============================================================================

rule filtered_to_abyss_input:
    """
    Convert filtered .assoc to ABySS FASTA input.
    Format: >pvalue kmer_sequence
    Column 2 = k-mer, column 9 = p-value (0-indexed: $2 and $9 in the .assoc).
    """
    input:
        assoc = os.path.join(RESULTS_DIR, "11_kmer_coverage", "filtered", "{sex}mers.assoc")
    output:
        abyss_input = os.path.join(RESULTS_DIR, "11_kmer_coverage", "assembly", "{sex}mers_abyss.input")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=2000,
        runtime=10
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "{sex}mers_to_abyss.log")
    shell:
        """
        mkdir -p $(dirname {output.abyss_input})
        mkdir -p $(dirname {log})
        awk -F'\\t' '{{print ">" $9 "\\n" $2}}' {input.assoc} > {output.abyss_input} 2> {log}
        echo "K-mers for assembly: $(grep -c '^>' {output.abyss_input})" >> {log}
        """


rule abyss_filtered:
    """
    Assemble coverage-filtered sex-specific k-mers with ABySS.
    sex wildcard: y or w
    """
    input:
        os.path.join(RESULTS_DIR, "11_kmer_coverage", "assembly", "{sex}mers_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "11_kmer_coverage", "assembly", "{sex}mers_abyss.output")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "abyss_{sex}mers.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        ABYSS -k25 -c0 -e0 {input} -o {output} &> {log}
        """


rule blast_filtered:
    """
    BLAST coverage-filtered assembled contigs against reference genome.
    Reuses the BLAST database created in step 06.
    sex wildcard: y or w
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "11_kmer_coverage", "assembly", "{sex}mers_abyss.output"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "11_kmer_coverage", "blast", "{sex}mers_blast.out")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=15000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "11_kmer_coverage", "blast_{sex}mers.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        blastn -query {input.contigs} -db {input.ref} \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
            -evalue 1e-5 \
            -num_threads {resources.cpus_per_task} \
            > {output.blast} 2> {log}
        """
