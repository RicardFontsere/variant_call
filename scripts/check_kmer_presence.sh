#!/usr/bin/env bash
# Quick presence/absence check for a set of k-mers against the per-sample KMC
# canonical databases from rule 06 (k=31, canonical, -ci2). Same k-mer
# definition the pipeline used, so a hit here means the same thing a hit meant
# in the association test.
#
#   scripts/check_kmer_presence.sh -d results/06_kmer -i male_specific_kmers.assoc
#
# Writes to the output dir:
#   summary.tsv          sample, sex, n_present, n_total, fraction
#   presence_matrix.tsv  k-mer x sample, 0/1
# and prints the sex breakdown to stdout.
#
# kmc / kmc_tools are taken from $KMC and $KMC_TOOLS, else from PATH.

set -euo pipefail

KMER_DIR=""
KMER_FILE=""
OUT="kmer_check"
N=1000
MALE_PAT="ML"
FEMALE_PAT="FL"
MIN_COUNT=1
THREADS=4

usage() {
    cat <<EOF
Usage: $(basename "$0") -d <kmer_dir> -i <kmer_file> [options]

  -d DIR    directory holding per-sample subdirs with output_kmc_canon.* (rule 06)
  -i FILE   k-mer source: a .assoc from rule 06, or a plain list of 31-mers
  -o DIR    output directory (default: $OUT)
  -n INT    max k-mers to test (default: $N; 0 = all)
  -m PAT    male sample-name pattern (default: $MALE_PAT)
  -f PAT    female sample-name pattern (default: $FEMALE_PAT)
  -c INT    min k-mer count to call present (default: $MIN_COUNT)
  -t INT    threads (default: $THREADS)
EOF
}

while getopts "d:i:o:n:m:f:c:t:h" opt; do
    case "$opt" in
        d) KMER_DIR=$OPTARG ;;
        i) KMER_FILE=$OPTARG ;;
        o) OUT=$OPTARG ;;
        n) N=$OPTARG ;;
        m) MALE_PAT=$OPTARG ;;
        f) FEMALE_PAT=$OPTARG ;;
        c) MIN_COUNT=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

[ -n "$KMER_DIR" ] && [ -n "$KMER_FILE" ] || { usage >&2; exit 2; }
[ -d "$KMER_DIR" ] || { echo "No such directory: $KMER_DIR" >&2; exit 1; }
[ -s "$KMER_FILE" ] || { echo "Missing or empty k-mer file: $KMER_FILE" >&2; exit 1; }

KMC=${KMC:-kmc}
KMC_TOOLS=${KMC_TOOLS:-kmc_tools}
command -v "$KMC" >/dev/null || { echo "kmc not found (set \$KMC)" >&2; exit 1; }
command -v "$KMC_TOOLS" >/dev/null || { echo "kmc_tools not found (set \$KMC_TOOLS)" >&2; exit 1; }

mkdir -p "$OUT/tmp" "$OUT/dump"
LOG="$OUT/kmc.log"
: > "$LOG"

# ---------------------------------------------------------------------------
# 1. k-mers -> FASTA. Takes the first 31-mer looking field on each line, which
#    covers both a bare list and the SNP column of a PLINK .assoc.
# ---------------------------------------------------------------------------
awk '
    { for (i = 1; i <= NF; i++)
          if (length($i) >= 31 && substr($i, 1, 31) ~ /^[ACGT]+$/) {
              print substr($i, 1, 31); break
          } }
' "$KMER_FILE" | sort -u | awk -v n="$N" 'n == 0 || NR <= n' > "$OUT/kmers.txt"

[ -s "$OUT/kmers.txt" ] || { echo "No 31-mers found in $KMER_FILE" >&2; exit 1; }
awk '{ print ">k" NR "\n" $0 }' "$OUT/kmers.txt" > "$OUT/query.fa"

# ---------------------------------------------------------------------------
# 2. Build the query database. -ci1 keeps every k-mer (each appears once);
#    KMC canonicalises it here, exactly as it did for the sample databases.
# ---------------------------------------------------------------------------
"$KMC" -t"$THREADS" -k31 -fm -ci1 "$OUT/query.fa" "$OUT/query_db" "$OUT/tmp" >> "$LOG" 2>&1

# Canonical forms are the row keys from here on: the intersect dumps report
# canonical k-mers, which need not match the input orientation.
"$KMC_TOOLS" transform "$OUT/query_db" dump "$OUT/query_dump.txt" >> "$LOG" 2>&1
cut -f1 "$OUT/query_dump.txt" | sort > "$OUT/kmers.canon.txt"
TOTAL=$(wc -l < "$OUT/kmers.canon.txt")

# ---------------------------------------------------------------------------
# 3. Intersect against each sample. -ocleft keeps the sample's counts; without
#    it kmc_tools takes min(sample, query) and every count collapses to 1.
# ---------------------------------------------------------------------------
: > "$OUT/samples.txt"
: > "$OUT/long.tsv"

shopt -s nullglob
DBS=("$KMER_DIR"/*/output_kmc_canon.kmc_pre)
shopt -u nullglob
[ ${#DBS[@]} -gt 0 ] || { echo "No output_kmc_canon.kmc_pre under $KMER_DIR" >&2; exit 1; }

for pre in "${DBS[@]}"; do
    sample=$(basename "$(dirname "$pre")")
    echo "$sample" >> "$OUT/samples.txt"

    "$KMC_TOOLS" simple "${pre%.kmc_pre}" "$OUT/query_db" \
        intersect "$OUT/tmp/$sample" -ocleft >> "$LOG" 2>&1
    "$KMC_TOOLS" transform "$OUT/tmp/$sample" dump "$OUT/dump/$sample.txt" >> "$LOG" 2>&1
    rm -f "$OUT/tmp/$sample.kmc_pre" "$OUT/tmp/$sample.kmc_suf"

    awk -v s="$sample" -v c="$MIN_COUNT" -v OFS='\t' \
        '$2 >= c { print $1, s }' "$OUT/dump/$sample.txt" >> "$OUT/long.tsv"
done

# ---------------------------------------------------------------------------
# 4. Summarise.
# ---------------------------------------------------------------------------
awk -v OFS='\t' -v total="$TOTAL" -v mp="$MALE_PAT" -v fp="$FEMALE_PAT" '
    NR == FNR { order[++n] = $1; next }
    { count[$2]++ }
    END {
        print "sample", "sex", "n_present", "n_total", "fraction"
        for (i = 1; i <= n; i++) {
            s = order[i]
            sex = (s ~ mp) ? "M" : (s ~ fp) ? "F" : "?"
            c = count[s] + 0
            printf "%s\t%s\t%d\t%d\t%.4f\n", s, sex, c, total, total ? c / total : 0
        }
    }
' "$OUT/samples.txt" "$OUT/long.tsv" > "$OUT/summary.tsv"

awk -v OFS='\t' '
    FILENAME == ARGV[1] { samples[++ns] = $1; next }
    FILENAME == ARGV[2] { kmers[++nk] = $1; next }
    { present[$1 SUBSEP $2] = 1 }
    END {
        printf "kmer"
        for (j = 1; j <= ns; j++) printf "\t%s", samples[j]
        printf "\n"
        for (i = 1; i <= nk; i++) {
            printf "%s", kmers[i]
            for (j = 1; j <= ns; j++)
                printf "\t%d", (kmers[i] SUBSEP samples[j]) in present ? 1 : 0
            printf "\n"
        }
    }
' "$OUT/samples.txt" "$OUT/kmers.canon.txt" "$OUT/long.tsv" > "$OUT/presence_matrix.tsv"

echo "Tested $TOTAL k-mers against ${#DBS[@]} samples (min count $MIN_COUNT)"
echo
awk -F'\t' 'NR > 1 { printf "  %-20s %s  %5d/%-5d  %.3f\n", $1, $2, $3, $4, $5 }' "$OUT/summary.tsv"
echo
awk -F'\t' -v total="$TOTAL" '
    NR > 1 && $2 == "M" { m++; mp += $3; if ($3 < total) mbad++ }
    NR > 1 && $2 == "F" { f++; fp += $3; if ($3 > 0)     fbad++ }
    END {
        printf "males:   %d samples, mean %.1f/%d present, %d with a missing k-mer\n",
               m, m ? mp / m : 0, total, mbad + 0
        printf "females: %d samples, mean %.1f/%d present, %d with any k-mer present\n",
               f, f ? fp / f : 0, total, fbad + 0
        if (mbad + fbad == 0)
            print "\nClean: present in every male, absent in every female."
        else
            print "\nThe filter does not hold on the raw databases - see summary.tsv."
    }
' "$OUT/summary.tsv"

echo
echo "Wrote $OUT/summary.tsv and $OUT/presence_matrix.tsv"
