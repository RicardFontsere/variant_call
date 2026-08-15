#!/usr/bin/env bash
# =============================================================================
# ONE-OFF: verify that male_specific_kmers.assoc really contains male k-mers
#
# The chain being audited (rules/06_Kmer_GWAS.smk):
#   phenotype.pheno   males -> 1, females -> 2          (kmers_table_to_bed)
#   kmers_to_plink.0  .fam phenotype column             (kmers_table_to_bed)
#   kmers.assoc       F_A / F_U columns                 (plink --assoc)
#   awk $5==1 && $6==0  -> male_specific_kmers.assoc    (kmer_association)
#
# Two conventions have to line up for that awk to mean what the comment says,
# and neither is checked anywhere:
#   1. PLINK reads phenotype 1 as CONTROL and 2 as CASE. Males are coded 1,
#      so F_A is the FEMALE column and F_U is the MALE column -- the opposite
#      of what the rule assumes. If that is what happened, the two output
#      files are simply swapped.
#   2. F_A/F_U are frequencies of allele A1. If A1 is the "k-mer absent"
#      allele rather than "k-mer present", the test inverts again.
#
# Two inversions cancel, so guessing is not safe -- this measures instead.
# Steps 1-6 are mechanical. Step 7 prints commands for an independent
# ground-truth check straight against the raw FASTQs, which depends on no
# PLINK or kmersGWAS convention at all.
#
# Usage:
#   module load PLINK/2.00a3.7-foss-2022a      # same plink used in the run
#   oneoff/check_kmer_orientation.sh --config master/config/config_RT.yaml
# =============================================================================

set -uo pipefail   # deliberately not -e: every check should report, not abort

CONFIG=""
NKMERS=20

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)  CONFIG="$2";  shift 2 ;;
        --nkmers)  NKMERS="$2";  shift 2 ;;
        -h|--help) sed -n '2,30p' "$0"; exit 1 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ -f "$CONFIG" ]] || { echo "--config <config.yaml> is required" >&2; exit 1; }

yaml_get() {
    sed -nE "s/^$1:[[:space:]]*\"?([^\"#]*[^\"# ])\"?[[:space:]]*(#.*)?$/\1/p" "$CONFIG" | head -1
}
RES="$(yaml_get results_dir)"
READS="$(yaml_get reads_dir)"
MP="$(yaml_get male_pattern)"
FP="$(yaml_get female_pattern)"

PLINKDIR="$RES/06_kmer/combined/plink"
INTER="$RES/06_kmer/combined/intermediate"
RESULTS="$RES/06_kmer/combined/results"
PHENO="$INTER/phenotype.pheno"
BFILE="$PLINKDIR/kmers_to_plink.0"
ASSOC="$PLINKDIR/kmers.assoc"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

hdr() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }
ok()  { printf '   [ok]   %s\n' "$1"; }
bad() { printf '   [FAIL] %s\n' "$1"; }
warn(){ printf '   [warn] %s\n' "$1"; }

echo "results_dir : $RES"
echo "male/female : /$MP/ , /$FP/"

# ---------------------------------------------------------------------------
hdr "1. phenotype.pheno -- did every sample get the label its name implies?"
# The rule builds this with ($2~mp?1:2), so ANY sample not matching the male
# pattern silently becomes a female, including a typo or an unmatched name.
# ---------------------------------------------------------------------------
if [[ ! -f "$PHENO" ]]; then
    bad "missing $PHENO"
else
    cat "$PHENO"
    awk -v mp="$MP" -v fp="$FP" 'NR>1{
        m = ($1 ~ mp); f = ($1 ~ fp)
        if (m && f)            { print "   [FAIL] " $1 " matches BOTH patterns"; e++ }
        else if (!m && !f)     { print "   [FAIL] " $1 " matches NEITHER pattern -> silently coded " $2; e++ }
        else if (m && $2 != 1) { print "   [FAIL] male "   $1 " coded " $2; e++ }
        else if (f && $2 != 2) { print "   [FAIL] female " $1 " coded " $2; e++ }
        n++
    } END { if(!e) print "   [ok]   all " n " samples labelled consistently with their names" }' "$PHENO"
fi

# ---------------------------------------------------------------------------
hdr "2. .fam -- what kmers_table_to_bed actually wrote as the phenotype"
# Column 6 must be literal 1/2 for PLINK to run a case/control test. If it
# was written as 1.000000/2.000000 PLINK treats it as a quantitative trait
# and --assoc emits completely different columns (see step 3).
# ---------------------------------------------------------------------------
if [[ ! -f "$BFILE.fam" ]]; then
    bad "missing $BFILE.fam"
else
    awk '{printf "   %-30s pheno=%s\n", $2, $6}' "$BFILE.fam"
    if awk '{if($6!="1" && $6!="2") exit 1}' "$BFILE.fam"; then
        ok "phenotype column is integer 1/2 -> PLINK will treat as case/control"
    else
        bad "phenotype column is NOT plain 1/2 -> PLINK treats it as QUANTITATIVE"
        bad "  --assoc then emits BETA/SE/R2/T/P and the \$5/\$6 awk reads the wrong columns"
    fi
fi

# ---------------------------------------------------------------------------
hdr "3. kmers.assoc -- case/control (10 cols) or quantitative (9 cols)?"
# ---------------------------------------------------------------------------
if [[ ! -f "$ASSOC" ]]; then
    bad "missing $ASSOC"
else
    head -1 "$ASSOC"
    NF_ASSOC=$(awk 'NR==1{print NF; exit}' "$ASSOC")
    case "$NF_ASSOC" in
        10) ok "10 columns: CHR SNP BP A1 F_A F_U A2 CHISQ P OR -- case/control, \$5/\$6 are F_A/F_U" ;;
        9)  bad "9 columns: quantitative output. \$5=BETA and \$6=SE, so the male/female"
            bad "  awk filters are meaningless and both result files are garbage." ;;
        *)  warn "unexpected column count: $NF_ASSOC" ;;
    esac
fi

# ---------------------------------------------------------------------------
hdr "4. PLINK log -- which samples did it call cases?"
# PLINK prints 'Among remaining phenotypes, N are cases and M are controls.'
# Cases are phenotype 2 = FEMALES under this pipeline's coding.
# ---------------------------------------------------------------------------
for L in "$PLINKDIR/kmers.log" "$RES/logs/06_kmer/kmer_association.log"; do
    [[ -f "$L" ]] && grep -iE 'case|control|quantitative|phenotype' "$L" | sed 's/^/   /'
done

# ---------------------------------------------------------------------------
hdr "5. bed batches -- is anything being silently ignored?"
# The rule only ever runs plink on kmers_to_plink.0. If kmers_table_to_bed
# split the table across batches, every k-mer past batch 0 was dropped.
# ---------------------------------------------------------------------------
mapfile -t BEDS < <(ls "$PLINKDIR"/kmers_to_plink.*.bed 2>/dev/null)
if [[ ${#BEDS[@]} -eq 0 ]]; then
    bad "no kmers_to_plink.*.bed found in $PLINKDIR"
elif [[ ${#BEDS[@]} -eq 1 ]]; then
    ok "single batch -- nothing dropped"
else
    bad "${#BEDS[@]} bed batches exist but only kmers_to_plink.0 is analysed:"
    printf '          %s\n' "${BEDS[@]}"
fi

# ---------------------------------------------------------------------------
hdr "6. THE TEST -- who actually carries the 'male-specific' k-mers?"
# Pull k-mers out of male_specific_kmers.assoc, re-extract their genotypes
# with --recode A, and cross-tabulate dosage against sex. This bypasses
# F_A/F_U entirely and reads the genotype matrix directly.
#
#   A1A1 = homozygous for A1, A2A2 = homozygous for the other allele.
#   Presence/absence k-mers are coded as opposing homozygotes, so every
#   sample should sit at one extreme or the other.
# ---------------------------------------------------------------------------
MALE_ASSOC="$RESULTS/male_specific_kmers.assoc"
FEMALE_ASSOC="$RESULTS/female_specific_kmers.assoc"
echo "   male_specific_kmers.assoc  : $( [[ -f $MALE_ASSOC   ]] && wc -l < "$MALE_ASSOC"   || echo MISSING) k-mers"
echo "   female_specific_kmers.assoc: $( [[ -f $FEMALE_ASSOC ]] && wc -l < "$FEMALE_ASSOC" || echo MISSING) k-mers"

if ! command -v plink >/dev/null; then
    warn "plink not on PATH -- skipping; module load the same PLINK used in the run"
elif [[ ! -s "$MALE_ASSOC" ]]; then
    warn "$MALE_ASSOC is empty -- nothing to test"
else
    cut -f2 "$MALE_ASSOC" | head -"$NKMERS" > "$WORK/kmers.txt"
    echo "   testing $(wc -l < "$WORK/kmers.txt") k-mers from male_specific_kmers.assoc"
    plink --noweb --bfile "$BFILE" --extract "$WORK/kmers.txt" \
          --recode A --out "$WORK/check" >/dev/null 2>&1

    if [[ ! -f "$WORK/check.raw" ]]; then
        bad "plink --recode A produced no output; run it by hand to see the error"
    else
        echo
        awk -v mp="$MP" -v fp="$FP" '
        NR==1 { n = NF - 6; next }
        {
            a1=0; a2=0; het=0; na=0
            for (i=7; i<=NF; i++) {
                if      ($i == 2)    a1++
                else if ($i == 0)    a2++
                else if ($i == 1)    het++
                else                 na++
            }
            sex = ($2 ~ mp) ? "M" : ($2 ~ fp) ? "F" : "?"
            printf "   %-28s %s  A1A1=%-4d A2A2=%-4d het=%-4d NA=%-4d\n", $2, sex, a1, a2, het, na
            if (sex == "M") { m_a1 += a1; m_tot += n }
            if (sex == "F") { f_a1 += a1; f_tot += n }
        }
        END {
            printf "\n   males  carrying A1: %d/%d\n", m_a1, m_tot
            printf "   females carrying A1: %d/%d\n", f_a1, f_tot
            print  ""
            if (m_a1 == m_tot && f_a1 == 0) {
                print "   -> A1 tracks MALES (fixed in males, absent in females)."
                print "      If A1 is the PRESENT allele, the file is correctly named."
            } else if (f_a1 == f_tot && m_a1 == 0) {
                print "   -> A1 tracks FEMALES (fixed in females, absent in males)."
                print "      If A1 is the PRESENT allele, the two result files are SWAPPED."
            } else {
                print "   -> NOT a clean split. These k-mers do not separate the sexes at all;"
                print "      the awk filter is not selecting what the rule claims."
            }
        }' "$WORK/check.raw"

        echo
        echo "   A1/A2 as stored in the .bim for the first tested k-mer:"
        # .bim is always 6 columns: CHR SNP CM BP A1 A2
        awk -v k="$(head -1 "$WORK/kmers.txt")" '
            $2 == k { if (NF == 6) printf "     A1=%s  A2=%s\n", $5, $6
                      else         printf "     [warn] unexpected .bim format (%d cols)\n", NF
                      found = 1; exit }
            END { if (!found) print "     [warn] k-mer not found in .bim" }' "$BFILE.bim"
        echo "   (A1 fixed in a sex only tells you that allele tracks that sex --"
        echo "    step 7 settles whether A1 means the k-mer is PRESENT or ABSENT.)"
    fi
fi

# ---------------------------------------------------------------------------
hdr "7. GROUND TRUTH -- grep the k-mer in the raw reads"
# Independent of PLINK, kmersGWAS, and every allele convention: if the k-mer
# is genuinely male-specific it is in the male FASTQs and not the female
# ones. Check both strands, since reads come from either.
# ---------------------------------------------------------------------------
if [[ -s "$MALE_ASSOC" ]]; then
    K=$(cut -f2 "$MALE_ASSOC" | head -1)
    KRC=$(echo "$K" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    M_SAMPLE=$(ls "$READS" | grep -m1 "$MP")
    F_SAMPLE=$(ls "$READS" | grep -m1 "$FP")
    cat <<EOF
   k-mer   : $K
   revcomp : $KRC

   Run these two (a few minutes each; head keeps it bounded):

     for S in $M_SAMPLE $F_SAMPLE; do
       echo -n "\$S: "
       zcat $READS/\$S/*_R1*.gz 2>/dev/null | head -n 40000000 \\
         | grep -cE "$K|$KRC"
     done

   Expect a nonzero count in $M_SAMPLE and zero in $F_SAMPLE.
   Reversed => the two result files are swapped.
   Zero in both => raise the head limit; coverage may be too low to catch it.
EOF
fi

echo
echo "Done."
