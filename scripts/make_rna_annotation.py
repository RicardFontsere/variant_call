#!/usr/bin/env python3
# =============================================================================
# Build the transcript annotation table edgeR attaches to the counts, from an
# NCBI RefSeq GFF and the transcriptome cut from it (rule rna_transcriptome).
#
# One row per transcript of the transcriptome that has an NCBI accession
# (transcript_id attribute: XM_, XR_, NM_, NR_):
#   gid       transcript accession, as in the transcriptome after the rna-
#             prefix is stripped
#   gname     gene symbol (gene attribute)
#   gproduct  product description (product attribute)
#   chr       Chr<n> from the chromosome attribute of the region lines,
#             ChrUnknown for unplaced scaffolds
#   chr_id    sequence accession (GFF column 1)
#   start/end genomic span of the transcript
#   length    spliced transcript length, from the transcriptome
#
# A transcript annotated at a second locus gets the ID rna-<accession>-2: it is
# left out, so every accession has one row and one location.
#
# Usage: make_rna_annotation.py genomic.gff transcriptome.fa out.tsv
# =============================================================================
import sys
from urllib.parse import unquote

if len(sys.argv) != 4:
    sys.exit("usage: make_rna_annotation.py genomic.gff transcriptome.fa out.tsv")
gff_file, fasta_file, out_file = sys.argv[1:]


def attributes(column9):
    attrs = {}
    for field in column9.rstrip(";").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


# Transcript lengths, keyed by the (prefix-stripped) fasta name
lengths = {}
name = None
with open(fasta_file) as fasta:
    for line in fasta:
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name is not None:
            lengths[name] += len(line.strip())

chromosome = {}
transcripts = {}
with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            continue
        seqid, feature, start, end, attrs = cols[0], cols[2], cols[3], cols[4], attributes(cols[8])
        if feature == "region":
            chrom = attrs.get("chromosome", "")
            if chrom and chrom != "Unknown":
                chromosome[seqid] = "Chr" + chrom
            elif attrs.get("genome") == "mitochondrion":
                chromosome[seqid] = "ChrMT"
            continue
        tid = attrs.get("transcript_id")
        # Only the primary placement: a second locus has the ID rna-<tid>-2
        if not tid or attrs.get("ID") != "rna-" + tid or tid in transcripts:
            continue
        transcripts[tid] = (attrs.get("gene", "NA"), attrs.get("product", "NA"),
                            seqid, start, end)

n_written = 0
with open(out_file, "w") as out:
    out.write("gid\tgname\tgproduct\tchr\tchr_id\tstart\tend\tlength\n")
    for tid, (gname, product, seqid, start, end) in transcripts.items():
        if tid not in lengths:
            continue  # not extracted by gffread, so never quantified
        out.write("\t".join([tid, gname, product, chromosome.get(seqid, "ChrUnknown"),
                             seqid, start, end, str(lengths[tid])]) + "\n")
        n_written += 1

print(f"Transcripts in the transcriptome: {len(lengths)}")
print(f"Transcripts with an accession in the GFF: {len(transcripts)}")
print(f"Annotated (in both): {n_written}")
if n_written == 0:
    sys.exit("No transcript of the GFF is in the transcriptome: do the names match?")
