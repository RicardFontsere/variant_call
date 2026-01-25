#!/usr/bin/env python3
"""
Generate Picard-style interval files from a reference genome .fai index and .dict file.

Chromosomes (defined by prefix in config) get individual interval files.
Unplaced contigs (defined by prefix in config) are grouped into one file.

Usage:
    python generate_intervals.py <reference.fai> <reference.dict> <output_dir> <chrom_prefix> <contig_prefix>

Example:
    python generate_intervals.py RT_genome.fna.fai RT_genome.dict intervals NC_ NW_
"""

import sys
import os


def read_dict_header(dict_path):
    """Read the sequence dictionary and return header lines."""
    header_lines = []
    with open(dict_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                header_lines.append(line.rstrip('\n'))
    return header_lines


def write_picard_interval_file(filepath, header_lines, intervals):
    """
    Write a Picard-style interval list file.
    
    Format:
    @HD	VN:1.6
    @SQ	SN:chr1	LN:12345678
    ...
    chr1	1	12345678	+	.
    """
    with open(filepath, 'w') as f:
        # Write header
        for line in header_lines:
            f.write(line + '\n')
        
        # Write intervals (tab-separated: contig, start, end, strand, name)
        for name, length in intervals:
            f.write(f"{name}\t1\t{length}\t+\t.\n")


def main():
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <reference.fai> <reference.dict> <output_dir> <chrom_prefix> <contig_prefix>")
        print(f"Example: {sys.argv[0]} RT_genome.fna.fai RT_genome.dict intervals NC_ NW_")
        sys.exit(1)

    fai_path = sys.argv[1]
    dict_path = sys.argv[2]
    output_dir = sys.argv[3]
    chrom_prefix = sys.argv[4]
    contig_prefix = sys.argv[5]

    # Validate input files exist
    if not os.path.exists(fai_path):
        print(f"Error: FAI file not found: {fai_path}")
        sys.exit(1)
    if not os.path.exists(dict_path):
        print(f"Error: Dict file not found: {dict_path}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Read the sequence dictionary header
    header_lines = read_dict_header(dict_path)
    if not header_lines:
        print(f"Error: No header lines found in {dict_path}")
        sys.exit(1)

    chromosomes = []
    contigs = []

    # Parse FAI file to get sequence names and lengths
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            length = parts[1]

            if name.startswith(chrom_prefix):
                chromosomes.append((name, length))
            elif name.startswith(contig_prefix):
                contigs.append((name, length))

    # Write individual chromosome interval files
    for name, length in chromosomes:
        filepath = os.path.join(output_dir, f"{name}.interval_list")
        write_picard_interval_file(filepath, header_lines, [(name, length)])

    # Write grouped contigs interval file
    if contigs:
        filepath = os.path.join(output_dir, "unplaced_contigs.interval_list")
        write_picard_interval_file(filepath, header_lines, contigs)

    print(f"Created {len(chromosomes)} chromosome interval files")
    print(f"Created 1 unplaced contigs file with {len(contigs)} contigs")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()