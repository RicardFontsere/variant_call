#!/usr/bin/env python3
"""
Compute per-sample genotype depth (DP) percentiles from bcftools query output.

Reads a tab-delimited file where:
  - First line: sample names (header)
  - Remaining lines: per-sample DP values (one column per sample, one row per site)
  - Missing values are represented as "."

Outputs:
  1. Per-sample percentiles table (TSV) - for notebook diagnostics
  2. Per-sample DP histograms (TSV) - for notebook plotting
  3. Global thresholds file (TSV) - for pipeline consumption
     Global lower = median of per-sample 5th percentiles
     Global upper = median of per-sample 99th percentiles

Usage:
    python compute_dp_thresholds.py <dp_table> <percentiles_out> <histograms_out> <thresholds_out>
"""

import sys
import os
import numpy as np


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <dp_table> <percentiles_out> <histograms_out> <thresholds_out>")
        sys.exit(1)

    dp_table_path = sys.argv[1]
    percentiles_out = sys.argv[2]
    histograms_out = sys.argv[3]
    thresholds_out = sys.argv[4]

    # Read header (sample names)
    with open(dp_table_path, 'r') as f:
        header = f.readline().rstrip('\n').split('\t')

    samples = [s for s in header if s]  # remove empty strings
    n_samples = len(samples)
    print(f"Found {n_samples} samples")

    # Read DP values into per-sample lists
    dp_values = {s: [] for s in samples}

    with open(dp_table_path, 'r') as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            for i, sample in enumerate(samples):
                val = fields[i] if i < len(fields) else "."
                if val != "." and val != "":
                    try:
                        dp_values[sample].append(int(val))
                    except ValueError:
                        pass

    # Compute per-sample statistics
    percentile_data = {}
    p5_list = []
    p99_list = []

    for sample in samples:
        values = np.array(dp_values[sample])
        if len(values) == 0:
            print(f"Warning: no DP values for sample {sample}")
            continue

        p5 = np.percentile(values, 5)
        p10 = np.percentile(values, 10)
        p25 = np.percentile(values, 25)
        p50 = np.percentile(values, 50)
        p75 = np.percentile(values, 75)
        p90 = np.percentile(values, 90)
        p95 = np.percentile(values, 95)
        p99 = np.percentile(values, 99)
        mean_dp = np.mean(values)
        n_sites = len(values)

        percentile_data[sample] = {
            'n_sites': n_sites, 'mean': mean_dp,
            'p5': p5, 'p10': p10, 'p25': p25, 'p50': p50,
            'p75': p75, 'p90': p90, 'p95': p95, 'p99': p99
        }

        p5_list.append(p5)
        p99_list.append(p99)
        print(f"  {sample}: n={n_sites}, mean={mean_dp:.1f}, p5={p5:.0f}, p99={p99:.0f}")

    # Write per-sample percentiles table (for notebook)
    with open(percentiles_out, 'w') as f:
        f.write("sample\tn_sites\tmean\tp5\tp10\tp25\tp50\tp75\tp90\tp95\tp99\n")
        for sample in samples:
            if sample in percentile_data:
                d = percentile_data[sample]
                f.write(f"{sample}\t{d['n_sites']}\t{d['mean']:.2f}\t"
                        f"{d['p5']:.1f}\t{d['p10']:.1f}\t{d['p25']:.1f}\t{d['p50']:.1f}\t"
                        f"{d['p75']:.1f}\t{d['p90']:.1f}\t{d['p95']:.1f}\t{d['p99']:.1f}\n")

    # Write per-sample histograms (for notebook plotting)
    # Binned from 0 to max_dp in integer bins
    max_dp = max(max(dp_values[s]) for s in samples if len(dp_values[s]) > 0)
    max_bin = min(max_dp, 500)  # cap at 500 for sanity

    with open(histograms_out, 'w') as f:
        f.write("sample\tdepth\tcount\n")
        for sample in samples:
            if len(dp_values[sample]) == 0:
                continue
            values = np.array(dp_values[sample])
            counts, _ = np.histogram(values, bins=range(0, max_bin + 2))
            for depth, count in enumerate(counts):
                if count > 0:
                    f.write(f"{sample}\t{depth}\t{count}\n")

    # Compute global thresholds: median of per-sample percentiles, rounded to integer
    global_lower = int(np.floor(np.median(p5_list)))
    global_upper = int(np.ceil(np.median(p99_list)))

    print(f"\nGlobal thresholds: DP_lower={global_lower}, DP_upper={global_upper}")

    with open(thresholds_out, 'w') as f:
        f.write(f"DP_lower\t{global_lower}\n")
        f.write(f"DP_upper\t{global_upper}\n")

    print("Done.")


if __name__ == "__main__":
    main()
