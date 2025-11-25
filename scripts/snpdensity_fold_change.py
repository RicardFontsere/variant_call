"This script reads in female and male snp density data and calculates and outputs M:F snp density fold change for each chromosome block"
#==============================================================================
import argparse
import sys
import os
import csv
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str,
                    help="A file with average snp density data from females")
parser.add_argument("males", type=str,
                    help="A file with average snp density data from males")
parser.add_argument("outfile", type=str,
                    help="Output file containing fold change values")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def extract_snps(source):
    snpsdict = {}
    with open(source, "r") as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 7:  # Ensure we have all fields
                chromosome = row[0]
                windowStart = row[1]
                windowEnd = row[2]
                # Use tuple as key to avoid delimiter issues
                chromosome_block = (chromosome, windowStart, windowEnd)
                
                # Extract the statistics (skip SnpDensityList at index 3)
                sum_snpdensity = row[4]
                average_snpdensity = row[5]
                logaverage_snpdensity = row[6]
                
                snpsdict[chromosome_block] = [sum_snpdensity, average_snpdensity, logaverage_snpdensity]
    return snpsdict

def combine_density(females_snps, males_snps):
    combined_dict = {}
    for block in females_snps:
        if block in males_snps:  # Only process blocks that exist in both
            fem_sum = float(females_snps[block][0])
            fem_average = float(females_snps[block][1])
            fem_logaverage = float(females_snps[block][2])
            
            mal_sum = float(males_snps[block][0])
            mal_average = float(males_snps[block][1])
            mal_logaverage = float(males_snps[block][2])
            
            # Calculate M:F log fold change
            combined_logaverage = mal_logaverage - fem_logaverage
            
            combined_dict[block] = [
                combined_logaverage,
                mal_average,
                mal_logaverage,
                fem_average,
                fem_logaverage
            ]
    
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    # Extract SNP density in males and females
    females_snps = extract_snps(args.females)
    print("Number of female chromosome blocks =", len(females_snps))
    males_snps = extract_snps(args.males)
    print("Number of male chromosome blocks =", len(males_snps))

    # Check for blocks that don't match
    females_only = set(females_snps.keys()) - set(males_snps.keys())
    males_only = set(males_snps.keys()) - set(females_snps.keys())
    
    if females_only:
        print(f"Warning: {len(females_only)} blocks found only in females file")
    if males_only:
        print(f"Warning: {len(males_only)} blocks found only in males file")

    # Combine density
    combined_density = combine_density(females_snps, males_snps)
    print("Number of filtered chromosome blocks =", len(combined_density))

    # Write output
    with open(args.outfile, "w") as outfile:
        header = "Chromosome,WindowStart,WindowEnd,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage\n"
        outfile.write(header)
        
        # Sort by chromosome and position for organized output
        for block in sorted(combined_density.keys()):
            # Unpack the tuple
            chromosome, windowStart, windowEnd = block
            
            # Write chromosome info
            outfile.write(f"{chromosome},{windowStart},{windowEnd},")
            
            # Write the statistics
            stats = combined_density[block]
            outfile.write(",".join(map(str, stats)))
            outfile.write("\n")

if __name__ == '__main__':
    main()