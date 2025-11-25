"This script reads in female and male coverage data and calculates and outputs M:F coverage fold change for each chromosome block"
#==============================================================================
import argparse
import sys
import csv
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("females", type=str,
                    help="A file with average coverage data from females")
parser.add_argument("males", type=str,
                    help="A file with average coverage data from males")
parser.add_argument("outfile", type=str,
                    help="An output file with M:F coverage fold change values")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def extract_depth(source):
    depthdict = {}
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
                
                # Extract the statistics (skip CovList at index 3)
                sumdepth = row[4]
                averagedepth = row[5]
                logaveragedepth = row[6]
                
                depthdict[chromosome_block] = [sumdepth, averagedepth, logaveragedepth]
    return depthdict

def combine_depth(females_depth, males_depth):
    combined_dict = {}
    for block in females_depth:
        if block in males_depth:  # Only process blocks that exist in both
            fem_sumdepth = float(females_depth[block][0])
            fem_averagedepth = float(females_depth[block][1])
            fem_logaveragedepth = float(females_depth[block][2])
            
            mal_sumdepth = float(males_depth[block][0])
            mal_averagedepth = float(males_depth[block][1])
            mal_logaveragedepth = float(males_depth[block][2])

            # Calculate M:F log fold change
            combined_logaveragedepth = mal_logaveragedepth - fem_logaveragedepth

            combined_dict[block] = [
                combined_logaveragedepth,
                mal_averagedepth,
                mal_logaveragedepth,
                fem_averagedepth,
                fem_logaveragedepth
            ]
    
    return combined_dict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    # Extract coverage in males and females
    females_depth = extract_depth(args.females)
    print("Number of chromosome blocks with coverage in females =", len(females_depth))
    males_depth = extract_depth(args.males)
    print("Number of chromosome blocks with coverage in males =", len(males_depth))

    # Check for blocks that don't match
    females_only = set(females_depth.keys()) - set(males_depth.keys())
    males_only = set(males_depth.keys()) - set(females_depth.keys())
    
    if females_only:
        print(f"Warning: {len(females_only)} blocks found only in females file")
    if males_only:
        print(f"Warning: {len(males_only)} blocks found only in males file")

    # Combine coverage
    combined_depth = combine_depth(females_depth, males_depth)
    print("Number of chromosome blocks with combined coverage =", len(combined_depth))

    # Write output
    with open(args.outfile, "w") as outfile:
        header = "Chromosome,WindowStart,WindowEnd,MFLogaverage,Maverage,Mlogaverage,Faverage,Flogaverage\n"
        outfile.write(header)
        
        # Sort by chromosome and position for organized output
        for block in sorted(combined_depth.keys()):
            # Unpack the tuple
            chromosome, windowStart, windowEnd = block
            
            # Write chromosome info
            outfile.write(f"{chromosome},{windowStart},{windowEnd},")
            
            # Write the statistics
            stats = combined_depth[block]
            outfile.write(",".join(map(str, stats)))
            outfile.write("\n")

if __name__ == '__main__':
    main()