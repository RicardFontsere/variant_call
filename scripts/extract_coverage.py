#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("coverage", type=str, nargs="+",
                    help="A folder of individual sample files containing normalized coverage data")
parser.add_argument("outfile", type=str,
                    help="Output file with average coverage across all samples")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder)]

def extract_depth(source, depthdict):
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if line:  # Skip empty lines
                parts = line.split('\t')  # Explicitly split by tab
                if len(parts) == 5:
                    chromosome = parts[0]
                    window_start = parts[1]
                    window_end = parts[2]
                    depth = float(parts[4])
                    
                    # Create a tuple key instead of string concatenation
                    # This completely avoids any delimiter issues
                    chromosome_block = (chromosome, window_start, window_end)
                    depthdict[chromosome_block].append(depth)
    return depthdict

def average_depth(coverage_depth):
    averagedict = {}
    for chromosome_block in coverage_depth:
        coverage_list = coverage_depth[chromosome_block]
        sum_coverage_list = sum(coverage_list)
        average = sum_coverage_list / len(coverage_list)
        coverage_values = ",".join(map(str, coverage_list))
        logaverage = np.log2(average + 1)
        
        # Store with the same tuple key
        averagedict[chromosome_block] = [coverage_values, sum_coverage_list, average, logaverage]
    return averagedict

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    # Get files
    if len(args.coverage) == 1:
        if os.path.isdir(args.coverage[0]):
            infiles = list_folder(args.coverage[0])
        else:
            infiles = args.coverage
    else:
        infiles = args.coverage

    # Extract coverage
    depthdict = defaultdict(list)
    for infile in infiles:
        print("Processing:", infile)
        depthdict = extract_depth(infile, depthdict)
    print("Number of chromosome blocks with coverage =", len(depthdict))

    # Calculate average coverage
    average = average_depth(depthdict)
    print("Number of chromosome blocks with average coverage =", len(average))

    # Write output
    with open(args.outfile, "w") as outfile:
        print("Writing output to:", args.outfile)
        header = "Chromosome,WindowStart,WindowEnd,CovList,Sum,Average,Logaverage\n"
        outfile.write(header)
        
        # Sort by chromosome and position for organized output
        for chromosome_block in sorted(average.keys()):
            # Unpack the tuple
            chromosome, window_start, window_end = chromosome_block
            
            # Get the statistics
            stats = average[chromosome_block]
            
            # Since CovList now contains commas, we need to quote it for proper CSV format
            cov_list = f'"{stats[0]}"'  # Quote the coverage list
            
            # Write the row
            outfile.write(f"{chromosome},{window_start},{window_end},{cov_list},{stats[1]},{stats[2]},{stats[3]}\n")

if __name__ == '__main__':
    main()