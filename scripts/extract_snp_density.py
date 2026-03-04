"This script reads snp density information from multiple samples and calculates and outputs average snp density for each chromosome block"
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
parser.add_argument("infolder", type=str, nargs="+",
                    help="A folder with individual snp density files")
parser.add_argument("outfile", type=str,
                    help="Output file with average SNP density across all samples")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".txt")]

def extract_snpdensity(source, snpdict):
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if line:  # Skip empty lines
                parts = line.split()  # Split by whitespace
                
                # Work with 7-column format
                if len(parts) == 7:
                    scaffold = parts[0]  # Chromosome name
                    window_start = parts[1]  # Start position
                    window_end = parts[2]  # End position
                    snpdensity = float(parts[6])  # SNP density value (seventh column)
                    
                    # Use tuple as key to avoid delimiter issues
                    scaffold_key = (scaffold, window_start, window_end)
                    snpdict[scaffold_key].append(snpdensity)
    
    return snpdict

def average_snpdensity(snpdensity):
    averagedict = {}
    for chromosome_block in snpdensity:
        snplist = snpdensity[chromosome_block]
        sum_snplist = sum(snplist)
        average = sum_snplist / len(snplist)
        # Use comma separator for consistency
        snp_values = ",".join(map(str, snplist))
        logaverage = np.log2(average + 1)
        
        # Store with the same tuple key
        averagedict[chromosome_block] = [snp_values, sum_snplist, average, logaverage]
    
    return averagedict

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    # Get files
    if len(args.infolder) == 1:
        if os.path.isdir(args.infolder[0]):
            infiles = list_folder(args.infolder[0])
        else:
            infiles = args.infolder
    else:
        infiles = args.infolder

    # Extract snp density
    snpdensitydict = defaultdict(list)
    for infile in infiles:
        print("Processing:", infile) 
        snpdensitydict = extract_snpdensity(infile, snpdensitydict)
    print("Number of chromosome blocks with snp density =", len(snpdensitydict))

    # Calculate average snp density
    average = average_snpdensity(snpdensitydict)
    print("Number of chromosome blocks with average snp density =", len(average))

    # Write output
    with open(args.outfile, "w") as outfile:
        print("Writing output to:", args.outfile)
        header = "Chromosome,WindowStart,WindowEnd,SnpDensityList,Sum,Average,Logaverage\n"
        outfile.write(header)
        
        # Sort by chromosome and position for organized output
        for chromosome_block in sorted(average.keys()):
            # Unpack the tuple
            chromosome, window_start, window_end = chromosome_block
            
            # Get the statistics
            stats = average[chromosome_block]
            
            # Since SnpDensityList now contains commas, quote it for proper CSV
            snp_list = f'"{stats[0]}"'
            
            # Write the row
            outfile.write(f"{chromosome},{window_start},{window_end},{snp_list},{stats[1]},{stats[2]},{stats[3]}\n")

if __name__ == '__main__':
    main()