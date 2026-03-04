import csv
import sys

fileout=(sys.argv[2])
fout= open(fileout, 'w')

with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        kmer = strLine[1]
        pval=strLine[8]
        fout.write(">"+pval+"\n"+kmer+"\n")
fout.close()