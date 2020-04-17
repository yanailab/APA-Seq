
from __future__ import print_function
import time
import sys
import re
import os
import glob
import numpy as np

def process_align_summary(infile):

    num="\d+[\.|\d|e|\-]*"
    print(infile)
    with open (infile,'r') as a:
	cline = a.readline()
	while cline.startswith("Warning"):
	    cline = a.readline()	
        cline = a.readline()
	cline = a.readline()
        regex = r"^\s+\d+\s+\(("+num+")\%\)\s+aligned 0 times$"
        m = re.search(regex,cline)
        if m:
            unaligned = m.groups()
        cline= a.readline()
	regex = r"^\s+\d+\s+\(("+num+")\%\)\s+aligned exactly 1 time$"
        m = re.search(regex,cline)
        if m:
            unique = m.groups()
        cline= a.readline()
        regex = r"^\s+\d+\s+\(("+num+")\%\)\s+aligned \>1 times$"
        m = re.search(regex,cline)
        if m:
            multi_aligned = m.groups()
        array = [str(unaligned[0]), str(unique[0]),str(multi_aligned[0])]
    return array


def get_stats(od, outfile):

    outfile = os.path.join(os.getcwd(),outfile)
    os.chdir(od)

    top = open (outfile,'w')
    header = "\t".join([ ' % unaligned', ' % aligned once', '% multi aligned']) + "\n"
    top.write(header)
    
    arr = []
    for file in glob.glob("*.stat"):
        topl = process_align_summary(file)
        tmp='\t'.join(topl)+'\n'
        top.write(tmp)
	arr.append(topl)

    top.close()
    
    l = np.asarray(arr).astype(np.float)
    means = l.mean(axis=0)
    print(means)

if __name__ == "__main__":
    outdir=sys.argv[1] # aligned/"
    outfile=sys.argv[2] # report.txt"
    get_stats(outdir,outfile)


