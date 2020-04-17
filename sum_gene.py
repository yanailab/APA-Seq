#!/usr/local/bin/python

import sys,os
import numpy as np
#import csv

filename = '../depth_oriented_test/C14B9.4.txt'
gene_np = np.genfromtxt(filename,delimiter="\t")
print gene_np.sum()

#gene_name = sys.argv[1]

