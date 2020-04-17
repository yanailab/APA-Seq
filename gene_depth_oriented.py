#!/usr/bin/python

import os
import sys
import numpy
import glob

samplesDir = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/samples/'
sampleDirs = os.listdir(samplesDir)
geneDir = '/illumina1/YanaiLab/eitan/UTRome/genes/'
genes = glob.glob(geneDir + "*.fasta")
depthDir = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/depth_oriented/'
gene_orientation = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/scripts/genes_orientation.txt'

with open(gene_orientation) as ori:
        content = ori.readlines()
orientation_dict = {}
for line in content:
        (gene, orientation) = line.split(' ')
        orientation_dict[gene] = orientation
        
#for i in range(20):
for i in range(len(genes)):
	fastaFile = genes[i] #get fasta of gene #i
	gName = fastaFile.replace(geneDir,'',1).replace('.fasta','',1)
	if i % 100 == 0:
		print(i)
	with open(fastaFile) as f:
		content = f.readlines()
        if '-' in orientation_dict[gName]:
             strand = 'minus'
             print '-'
        elif '+' in orientation_dict[gName]:
             strand = 'plus'
             print '+'
        else:
             strand = 'notfound'  
             print(gName + ' strand not known')
	
	# now we declare a vector of zeros the same length as the length of the region in the genome we are trying to map to:
	vector_of_locations=content[0].split(":")[1].split("-")
	length_of_vector_of_zeros=int(max(vector_of_locations))-int(min(vector_of_locations))
	
	# generating a position * sample two dimentional array
	matrix = [ [ 0 for m in range(len(sampleDirs)) ] for x in range(length_of_vector_of_zeros)]
        Matrix = numpy.matrix(matrix)
	
	sampleDirs = sorted(sampleDirs)
	for x,sampleDir in enumerate(sampleDirs):
		sampleDir = sampleDirs[x]
    		samFile = samplesDir + sampleDir + "/plus/" + gName + '.sam'
		if os.path.isfile(samFile):
        		with open(samFile) as sam:
    				content = sam.readlines()
    			for line in content:
           			temporal=line.split("\t")
           			if ( temporal[2]!="*" ) and ( len(temporal)>8 ) and ( strand!='notfound' ): # the read did map to some area of the genome, its a valid read
                                    align_location = int(temporal[3])
                                    seq_length = int(len(temporal[9]))
                                    if strand == 'minus':
                                       ind_in_vect = align_location
                                    elif strand == 'plus':
                                       ind_in_vect = align_location + seq_length

                                    if (length_of_vector_of_zeros < ind_in_vect):
                                        print("index out of bound at " + str(ind_in_vect) + " at gene:" + gName + " at sample:" + sampleDir)
                                    else:
                                        Matrix[ind_in_vect-1,x] += 1
                                                
	# write matrix to file
	outfile = depthDir + gName + ".txt"
	numpy.savetxt(outfile, Matrix, fmt='%d', delimiter='\t', newline='\n')
      #  print 'sum of matrix for ' + gName + ' is ' + Matrix.sum().sum()
#SAM:
# 0      0       CHROMOSOME_III:8917853-8923350  4471    0       66M     *       0       0       CAAAAACTTTTTTTTTTTTGTACAAAAACGTGTCCCAGACGTGTAAATGTTATAAGACACACACGT      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      AS:i:-36        XN:i:0  XM:i:18 XO:i:0  XG:i:0  NM:i:18 MD:Z:0A3T5G3A16T3A3G0A0G0A3T0C0A6T2T2A1G0C1     YT:Z:UU 
