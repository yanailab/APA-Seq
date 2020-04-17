#!/usr/bin/python

import os
import sys
import numpy
import glob

samplesDir = '/illumina1/YanaiLab/eitan/UTRome/samples/'
sampleDirs = os.listdir(samplesDir)
geneDir = '/illumina1/YanaiLab/eitan/UTRome/genes/'
genes = glob.glob(geneDir + "*.fasta")
depthDir = '/illumina1/YanaiLab/eitan/UTRome/depth/'


for i in range(len(genes)):
	fastaFile = genes[i] #get fasta of gene #i
	gName = fastaFile.replace(geneDir,'',1).replace('.fasta','',1)
	if i % 100 == 0:
		print(i)
	with open(fastaFile) as f:
		content = f.readlines()
	
	# now we declare a vector of zeros the same length as the length of the region in the genome we are trying to map to:
	vector_of_locations=content[0].split(":")[1].split("-")
	length_of_vector_of_zeros=int(max(vector_of_locations))-int(min(vector_of_locations))
	
	# generating a position * sample two dimentional array
	Matrix = [ [ 0 for m in range(len(sampleDirs)) ] for x in range(length_of_vector_of_zeros)]
	
	sampleDirs = sorted(sampleDirs)
	for x,sampleDir in enumerate(sampleDirs):
		sampleDir = sampleDirs[x]
    		samFile = samplesDir + sampleDir + "/plus/" + gName + '.sam'
		if os.path.isfile(samFile):
        		with open(samFile) as sam:
    				content = sam.readlines()    			
    			for line in content:
           			temporal=line.split("\t")
           			if temporal[2]!="*" and len(temporal)>8: # the read did map to some area of the genome, its a valid read
						ind_in_vect=int(temporal[3]) # alignment start position
						if length_of_vector_of_zeros >= ind_in_vect:
							Matrix[ind_in_vect-1][x] = Matrix[ind_in_vect-1][x] + 1
						else:
							print("index out of bound at " + str(ind_in_vect) + " at gene:" + gName + " at sample:" + sampleDir)

	# write matrix to file
	outfile = depthDir + gName + ".txt"
	numpy.savetxt(outfile, Matrix, fmt='%d', delimiter='\t', newline='\n')
                 


     
