#!/usr/bin/python
import re
import glob
import os
import time
import subprocess
    

baseDir = '/illumina1/YanaiLab/eitan/UTRome/samples/'
geneDir = '/illumina1/YanaiLab/eitan/UTRome/genes/'
reference_fasta="c_elegans.WS230_spikein.genomic.fa"

# Generating Fasta and Bowtie-2 indexes for all genes! Should only run once for this species!
if os.path.isdir(geneDir):
	print("gene directory exists: " + geneDir)
else: 
	PLUS = open('new_g_p_sorted.txt').readlines()
	MINUS = open('new_g_m_sorted.txt').readlines()
	
	CHROMOSOMES = ['CHROMOSOME_I', 'CHROMOSOME_II', 'CHROMOSOME_III' , 'CHROMOSOME_IV', 'CHROMOSOME_V' , 'CHROMOSOME_X']
	
	minAddresses = dict()
	plusAddresses = dict()
	
	for chrom in CHROMOSOMES:
		print(chrom)
		minusRecords = filter(lambda x:re.search(chrom + '\t', x), MINUS)
		plusRecords = filter(lambda x:re.search(chrom + '\t', x), PLUS)
		minAddressesArr = []
		plusAddressesArr = []
		for minRec in minusRecords:
			tokens = minRec.split('\t')
			xf = tokens[-1].split(':')[-1].replace('\n','')
			minAddressesArr.append([xf, tokens[3], tokens[4]]) # gene name, start, end
		for plusRec in plusRecords:
			tokens = plusRec.split('\t')
			xf = tokens[-1].split(':')[-1].replace('\n','')
			plusAddressesArr.append([xf, tokens[3], tokens[4]])	# gene name, start, end	
		
		minAddresses[minAddressesArr[0][0]] = [chrom,0,minAddressesArr[0][-1]]
		
		for i,minAdd in enumerate(minAddressesArr):
			if (i == 0):
				continue
			gene_name = minAdd[0]
			start_position = max(int(minAdd[1])-5000, int(minAddressesArr[i-1][-1]))
			end_position = minAdd[-1]
			fasta_path = geneDir + gene_name + '.fasta'
			#print ("gene name " + gene_name + " start:" + str(start_position) + " end:" + str(end_position) + " fasta_path:" + fasta_path + " reference_fasta:" + reference_fasta)
			command_line = "samtools faidx " + reference_fasta + " " + chrom + ":" + str(start_position) + "-" + str(end_position) + " > " + fasta_path
			os.system(command_line)
			minAddresses[minAdd[0]] = [ chrom, start_position, end_position ]
			
		for i,plusAdd in enumerate(plusAddressesArr[:-1]):
			gene_name = plusAdd[0]
			start_position = plusAdd[1]
			end_position = min(int(plusAdd[-1])+5000, int(plusAddressesArr[i+1][1]))
			fasta_path = geneDir + gene_name + '.fasta'
			command_line = "samtools faidx " + reference_fasta + " " + chrom + ":" + str(start_position) + "-" + str(end_position) + " > " + fasta_path
			os.system(command_line)
			plusAddresses[plusAdd[0]] = [ chrom, start_position, end_position ]
			
		plusAddresses[plusAddressesArr[-1][0]] = [chrom, minAddressesArr[-1][1], int(minAddressesArr[-1][-1]) + 5000]
	
	
	#generating bowtie2 indexes for all genes
	genes = glob.glob(geneDir + '*.fasta')		
	for i in range(len(genes)):	
		fasta_path = genes[i]
		gName = fasta_path.replace(geneDir,'',1)
		gName = gName.replace('.fasta','',1)
		btFilesPathPrefix= geneDir + gName
		command_line = "bowtie2-build " + fasta_path + " " + btFilesPathPrefix + " " + ">/dev/null" + "  &"
		os.system(command_line)
		if i % 20 == 0:
			print(i)
			time.sleep(10)			
	

# run bowtie for all samples reads			
print ("About to run bowtie for each sample")	
sampleDirs = glob.glob(baseDir + '*')

for sampleDir in sampleDirs:
	genes = glob.glob(sampleDir + "/plus/*")
	print (sampleDir)
	for i in range(len(genes)):	
		plusFile = genes[i]
		tmpdir = sampleDir + "/plus/"
		gName = plusFile.replace(tmpdir,'',1)
		btFilesPathPrefix = geneDir + gName
		samFile = sampleDir + "/plus/" + gName + '.sam'
		samStat = samFile + '.stat'
		command_line1 = "bowtie2 --mp 2,1 -D 15 -R 2 -N 1 -L 10 -I S,1,0.75 " + " -x " + btFilesPathPrefix + " -r " + plusFile + " > " + samFile + " 2>" + samStat + "  &"
		if os.path.isfile(btFilesPathPrefix+".1.bt2"):
			os.system(command_line1)	
		else:
			print(gName + " does not have an index")			
		if i % 20 == 0:
			print(i)
			time.sleep(10)

exit()
		
