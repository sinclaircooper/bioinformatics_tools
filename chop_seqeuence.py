#! /usr/bin/python 

import sys 
import re 

'''
Takes any seqeunce and chops it in half exactly.
Outputs each sequnce in separate files in the format
>seq_id_1st_half
CATACGTCATCATCAGCAT
>seq_id_2nd_half
ATCAGCATCGACTCAGATC

These seqeuneces can then be batch processed (i.e. using cap3 to check for circularity).
Expects sequnces in 1 line fasta format.

Usage:
cat seqsin.fa |fasta_1line| chop_sequence.py

'''

headers = []
currentseq = ''

for line in sys.stdin:
	if '>' in line:currentseq = line[:-1]
	
	else:
		currentseq2 = currentseq.replace(">","")+".chp"
		with open(currentseq2,'w') as fileout:		
			fileout.write(currentseq+'_1st_half'+'\n'+line[:(len(line)/2)]+'\n'+currentseq+'_2nd_half'+'\n'+line[-len(line)/2:-1])
