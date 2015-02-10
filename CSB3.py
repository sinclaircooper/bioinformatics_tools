#!/usr/bin/python

#searches for minicircle specific Conserved seqeunce 3 in a fasta seqeunces 
#only returns sequences that contain csb3

import sys 
import re


current_name=''
current_seq=''
mydict = {}

for line in sys.stdin:
	if any ('>' in l for l in line):
		if current_name != '':
			mydict[current_name] = current_seq
		current_name = line[0:].rstrip('\n')
		current_seq = ''
	else:
		current_seq = current_seq + line.rstrip('\n')

#csb3 = '(G|T)G(G|A)TT(G|T)(T|G|A)(T|C)G'
csb3 = 'GGGGTTGGTGTA'
contig_sequences =mydict.values()

current_key = ''
current_seq = ''
new_mini_circles={}

for current_key, current_seq in mydict.iteritems():
	if re.search(csb3,current_seq):# in current_seq:
		print current_key,'\n',current_seq
		new_mini_circles[current_key] = current_seq 
