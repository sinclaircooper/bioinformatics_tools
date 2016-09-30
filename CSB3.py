#!/usr/bin/python

#searches for minicircle specific Conserved seqeunce 3 in fasta seqeunces 
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
		mydict[current_name] = current_seq

#csb3 = '(G|T)G(G|A)TT(G|T)(T|G|A)(T|C)G'
#csb3 = 'GGGGTTGGTGTA'
#csb3 = 'GGG[G|A]TT[G|T]GTGT|ACAC[A|C]AA[T|C]CCC|CCC[C|T]AA[C|A]CACA|TGTG[T|G]TT[A|G]GGG'
csb3= 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTAGTGT|ACACTAACCCC'
#csb3='GGGATTTGTGT|ACACAAATCCC'
contig_sequences =mydict.values()

#print mydict

current_key = ''
current_seq = ''
new_mini_circles={}

for current_key, current_seq in mydict.iteritems():
	#print current_key,current_seq
	if re.findall(csb3,current_seq):# in current_seq:
		print current_key,'\n',current_seq
		new_mini_circles[current_key] = current_seq 
	
