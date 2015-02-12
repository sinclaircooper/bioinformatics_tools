#! /usr/bin/python

import sys 
import re 

'''
A script to annotate minicircle seqeunces that have been tested for circularity and re-oriented according to the sense
and position of CSB3. It will output a gff3 format feature file. 
The features annotated with this script are; circularity, CSB3 sequence and inverted repeats

USEAGE 
This script reads fasta seqeunces from stdin, sequences must be all on oneline.
it expects a list of contig headers (one header per line) that have passed the circular test.

cat minicircle_sequences.fa |fasta_1line| annotate_minicircles.py circular_contigs.txt
'''

circular_2 = []
with open(sys.argv[1]) as f:
    circular = f.readlines()

mini_id = ''
mini_id_l = []
out_string = ''
csb3_G = re.compile('GGGGTTGGTGT')
csb3_A= re.compile('GGGGTTGATGT') 
alt_csb3_RC=re.compile('ACACAAATCCC')
alt_csb3=re.compile('GGGATTTGTGT')
csb1 = re.compile('[A|G][G|A|T]GGGCGT[T|G]C')
csb2 = re.compile('[T|C]C[C|A]CGT[T|G|A]C')
invt_rpt_fwd = re.compile('TAATA[G|A]ATA')
invt_rpt_rev = re.compile('TAT[T|C]TATTA') 
unusual_CSB3=re.compile('AGGGTTGGTGT')
 
for s in circular:
	circular_2.append(s.strip()[1:])


print '##gff-version 3'
for line in sys.stdin:
	if '>' in line: 
		mini_id = line[1:-1]
		mini_id_l.append(mini_id[:-1])
	else: 
		if mini_id in circular_2:
			out_string='%s\t.\tregion\t1\t%s\t.\t+\t.\tID=%s;Is_circular=true;'%(mini_id,len(line)-1,mini_id)
		else:
			out_string='%s\t.\tregion\t1\t%s\t.\t+\t.\tID=%s'%(mini_id,len(line),mini_id)
		print out_string		
		for m in csb3_G.finditer(line):
			print '%s\t.\tmisc_feature\t%s\t%s\t.\t.\t.\tName=CSB3_G'%(mini_id,m.start()+1,m.end()+1)
		for m in unusual_CSB3.finditer(line):
			print '%s\t.\tmisc_feature\t%s\t%s\t.\t.\t.\tName=unusual_CSB3'%(mini_id,m.start()+1,m.end()+1)
		for m in csb3_A.finditer(line):
			print '%s\t.\tmisc_feature\t%s\t%s\t.\t.\t.\tName=CSB3_A'%(mini_id,m.start()+1,m.end()+1)
		for m in csb1.finditer(line):
			print '%s\t.\tmisc_feature\t%s\t%s\t.\t.\t.\tName=CSB1'%(mini_id,m.start()+1,m.end()+1)
		for m in csb2.finditer(line):
			print '%s\t.\tmisc_feature\t%s\t%s\t.\t.\t.\tName=CSB2'%(mini_id,m.start()+1,m.end()+1)
		for m in invt_rpt_fwd.finditer(line):
			print '%s\t.\trpt_type\t%s\t%s\t.\t.\t.\tName=Inverted_repeat_fwd'%(mini_id,m.start()+1,m.end()+1)
		for m in invt_rpt_rev.finditer(line):
			print '%s\t.\trpt_type\t%s\t%s\t.\t.\t.\tName=Inverted_repeat_rev'%(mini_id,m.start()+1,m.end()+1)


