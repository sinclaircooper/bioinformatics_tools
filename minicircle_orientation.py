#! /usr/bin/python
#GGGGTTGGTGTA
#take a minicircle and put csb3 at the start. also changes them to to the correct sense

'''
Reads input from stdin prints to stdout.
-Will complelent/reverse_complement/reverse the seqeunces according the detected orientation 
of CSB3.
-Will also chop and re-orientate sequeces with CSB3 at the start if they have passed the circular
test. 
-Expects a list of contig IDs that have passed the circularity test.

Useage

cat minicircles.fa | fasta_1line | minicircle_orintation.py passed_circular_test.txt > minicircles.oriented
''' 


import string 
import sys 
import re 

seq = []
mini_id = []
new_seq = []
circ_2 = []
dict_seq = {}
csb3_rev_comp  = re.compile(r'ACACCAACCCC|ACATCAACCCC')
csb3 =           re.compile(r'GGGGTTGGTGT|GGGGTTGATGT|AGGGTTGGTGT')
csb3_comp =      re.compile(r'CCCCAACCACA|CCCCAACTACA')
csb3_rev =       re.compile(r'TGTGGTTGGGG|TGTAGTTGGGG')
alt_csb3_RC =    re.compile(r'ACACAAATCCC')
alt_csb3 =       re.compile(r'GGGATTTGTGT')

big_csb3=re.compile(r'GGGGTTGGTGT|GGGGTTGATGT|AGGGTTGGTGT')#GGGGTTGATGT

'''
csb3_rev_comp  = re.compile(r'ACACCAACCCC|ACACAAATCCC')
csb3 =           re.compile(r'GGGGTTGGTGT|GGGATTTGTGT')
csb3_comp =      re.compile(r'CCCCAACCACA|CCCTAAACACA')
csb3_rev =       re.compile(r'TGTGGTTGGGG|TGTGTTTAGGG')
'''

def complement(sequence):
	complement = {'A':'T','C':'G','G':'C','T':'A', '-':'-'}
	return "".join([complement.get(nt.upper(), '') for nt in sequence])


with open(sys.argv[1]) as f:
    circular = f.readlines()
norm=0
rev=0
rc=0
c=0

current_name=''
current_seq=''
mydict = {}

for line in sys.stdin:
	if any ('>' in l for l in line):
		if current_name != '':
			mydict[current_name] = current_seq
		current_name = line[1:].rstrip('\n')
		current_seq = ''
	else:
		current_seq = current_seq + line.rstrip('\n')
		mydict[current_name] = current_seq

list_seq=mydict.items()


for ls in list_seq:
	if re.findall(csb3,ls[1]): 
		mini_id.append(ls[0])
		seq.append(ls[1])
	elif re.findall(csb3_rev,ls[1]):
		mini_id.append(ls[0]) 
		seq.append(ls[1][::-1])
	elif re.findall(csb3_rev_comp,ls[1]): 
		mini_id.append(ls[0])
		seq.append(complement(ls[1][::-1]))
	elif re.findall(csb3_comp,ls[1]): 
		mini_id.append(ls[0])
		seq.append(complement(ls[1][:-1]))
	elif re.findall(alt_csb3,ls[1]):
		mini_id.append(ls[0])
		seq.append(complement(ls[1][::-1]))
	elif re.findall(alt_csb3_RC,ls[1]):
		mini_id.append(ls[0])
		seq.append(ls[1][:-1])
	#else:
		#print 'fail',ls[0]


list_seq_2=zip(mini_id,seq) 



for s in circular:
	circ_2.append(s.strip())



for ls in list_seq_2:
	if ls[0] in circ_2:
			match=re.search(big_csb3,ls[1])
			chop=match.span()[0]
			print '>'+ls[0]+'\n'+\
				ls[1][chop:]+ls[1][:chop]
	else:
		print  '>'+ls[0]+'\n'+ls[1]


