#! /usr/bin/python

#This script takes a fasta sequence and returns the complemented version of the sequence

import sys

def Complement(sequence):

	complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}

	return "".join([complement.get(nt.upper(), '') for nt in sequence[::1]])

headers = []
reversed_seqs = []

for line in sys.stdin:
	if '>' in line:
		headers.append(line)
	else:
		reversed_seqs.append(Complement(line))

for h,r in zip(headers,reversed_seqs):
	print h+r
