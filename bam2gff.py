#! /usr/bin/python

import sys 
import re 
import os 

'''
Reads alingments of gRNAs made to minicircles in sam format and outputs a feature file (GFF3)
Expects you to tell it weather these gRNA features were predicted from assembled minicircles or
from NGS reads; formats the output accordingly.

USEAGE
samtools view in.bam | bam2gff.py [read|align] > out.gff 
'''



read_or_align=(sys.argv[1])

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


#print 'print ##gff-version 3'
for line in sys.stdin:
	line = line.split()
	if '*' in line[2]: continue
	else:
		if int(line[1])==16:
			sense='+'
		elif int(line[1])==0:
			sense='-'
		elif int(line[1])==272:
			sense='-'
		if read_or_align=='read':
			cigar='_'+line[5]
		elif read_or_align=='align':
			cigar=''
		elif read_or_align=='NB':
			cigar='_NB'
		if read_or_align=='RNA':
			print line[2]+"\t.\tmRNA\t"+str(line[3])+"\t"+str(int(line[3])\
				+len(line[9][:-1]))+"\t.\t"+sense+"\t.\t"+"Name="+line[0]
		#if line[0].split("_")[1]+"_"+line[0].split("_")[2]!=line[2]:continue
		else:
			print line[2]+"\t.\tmRNA\t"+str(line[3])+"\t"+str(int(line[3])+len(line[9][:-1]))\
				+"\t.\t"+sense+"\t.\t"+"Name="+find_between(line[0],"Tb","ed")+cigar+"_rd"



