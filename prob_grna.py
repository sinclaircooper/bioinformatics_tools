#! /usr/bin/python 

import sys
import re 
import math 
import operator
import pickle
from pypeaks import Data, Intervals
import matplotlib.pyplot as plt
import random

'''
Takes a list of annotated minicircles and a set of scores for upstream and downstream nucleotide
scores (provided but nucleotide scores specific to your data set can be easily generated) and 
predicts the start and stop positions for gRNAs.

Useage
cat minicricles.fa | fasta_1line | prob_grna.py features_list.gff > grna_start_stop.gff

'''


up_dict={}
down_dict={}
header=''
list_probs_shuf={}
list_probs_shuf_2={}
background=[]


#calculate a 2 dimensional overlap 
def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

#get a substring
def getsubstring(start, stop, string):
	return string[start:stop]


#shuffle a string
def Shuffle(string):
	chars  = list(string)
	random.shuffle(chars)
	return ''.join(chars)

#calulate scores for the minicircle shuffled maxy times
def score_shuffled(seqs,maxy, jj):
	for line in seqs:
		s=0
		while s<=len(line)-160:
			cassette=line[s:s+160]
			score=score_cassette(cassette, jj)
			if s in list_probs_shuf:
				list_probs_shuf[s].append(score)
			else:
				list_probs_shuf[s]=[score]
			s+=1
		for key in list_probs_shuf:
			list_probs_shuf_2[key] = sum(list_probs_shuf[key])/(maxy+1)
	return list_probs_shuf_2.values()

#mean function
def m(values):  
        size = len(values)  
        s = 0.0  
        for n in range(0, size):  
                s += values[n]  
        return s / size;  

#std deviation
def SD(values, mean):  
        size = len(values)  
        s = 0.0  
        for n in range(0, size):  
		s += ((values[n] - mean)**2)
	return math.sqrt((1.0/(size-1))*(s))

#get the upstream+downstream scores (calculated from the alignments)
with open('nucleotide_scores.up') as n:
	up_scores=n.readlines()
with open('nucleotide_scores.down') as n:
	down_scores=n.readlines()
#get the annotations from alignments for plotting later 
with open(sys.argv[1]) as n:
	features=n.readlines()

f1=[]
f2=[]
f3=[]
f4=[]
inverted_repeats1=[]
inverted_repeats2=[]
inverted_repeats3=[]
inverted_repeats4=[]
#get the feature annotations
for f in features:
	if '#' in f:continue
	if 'CSB' in f:continue
	#if 'Inverted' in f:continue 
		#f=f.split()
		#inverted_repeats1.append(f[0])
		#inverted_repeats2.append(f[3])
		#inverted_repeats3.append(f[4])
		#inverted_repeats4.append(f[8])
	if 'circular' in f:continue 
	if 'region' in f:continue 
	else:
		f=f.split()
		#feature_dict[f[0].strip('\n')+'{'+str(count)+'}']=f[3],f[4],f[8]
		f1.append(f[0])
		f2.append(f[3])
		f3.append(f[4])
		f4.append(f[8])

Features=zip(f1,f2,f3,f4)

#get the scores from the appropriate files and populate a dictionary
for s in up_scores:
	s=s.split()
	up_dict[int(s[0])]=[math.log(float(i)+0.001) for i in s[1:]]
for s in down_scores:
	s=s.split()
	down_dict[int(s[0])]=[math.log(float(i)+0.001) for i in s[1:]]



#this is the function that applies the scoring matrices to a 160 bp section of the minicircle
def score_cassette(cassette, type):
        base = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        if type == 0: 
                return sum([up_dict[i][base[b]] for i,b in enumerate(cassette)])
        elif type == 1: 
                return sum([down_dict[i][base[b]] for i,b in enumerate(cassette)])
        elif type == 2: 
                return -sum([up_dict[i][base[b]] for i,b in enumerate(cassette)]) + \
                        sum([down_dict[i][base[b]] for i,b in enumerate(cassette)])









maxy=50
peaks1=[]
peaks2=[]
peaks3=[]
peaks4=[]
genes=re.compile('rps12|nd9|nd8|nd7|nd3|a6|murf2|cyb|cr4|cr3|co3|co2')
NB_scores_for_aln_up=[]
NB_scores_for_aln_dwn=[]
bk_up_pos,bk_down_pos=[],[]

for line in sys.stdin:
	cass=[]
	repeats=[]
	data2_pos=[]
	data3_pos=[]
	if '>' in line:header=line.strip('>')
	else:
		x=0
		data1,data2,data3=[],[],[]
		while x<len(line)-160:
			cassette=line[x:x+160]
			score1=score_cassette(cassette,0)
			score2=score_cassette(cassette,1)
			data1.append(x+80)
	               	data2.append(score1)
			data3.append(score2)
			x+=1
		for d in data2:
			d2_min=min(data2)
			data2_pos.append(d-d2_min)
		for d in data3:
			d3_min=min(data3)
			data3_pos.append(d-d3_min)
		header=header.strip('\n')
			
		
		
		for f in Features:
			if header==f[0]:
				if 'Inverted' in f[3]:
					repeats.append([f[0],f[3],f[1],f[2]])

		for position, item in enumerate(repeats):
			if int(position)+2>len(repeats):continue
			elif item[1]=='Name=Inverted_repeat_fwd' and repeats[position+1][1]=='Name=Inverted_repeat_rev':
				cass.append([item[0],repeats[position+1][2],item[3]])
		upstream_data=zip(data1,data2_pos)
		downstream_data=zip(data1,data3_pos)

		start_prob_c=[]
		stop_prob_c=[]

		for c in cass:
			if header==c[0]:
				cass_up=[]
				for u in upstream_data:
					if int(u[0])>int(c[1]) and int(u[0])<int(c[2]):
						cass_up.append(u)
				if len(cass_up)<1:continue 
				else:
					start_prob_c.append(max(cass_up,key=lambda x:x[1]))
				cass_down=[]					
				for d in downstream_data:
					if int(d[0])>int(c[1]) and int(d[0])<int(c[2]) and int(d[0])>int(max(cass_up,key=lambda x:x[1])[0]):
						cass_down.append(d)
				if len(cass_down)<1:continue
				else:
					stop_prob_c.append(max(cass_down,key=lambda x:x[1]))
					print header+'\t.\tmRNA\t'+str(max(start_prob_c,key=lambda x:x[1])[0])+'\t'+\
								str(max(stop_prob_c,key=lambda x:x[1])[0])+'\t.\t.\t.\tName=NB_prob'+'_'+\
								str(int(max(start_prob_c,key=lambda x:x[1])[1]))+'_'+\
								str(int(max(stop_prob_c,key=lambda x:x[1])[1]))
				start_prob_c=[]
				stop_prob_c=[]

			


