#! /usr/bin/python

# this is a script to take a set of contigs and generate a set of 'reads' from them
#the general approach is to start at base 1 and take a 100 bp snapshot of the 1kb sequence
# the window will then move 1 bp up the seq and take another snapshot

import string
import sys
import re 
from string import zfill 
#import random

current_header = ''


def getSubstring(string, start, stop):
	return string[start:stop]



for line in sys.stdin:
	if '>' in line: current_header = line[:-1] 
	else:
		x=0;
		counter = 0;
		while x<=len(line)-100:
			if getSubstring(line,x,x+100)=='':continue
			else:
				print current_header+'_'+str(counter)+"\n"+getSubstring(line,x,x+100)
				x=x+5;
				counter +=1;



		
				
				 


