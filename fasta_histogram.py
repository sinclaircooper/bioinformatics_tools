#!/usr/bin/python

import string
import sys
import re  
import numpy as np
import matplotlib.pyplot as plt

lengthseq = []

for contig in sys.stdin:
	if '>' in contig: continue
	else:
		contig_length = len(contig)
		lengthseq.append(contig_length)
count=0

for l in lengthseq:
	if l > 900:
		count=count+1
print min(lengthseq)
print max(lengthseq)
#plot_title = sys.argv[1]+'\ncontigs > 900bp='+str(count) 
y=range(0,6000,50)
plt.xlim(0,6000)
plt.hist(lengthseq, bins = y)
#plt.title(plot_title)
plt.xlabel('Contig length (nt)')
plt.ylabel('Frequency')
plt.show()

