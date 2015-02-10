#!/usr/bin/python

#geneerates a contig converage histogram from SAM output 

import string
import sys
import re  
import numpy as np
import matplotlib.pyplot as plt
import math

def m(values):  
        size = len(values)  
        s = 0.0  
        for n in range(0, size):  
                s += values[n]  
        return s / size; 


def SD(values, mean):  
        size = len(values)  
        s = 0.0  
        for n in range(0, size):  
		s += ((values[n] - mean)**2)
	return math.sqrt((1.0/(size-1))*(s))

number_of_reads = []

for line in sys.stdin: 
	number_of_reads.append(int(line[:-1]))

mean = m(number_of_reads)
S_D = SD(number_of_reads,m(number_of_reads))

#plot_title = sys.argv[1]
y=range(0,100,1)

plt.hist(number_of_reads, bins = y)
#plt.title(plot_title)
plt.xlabel('Number of reads')
plt.ylabel('Contig frequency')
plt.axvline(x=mean, color='r')
#plt.axvline(x=2456, color='g')
#plt.axvline(x=29022,color='g')
plt.axvline(x=mean+S_D*2, color='c')
#lines = plt.plot([15739,10],[10,15739])
plt.show()

#int(round(max(number_of_reads),3))
