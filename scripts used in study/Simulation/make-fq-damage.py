#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
import random
import os

usage = "usage: python %prog <fastq file> [options] "
parser = OptionParser(usage=usage, version="%prog v0.10")
parser.add_option("-f", "--filename", action="store", type="string", dest="filename",help="processed file")
parser.add_option("-o", "--output", action="store", type="string", dest="output",help="output filename")
parser.add_option("-r", "--rate", action="store", type="float", dest="rate",help="substitution rate of the end of reads")
(options, args) = parser.parse_args()

wdpath = os.getcwd()
outputfile = wdpath + '/' + options.output
#   calculate the substitution rate
c2t_rate = []
p = options.rate
for i in range(60):
    D = p * pow(1 - p, i) + 0.001
    c2t_rate.append(D)
g2a_rate = c2t_rate[:]

f2 = open(options.filename, 'r')
fout = open(outputfile, 'w')
while f2:
    line1 = f2.readline()
    if not line1:
        break
    fout.write(line1)
    if line1[0] == '@':
        line2 = list(f2.readline().strip('\n'))
        for i in range(len(line2)):
            if line2[i] == 'C':
                if random.uniform(0, 1) < c2t_rate[i]:
                    line2[i] = 'T'

            elif line2[i] == 'G':
                if random.uniform(0, 1) < c2t_rate[len(line2) - i - 1]:
                    line2[i] = 'A'
        line2 = ''.join(line2)
        line3 = f2.readline()
        line4 = f2.readline()
        fout.write(line2 + '\n' + line3 + line4)

f2.close()
fout.close()




