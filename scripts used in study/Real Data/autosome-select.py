#!/usr/bin/env python
from optparse import OptionParser
import os

usage = "usage: python %prog <SAM formatted data with MD field present from stdin> [options] "
parser = OptionParser(usage=usage, version="%prog v0.10")
parser.add_option("-f", "--filename", action="store", type="string", dest="filename",help="filename of file that is processed")
parser.add_option("-o", "--out", action="store", type="string", dest="out",help="filename of output")
(options, args) = parser.parse_args()

wdpath = os.getcwd()
fname = wdpath + '/' + options.filename
f = open(fname, 'r')

f2name = wdpath + '/' + options.out
f2 = open(f2name, 'w')

for line in f:
    if line[0] == '#':
        f2.write(line)
        continue

    col = line.split("\t")
    chr = col[0]
    if chr.isdigit():
        f2.write(line)

f.close()
f2.close()