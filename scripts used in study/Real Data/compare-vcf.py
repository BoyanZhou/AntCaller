#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
import os

usage = "usage: python %prog  [options] "
parser = OptionParser(usage=usage, version="%prog v0.10")
parser.add_option("-f", "--filename", action="store", type="string", dest="filename",help="filename of file that is processed")
parser.add_option("-d", "--dbsnp", action="store", type="string", dest="dbsnp",help="dbsnpvcf")
parser.add_option("-p", "--prefix", action="store", type="string", dest="prefix",help="prefix of output files")
(options, args) = parser.parse_args()

wdpath = os.getcwd()
fname = wdpath + '/' + options.filename
f = open(fname, 'r')

f2name = wdpath + '/' + options.dbsnp
f2 = open(f2name, 'r')

f3name = wdpath + '/' + options.prefix + '.inright'
f3 = open(f3name, 'w')

f4name = wdpath + '/' + options.prefix + '.inwrong'
f4 = open(f4name, 'w')

f5name = wdpath + '/' + options.prefix + '.out'
f5 = open(f5name, 'w')

total = 0
inright = 0
inwrong = 0
out = 0

chr2 = '0'
pos2 = '0'
ref2 = 'N'
alt2 = 'N'

for line in f:
    if line[0] == '#':
        continue
    total += 1
    col = line.split('\t')
    chr = col[0]
    if not chr.isdigit():
        break
    pos = int(col[1])
    ref = col[3]
    alt = col[4][0]
    if not chr2.isdigit():
        out += 1
        f5.write(line)
        continue
    if int(chr) > int(chr2):
        while int(chr) > int(chr2):
            line2 = f2.readline()
            if line2[0] == '#':
                continue
            col2 = line2.split('\t')
            chr2 = col2[0]
            pos2 = int(col2[1])
            ref2 = col2[3]
            alt2 = col2[4].split(',')
    if int(chr) < int(chr2):
        out += 1
        f5.write(line)
    if int(chr) == int(chr2):
        if pos == pos2:
            if alt in alt2:
                inright += 1
                f3.write(line)
            else:
                inwrong += 1
                f4.write(line)
            continue
        if pos < pos2:
            out += 1
            f5.write(line)

        while pos > pos2:
            line2 = f2.readline()
            if line2[0] == '#':
                continue
            col2 = line2.split('\t')
            chr2 = col2[0]
            if not chr2.isdigit():
                out += 1
                f5.write(line)
                break
            if int(chr) < int(chr2):
                out += 1
                f5.write(line)
                break

            pos2 = int(col2[1])
            ref2 = col2[3]
            alt2 = col2[4].split(',')
            if pos == pos2:
                if alt in alt2:
                    inright += 1
                    f3.write(line)
                else:
                    inwrong += 1
                    f4.write(line)
            if pos < pos2:
                out += 1
                f5.write(line)

print(total, inright, inwrong, out)
f.close()
f2.close()
f3.close()
f4.close()
f5.close()


















