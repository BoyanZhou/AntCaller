#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
import os
import operator
tianyuan = [73, 263, 750, 1438, 2706, 4769, 5348, 5836, 7028, 8860, 11257, 11719, 14766, 15326, 16189, 16293]
wdpath = os.getcwd()

usage = "usage: python %prog <SAM formatted data with MD field present from stdin> [options] "
parser = OptionParser(usage=usage, version="%prog v0.10")
parser.add_option("-n", "--number", action="store", type="int", dest="number",help="number of snp be called")
(options, args) = parser.parse_args()

def compare(x):
    data = []
    true = 0
    for line in x:
        if line[0] == "#":
            continue
        col = line.split("\t")
        site = []
        site.append(int(col[1]))
        site.append(float(col[5]))
        data.append(site)
    data.sort(key=operator.itemgetter(1), reverse=True)
    if len(data) >= options.number:
        n = options.number
        for i in range(options.number):
            if data[i][0] in tianyuan:
                true += 1
    else:
        n = len(data)
        for i in data:
            if i[0] in tianyuan:
                true += 1
    return [n, true]

fsam_name = wdpath + '/' + 'simu_Tianyuan_MT.damage.sort.rmdup.samtools.vcf'
fsam = open(fsam_name, 'r')
sam = compare(fsam)
fsam.close()

fgatk_name = wdpath + '/' + 'simu_Tianyuan_MT.damage.sort.rmdup.group.gatk.vcf'
fgatk = open(fgatk_name, 'r')
gatk = compare(fgatk)
fgatk.close()

fmy_name = wdpath + '/' + 'simu_Tianyuan_MT.ant.vcf'
fmy = open(fmy_name, 'r')
my = compare(fmy)
fmy.close()

name = wdpath + '/' + 'result' +str(options.number) + '.txt'
f = open(name, 'a')
f.write(str(sam[0]) + "\t" + str(sam[1]) + "\t" +str(gatk[0]) + "\t" +str(gatk[1]) + "\t" +str(my[0]) + "\t" +str(my[1]) + "\n")
f.close()


