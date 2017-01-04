#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
import numpy as np
import math
import re
import os

usage = "usage: python %prog <pileup file and damageinformation file> [options] "
parser = OptionParser(usage=usage, version="%prog v1.0")
parser.add_option("-f", "--filename", action="store", type="string", dest="filename",help="filename of file that is processed")
parser.add_option("-d", "--damageinfo", action="store", type="string", dest="damagefilename",help="filename of file that contains damage infomation")
parser.add_option("-o", "--outputprefix", action="store", type="string", dest="outputprefix",help="the prefix of output file")
parser.add_option("-a", "--allsite", action="store_true", dest="allsite", help="output snp of all site", default=False)
parser.add_option("-q", "--quality", action="store", type="float", dest="quality",help="threshold of variant quality",default=20.0)

(options, args) = parser.parse_args()

wdpath = os.getcwd()

rr = 0.001
homo = (1 - 0.001)/2.0
FORMAT = 'GT:AD:DP:GQ:PL'

#   calculate the substitution rate
fname = wdpath + '/' + options.damagefilename
f = open(fname, 'r')
c2t_rate = []
g2a_rate = []
for line in f:
    col = line.split('\t')

    if not col[0].isdigit():
        continue
    if int(col[0]) == 0:
        break
    if int(col[1]) != 0:
        c2t_rate.append(int(col[4])/int(col[1]) + 0.0001)
    else:
        c2t_rate.append(0.0001)
    if int(col[7]) != 0:
        g2a_rate.append(int(col[10])/int(col[7]) + 0.0001)
    else:
        g2a_rate.append(0.0001)
f.close()

#   function of calculate error rate
def error(b):
    return pow(10.0, (33.0 - ord(b))/10.0)

#   function of calculating probability of 1-genotype
def prob_cal(p_in):
    return sum(p_GD[:p_in]) + sum(p_GD[p_in+1:])

f1name = wdpath + '/' + options.filename
f1 = open(f1name, 'r')

f2name = wdpath + '/' + options.outputprefix + '.vcf'
f2 = open(f2name, 'w')

for line in f1:
    col = line.split('\t')
    chromosome = col[0]
    position = col[1]
    ref = col[2]
    count = int(col[3])
    if count < 2:
        continue

    Nt = col[4]
    quals = col[5]
    pos5 = re.findall('\d+', col[6])
    pos3 = re.findall('\d+', col[7])
    Mquals = re.findall('\d+', col[8])
    p = 0.01

    genotype = ['AA', 'AG', 'AC', 'AT', 'GG', 'GC', 'GT', 'CC', 'CT', 'TT']
#   calculate prior probability of a genotype p(D)
    p_G = np.array([homo, rr, rr, rr, homo, rr, rr, homo, rr, homo])

    A_count = 0
    G_count = 0
    C_count = 0
    T_count = 0

#   calculate p(D|G)
    p_DG = np.array([1.0]*10)
    for n in range(count):
        r = error(quals[n]) + 0.01 + 1.0/pow(10.0, int(Mquals[n])/10.0)
        q = g2a_rate[int(pos3[n]) - 1]
        p = c2t_rate[int(pos5[n]) - 1]
        if Nt[n] == 'A':
            A_count += 1
            pA = 1.0 - r
            pG = (1.0 - r)*q + (1.0 - q)*r/3.0
            pC = r/3.0
            pT = r/3.0
            p_DG *= np.array([pA, 0.5*pA + 0.5*pG, 0.5*pA + 0.5*pC, 0.5*pA + 0.5*pT, pG, 0.5*pG + 0.5*pC, 0.5*pG + 0.5*pT, pC, 0.5*pC + 0.5*pT, pT])
        elif Nt[n] == 'G':
            G_count += 1
            pA = r/3.0
            pG = (1.0 - q)*(1.0 - r) + r*q/3.0
            pC = r/3.0
            pT = r/3.0
            p_DG *= np.array([pA, 0.5*pA + 0.5*pG, 0.5*pA + 0.5*pC, 0.5*pA + 0.5*pT, pG, 0.5*pG + 0.5*pC, 0.5*pG + 0.5*pT, pC, 0.5*pC + 0.5*pT, pT])
        elif Nt[n] == 'C':
            C_count += 1
            pA = r/3.0
            pG = r/3.0
            pC = (1.0 - p)*(1.0 - r) + r*p/3.0
            pT = r/3.0
            p_DG *= np.array([pA, 0.5*pA + 0.5*pG, 0.5*pA + 0.5*pC, 0.5*pA + 0.5*pT, pG, 0.5*pG + 0.5*pC, 0.5*pG + 0.5*pT, pC, 0.5*pC + 0.5*pT, pT])
        elif Nt[n] == 'T':
            T_count += 1
            pA = r/3.0
            pG = r/3.0
            pC = p*(1.0 - r) + r*(1.0 - p)/3.0
            pT = 1.0 - r
            p_DG *= np.array([pA, 0.5*pA + 0.5*pG, 0.5*pA + 0.5*pC, 0.5*pA + 0.5*pT, pG, 0.5*pG + 0.5*pC, 0.5*pG + 0.5*pT, pC, 0.5*pC + 0.5*pT, pT])
        else:
            continue

    p_GD = (p_G*p_DG)/sum(p_G*p_DG)
    if 0 in p_GD:
        for i in range(len(p_GD)):
            if p_GD[i] == 0:
                p_GD[i] = 1.000000000000017e-308

    max_pos = p_GD.argmax()
    alt = genotype[max_pos]

    if ref == 'A':
        if alt == 'AA' and options.allsite:
            ALT = '.'
            Qual = round(-10.0 * math.log(prob_cal(max_pos), 10.0),2)
            INFO = 'AC=0;AF=0.00;AN=2;DP=' + str(count)
            PL = int(round(-10.0*math.log(p_GD[max_pos], 10.0)))
            Sample1 = '0/0:' + str(A_count) + ',' + str(G_count + C_count + T_count) + ':' + str(count) + ':' +str(Qual)+ ':' + str(PL)
            f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
        if alt != 'AA':
            Qual = round(-10.0 * math.log(p_GD[0], 10.0),2)
            if Qual < options.quality:
                continue
            if ref in alt:
                INFO = 'AC=1;AF=0.500;AN=2;DP=' + str(count)
                if alt == 'AG':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str((int(round(-10.0 * math.log(prob_cal(1), 10.0)))))
                    Sample1 = '0/1:' + str(A_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'AC':
                    ALT = 'C'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(2), 10.0))))
                    Sample1 = '0/1:' + str(A_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'AT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(3), 10.0))))
                    Sample1 = '0/1:' + str(A_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            elif alt[0] == alt[1]:
                INFO = 'AC=2;AF=1.00;AN=2;DP=' + str(count)
                if alt == 'GG':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(4), 10.0))))
                    Sample1 = '1/1:' + str(A_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CC':
                    ALT = 'C'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(7), 10.0))))
                    Sample1 = '1/1:' + str(A_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'TT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(9), 10.0))))
                    Sample1 = '1/1:' + str(A_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            else:
                INFO = 'AC=2;AF=0.500,0.500;AN=2;DP=' + str(count)
                if alt == 'GC':
                    ALT = 'G,C'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(5), 10.0))))
                    Sample1 = '1/2:' + str(A_count) + ',' + str(G_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GT':
                    ALT = 'G,T'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(6), 10.0))))
                    Sample1 = '1/2:' + str(A_count) + ',' + str(G_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CT':
                    ALT = 'C,T'
                    PL = str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(8), 10.0))))
                    Sample1 = '1/2:' + str(A_count) + ',' + str(C_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

    elif ref == 'G':
        if alt == 'GG' and options.allsite:
            ALT = '.'
            Qual = round(-10.0 * math.log(prob_cal(max_pos), 10.0),2)
            INFO = 'AC=0;AF=0.00;AN=2;DP=' + str(count)
            PL = int(round(-10.0*math.log(p_GD[max_pos], 10.0)))
            Sample1 = '0/0:' + str(G_count) + ',' + str(A_count + C_count + T_count) + ':' + str(count) + ':' + str(Qual) + ':' + str(PL)
            f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
        if alt != 'GG':
            Qual = round(-10.0 * math.log(p_GD[4], 10.0),2)
            if Qual < options.quality:
                continue
            if ref in alt:
                INFO = 'AC=1;AF=0.500;AN=2;DP=' + str(count)
                if alt == 'AG':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(1), 10.0))))
                    Sample1 = '0/1:' + str(G_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GC':
                    ALT = 'C'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(5), 10.0))))
                    Sample1 = '0/1:' + str(G_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(6), 10.0))))
                    Sample1 = '0/1:' + str(G_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            elif alt[0] == alt[1]:
                INFO = 'AC=2;AF=1.00;AN=2;DP=' + str(count)
                if alt == 'AA':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(0), 10.0))))
                    Sample1 = '1/1:' + str(G_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CC':
                    ALT = 'C'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(7), 10.0))))
                    Sample1 = '1/1:' + str(G_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'TT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(9), 10.0))))
                    Sample1 = '1/1:' + str(G_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            else:
                INFO = 'AC=2;AF=0.500,0.500;AN=2;DP=' + str(count)
                if alt == 'AC':
                    ALT = 'A,C'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(2), 10.0))))
                    Sample1 = '1/2:' + str(G_count) + ',' + str(A_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'AT':
                    ALT = 'A,T'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(3), 10.0))))
                    Sample1 = '1/2:' + str(G_count) + ',' + str(A_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CT':
                    ALT = 'C,T'
                    PL = str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(8), 10.0))))
                    Sample1 = '1/2:' + str(G_count) + ',' + str(C_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

    elif ref == 'C':
        if alt == 'CC' and options.allsite:
            ALT = '.'
            Qual = round(-10.0 * math.log(prob_cal(max_pos), 10.0),2)
            INFO = 'AC=0;AF=0.00;AN=2;DP=' + str(count)
            PL = int(round(-10.0*math.log(p_GD[max_pos], 10.0)))
            Sample1 = '0/0:' + str(C_count) + ',' + str(A_count + G_count + T_count) + ':' + str(count) + ':' + str(Qual) + ':' + str(PL)
            f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
        if alt != 'CC':
            Qual = round(-10.0 * math.log(p_GD[7], 10.0),2)
            if Qual < options.quality:
                continue
            if ref in alt:
                INFO = 'AC=1;AF=0.500;AN=2;DP=' + str(count)
                if alt == 'AC':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(2), 10.0))))
                    Sample1 = '0/1:' + str(C_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GC':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(5), 10.0))))
                    Sample1 = '0/1:' + str(C_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(8), 10.0))))
                    Sample1 = '0/1:' + str(C_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            elif alt[0] == alt[1]:
                INFO = 'AC=2;AF=1.00;AN=2;DP=' + str(count)
                if alt == 'AA':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(0), 10.0))))
                    Sample1 = '1/1:' + str(C_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GG':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(4), 10.0))))
                    Sample1 = '1/1:' + str(C_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'TT':
                    ALT = 'T'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(9), 10.0))))
                    Sample1 = '1/1:' + str(C_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
            else:
                INFO = 'AC=2;AF=0.500,0.500;AN=2;DP=' + str(count)
                if alt == 'AG':
                    ALT = 'A,G'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(1), 10.0))))
                    Sample1 = '1/2:' + str(C_count) + ',' + str(A_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'AT':
                    ALT = 'A,T'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(3), 10.0))))
                    Sample1 = '1/2:' + str(G_count) + ',' + str(A_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GT':
                    ALT = 'G,T'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(6), 10.0))))
                    Sample1 = '1/2:' + str(C_count) + ',' + str(G_count) + ',' + str(T_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

    elif ref == 'T':
        if alt == 'TT' and options.allsite:
            ALT = '.'
            Qual = round(-10.0 * math.log(prob_cal(max_pos), 10.0),2)
            INFO = 'AC=0;AF=0.00;AN=2;DP=' + str(count)
            PL = int(round(-10.0*math.log(p_GD[max_pos], 10.0)))
            Sample1 = '0/0:' + str(T_count) + ',' + str(A_count + G_count + C_count) + ':' + str(count) + ':' +str(Qual)+ ':' + str(PL)
            f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
        if alt != 'TT':
            Qual = round(-10.0 * math.log(p_GD[9], 10.0),2)
            if Qual < options.quality:
                continue
            if ref in alt:
                INFO = 'AC=1;AF=0.500;AN=2;DP=' + str(count)
                if alt == 'AT':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(3), 10.0))))
                    Sample1 = '0/1:' + str(T_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GT':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(-10.0*round(math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(6), 10.0))))
                    Sample1 = '0/1:' + str(T_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CT':
                    ALT = 'C'
                    PL = str(int(-10.0*round(math.log(p_GD[9], 10.0)))) + ',' + str(int(-10.0*round(math.log(p_GD[8], 10.0)))) + ',' + str(int(-10.0*round(math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(8), 10.0))))
                    Sample1 = '0/1:' + str(T_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

            elif alt[0] == alt[1]:
                INFO = 'AC=2;AF=1.00;AN=2;DP=' + str(count)
                if alt == 'AA':
                    ALT = 'A'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(0), 10.0))))
                    Sample1 = '1/1:' + str(T_count) + ',' + str(A_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GG':
                    ALT = 'G'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(4), 10.0))))
                    Sample1 = '1/1:' + str(T_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'CC':
                    ALT = 'C'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(7), 10.0))))
                    Sample1 = '1/1:' + str(T_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                else:
                    continue
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
            else:
                INFO = 'AC=2;AF=0.500,0.500;AN=2;DP=' + str(count)
                if alt == 'AG':
                    ALT = 'A,G'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[1], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(1), 10.0))))
                    Sample1 = '1/2:' + str(T_count) + ',' + str(A_count) + ',' + str(G_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'AC':
                    ALT = 'A,C'
                    PL = str(int(round(-10.0*math.log(p_GD[9], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[3], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[0], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[2], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[7], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(2), 10.0))))
                    Sample1 = '1/2:' + str(T_count) + ',' + str(A_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                elif alt == 'GC':
                    ALT = 'G,C'
                    PL = str(int(round(-10.0*math.log(p_GD[7], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[5], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[4], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[8], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[6], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[9], 10.0))))
                    GQ = str(int(round(-10.0 * math.log(prob_cal(5), 10.0))))
                    Sample1 = '1/2:' + str(T_count) + ',' + str(G_count) + ',' + str(C_count) + ':' + str(count) + ':' + GQ + ':' + str(PL)
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
f1.close()
f2.close()

