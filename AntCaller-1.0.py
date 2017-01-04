#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
from optparse import OptionGroup
from multiprocessing.dummy import Pool as ThreadPool
import numpy as np
import time
import sys
import math
import os
import re

usage = "usage: python %prog <SAM formatted data with MD field present from stdin> [options] "
parser = OptionParser(usage=usage, version="%prog v1.0")
group = OptionGroup(parser,"pileup options","parameter setting for pileup reads")
group.add_option("--pileup", action="store_true", dest="pileup",help="pileup reads from sam file",default=False)
parser.add_option_group(group)

group = OptionGroup(parser,"extract damage information options","parameter setting for extracting damage information")
group.add_option("--extract", action="store_true", dest="extract",help="extract damage information from sam/bam file",default=False)
group.add_option("-l", "--length", action="store", type="int", dest="length",help="the length of the longest reads", default=300)
parser.add_option_group(group)

group = OptionGroup(parser,"snp calling options","parameter setting for snp calling")
group.add_option("--snpcalling", action="store_true", dest="snpcalling",help="call snp after pileup and damageinfo",default=False)
group.add_option("-f", "--pileupfile", action="store", type="string", dest="pileupfile",help="input pileup file for snp calling")
group.add_option("-d", "--damageinfo", action="store", type="string", dest="damageinfo",help="input the file that stores the damage information")
group.add_option("-p", "--ploidy", action="store", type="int", dest="ploidy",help="ploidy of the genome",default=2)
group.add_option("-a", "--allsite", action="store_true", dest="allsite",help="output snp of all site",default=False)
group.add_option("-q", "--quality", action="store", type="float", dest="quality",help="threshold of variant quality",default=20.0)
parser.add_option_group(group)

parser.add_option("-o", "--outputprefix", action="store", type="string", dest="outputprefix",help="the prefix of output file")
parser.add_option("-t", "--thread", action="store", type="int", dest="thread",help="the number of thread (not for extracting damage information )", default=5)
(options, args) = parser.parse_args()

time_start = time.time()
wdpath = os.getcwd()

#   program of pileup
if options.pileup:
    headname = wdpath + '/' + options.outputprefix + '.header'
    fhead = open(headname, 'w')
    chr_list = []
    fsplit = []
#   split the sam/bam file
    for line in sys.stdin:
        if line[0] == '@':
            fhead.write(line)
            continue
        col = line.split('\t')
        chromosome = col[2]
        if chromosome == '*':
            continue
        if not chr_list:
            chr_list.append(chromosome)
            splitname = wdpath + '/' + options.outputprefix + '.temporary.1'
            fsplit = open(splitname, 'w')
            fsplit.write(line)
            continue
        if chromosome != chr_list[-1]:
            fsplit.close()
            chr_list.append(chromosome)
            splitname = wdpath + '/' + options.outputprefix + '.temporary.' + str(len(chr_list))
            fsplit = open(splitname, 'w')
            fsplit.write(line)
        else:
            fsplit.write(line)
    fhead.close()
    fsplit.close()
    temp_list = range(1,len(chr_list)+1)
    temp_list.reverse()

#   end of spliting the sam/bam file
    def split(para):
        os.system('cat ' + options.outputprefix + '.temporary.' + str(para) + ' | python AntCaller_pileup.py > ' + options.outputprefix + '.pileup.' + str(para))

    pool = ThreadPool(options.thread)
    pool.map(split, temp_list)
    pool.close()
    pool.join()
    temp_list.reverse()

#   combine files
    fpilename = wdpath + '/' + options.outputprefix + '.AntCaller.pileup'
    fpile = open(fpilename, 'w')

#   add header
    fhead = open(headname, 'r')
    for line in fhead:
        if line:
            fpile.write(line)
    fhead.close()

    colname = '@chr' + '\t' + 'position' + '\t' + 'ref' + '\t' + 'num' + '\t' + 'base' + '\t' + 'quality' + '\t' + '3-end' +'\t' + '5-end' + '\t' + 'MQ' + '\n'
    fpile.write(colname)

    for i in temp_list:
        fopenname = wdpath + '/' + options.outputprefix + '.pileup.' + str(i)
        fopen = open(fopenname, 'r')
        for line in fopen:
            if line:
                fpile.write(line)
        fopen.close()
        delete_com = options.outputprefix + '.temporary.' + str(i)
        delete_pile = options.outputprefix + '.pileup.' + str(i)
        if os.path.exists(delete_com):
            os.remove(delete_com)
        if os.path.exists(delete_pile):
            os.remove(delete_pile)
    fpile.close()

#   delete header
    headname2 = options.outputprefix + '.header'
    if os.path.exists(headname2):
        os.remove(headname2)

#   program of extracting damage information
if options.extract:
# define the cigar_mode
    cigar_mode = re.compile("([0-9]*)([A-Z])")

    def damageinfo(inread, inmd):
        read_ref = ''
        readlen = len(inread)
        pos_five = 0
        c2t = 0
        g2a = 0
        subinfo = ''
        mdlist = re.findall('(\d+|\D+)', inmd)
        for e in mdlist:
            if e.isdigit():
                read_ref += inread[pos_five:pos_five + int(e)]
                pos_five += int(e)

            elif '^'in e:
                continue

            elif e.isalpha():
                read_ref += e
                ncount_other_five[pos_five] += 1
                ncount_other_three[readlen - pos_five - 1] += 1

                if e == 'C' and inread[pos_five] == 'T':
                    c2t += 1
                    ncount_c2t_five[pos_five] += 1
                    ncount_c2t_three[readlen - pos_five - 1] += 1
                    pos_five += 1
                    pos_three = readlen - pos_five + 1
                    subinfo = subinfo + '|' + str(pos_five) + 'T' + str(pos_three)

                elif e == 'G' and inread[pos_five] == 'A':
                    g2a += 1
                    ncount_g2a_five[pos_five] += 1
                    ncount_g2a_three[readlen - pos_five - 1] += 1
                    pos_five += 1
                    pos_three = readlen - pos_five + 1
                    subinfo = subinfo + '|' + str(pos_five) + 'A' + str(pos_three)
                else:
                    pos_five += 1
        read_reflen = len(read_ref)
        for i in range(read_reflen):
            n_count[i] += 1
            if read_ref[i] == 'C':
                ncount_c_five[i] += 1
                ncount_c_three[read_reflen - i - 1] += 1

            elif read_ref[i] == 'G':
                ncount_g_five[i] += 1
                ncount_g_three[read_reflen - i - 1] += 1

        subinfo = repr(c2t) + 'T' + repr(g2a) + 'A' + subinfo
        return subinfo

    f2name = wdpath + '/' + options.outputprefix + '.damageinfo'
    f2 = open(f2name, 'w')

#   record the count from five end
    n_count = [0] * options.length

    ncount_c2t_five = [0] * options.length
    ncount_g2a_five = [0] * options.length
    ncount_c2t_three = [0] * options.length
    ncount_g2a_three = [0] * options.length

    ncount_c_five = [0] * options.length
    ncount_c_three = [0] * options.length
    ncount_g_five = [0] * options.length
    ncount_g_three = [0] * options.length

    ncount_other_five = [0] * options.length
    ncount_other_three = [0] * options.length

#   read each line and process

    for line in sys.stdin:
        if '@' in line[0]:
            continue

        col = line.split('\t')
        read = col[9]
        read_325 = read[::-1]
        cigar = col[5]

        try:
            MD = line.split('MD:Z:')[1].split()[0].rstrip('\n')
        except IndexError:
            print ("This file has sequence without MD filed")
            continue

# if I and S exist, create new read for damage information or reads directly
        if 'I' in cigar or 'S' in cigar:
            newread = ''
            cigarcount = 0
            cigarcomponents = cigar_mode.findall(cigar)
            for p in cigarcomponents:
                if 'I' in p[1]:
                    cigarcount += int(p[0])
                elif 'S' in p[1]:
                    cigarcount += int(p[0])
                else:
                    cigarcount2 = cigarcount + int(p[0])
                    newread += read[cigarcount:cigarcount2]
                    cigarcount = cigarcount2
            damagetag = damageinfo(newread, MD)

        else:
            damagetag = damageinfo(read, MD)

    f2.write("total_count" + '\t' + "c_5end" + '\t' + "g_5end" + '\t' + "other_5end" + '\t' + "c2t_5end" + '\t' + "g2a_5end" + '\t' +
            "c_3end" + '\t' + "g_3end" + '\t' + "other_3end" + '\t' + "c2t_3end" + '\t' + "g2a_3end" + '\n')

    for j in range(len(n_count)):
        f2.write(str(n_count[j]) + '\t' + str(ncount_c_five[j]) + '\t' + str(ncount_g_five[j]) + '\t' + str(ncount_other_five[j]) + '\t' + str(ncount_c2t_five[j]) + '\t' + str(ncount_g2a_five[j]) + '\t' +
                str(ncount_c_three[j]) + '\t' + str(ncount_g_three[j]) + '\t' + str(ncount_other_three[j]) + '\t' + str(ncount_c2t_three[j]) + '\t' + str(ncount_g2a_three[j]) + '\n')
    f2.close()

if options.snpcalling:
    quality = options.quality
    if options.allsite:
        quality = 0.0

    if options.ploidy == 2:
        fpileupname = wdpath + '/' + options.pileupfile
        fpileup = open(fpileupname, 'r')
        nline = 0
        for line in fpileup:
            if line[0] != '@':
                nline += 1
        fpileup.close()

        line_each = int(nline/options.thread) + 1

        fpileup = open(fpileupname, 'r')

        filenum = 1
        filesplitname = wdpath + '/'+ options.outputprefix + '.temporary.pileup.' + str(filenum)
        filesplit = open(filesplitname, 'w')
        linecount = 0
        for line in fpileup:
            if line[0] == '@':
                continue
            filesplit.write(line)
            linecount += 1
            if linecount == line_each:
                filesplit.close()
                filenum += 1
                filesplitname = wdpath + '/'+ options.outputprefix + '.temporary.pileup.' + str(filenum)
                filesplit = open(filesplitname, 'w')
                linecount = 0
        filesplit.close()
        fpileup.close()

        if options.allsite:
            def split2(para):
                os.system('python AntCaller_snpcalling.py -f ' + options.outputprefix + '.temporary.pileup.' + str(para) + ' -d ' + options.damageinfo + ' -o ' + options.outputprefix + str(para) + ' -q ' + str(quality) + ' -a')
        else:
            def split2(para):
                os.system('python AntCaller_snpcalling.py -f ' + options.outputprefix + '.temporary.pileup.' + str(para) + ' -d ' + options.damageinfo + ' -o ' + options.outputprefix + str(para) + ' -q ' + str(quality))

        pool = ThreadPool(options.thread)
        pool.map(split2, range(1, options.thread + 1))
        pool.close()
        pool.join()

        fvcfname = wdpath + '/'+ options.outputprefix + '.vcf'
        fvcf = open(fvcfname, 'w')
        fvcf.write('#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + options.outputprefix + '\n')
        for i in range(1, options.thread + 1):
            f2splitname = wdpath + '/'+ options.outputprefix + str(i) + '.vcf'
            f2split = open(f2splitname, 'r')
            for line in f2split:
                fvcf.write(line)
            f2split.close()
            filesplitname = options.outputprefix + '.temporary.pileup.' + str(i)
            os.remove(f2splitname)
            os.remove(filesplitname)
        fvcf.close()

    if options.ploidy == 1:
        rr = 0.001
        homo = (1 - 0.001)/2.0
        FORMAT = 'GT:AD:DP:GQ:PL'

#   calculate the substitution rate
        fname = wdpath + '/' + options.damageinfo
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

#   function of calculating probability of 10 genotype
        def prob_cal(p_in):
            return sum(p_GD[:p_in]) + sum(p_GD[p_in+1:])

#   function of which
        def which(a, b):
            for j in range(len(a)):
                if a[j] == b:
                    return j

        f1name = wdpath + '/' + options.pileupfile  # read the pileup file
        f1 = open(f1name, 'r')
        f2name = wdpath + '/' + options.outputprefix + '.vcf'
        f2 = open(f2name, 'w')

        f2.write('#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + options.outputprefix + '\n')

#   process each line of sam file
        for line in f1:
            if '@' in line[0]:
                continue

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

            genotype = ['A', 'G', 'C', 'T']

            A_count = 0
            G_count = 0
            C_count = 0
            T_count = 0

            p_A = 1
            p_G = 1
            p_C = 1
            p_T = 1

            for n in range(count):
                r = error(quals[n]) + 0.01 + 1.0/pow(10.0, int(Mquals[n])/10.0)
                q = g2a_rate[int(pos3[n]) - 1]
                p = c2t_rate[int(pos5[n]) - 1]
                if Nt[n] == 'A':
                    A_count += 1
                    p_A *= 1.0 - r
                    p_G *= (1.0 - r)*q + r*(1.0-q)/3.0
                    p_C *= r/3.0
                    p_T *= r/3.0
                elif Nt[n] == 'G':
                    G_count += 1
                    p_A *= r/3.0
                    p_G *= (1.0 - q)*(1.0 - r) + r*q/3.0
                    p_C *= r/3.0
                    p_T *= r/3.0
                elif Nt[n] == 'C':
                    C_count += 1
                    p_A *= r/3.0
                    p_G *= r/3.0
                    p_C *= (1.0 - p)*(1.0 - r) + r*p/3.0
                    p_T *= r/3.0
                elif Nt[n] == 'T':
                    T_count += 1
                    p_A *= r/3.0
                    p_G *= r/3.0
                    p_C *= p*(1.0 - r) + r*(1.0 - p)/3.0
                    p_T *= 1.0 - r
                else:
                    continue

            p_DG = np.array([p_A, p_G, p_C, p_T])
            p_GD = p_DG/sum(p_DG)

            if 0 in p_GD:
                for i in range(len(p_GD)):
                    if p_GD[i] == 0:
                        p_GD[i] = 1.000000000000017e-308

            max_pos = p_GD.argmax()
            alt = genotype[max_pos]
            count_list = np.array([A_count, G_count, C_count, T_count])
            if options.allsite:
                if ref == alt:
                    ALT = '.'
                    Qual = round(-10.0 * math.log(prob_cal(max_pos), 10.0), 2)
                    INFO = 'AC=0;AF=0.00;AN=1;DP=' + str(count)
                    PL = int(round(-10.0*math.log(p_GD[max_pos], 10.0)))
                    Sample1 = '0:' + str(count_list[max_pos]) + ',' + str(sum(count_list)-count_list[max_pos]) + ':' + str(count) + ':' +str(Qual)+ ':' + str(PL)
                    f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
                else:
                    ALT = alt
                    Qual = round(-10.0 * math.log(p_GD[which(genotype, ref)], 10.0), 2)
                    INFO = 'AC=1;AF=1.00;AN=1;DP=' + str(count)
                    PL = str(int(round(-10.0*math.log(p_GD[which(genotype, ref)], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[which(genotype, alt)], 10.0))))
                    Sample1 = '1:' + str(count_list[which(genotype, ref)]) + ',' + str(sum(count_list)-count_list[which(genotype, ref)]) + ':' + str(count) + ':' +str(Qual)+ ':' + str(PL)
                    f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')
            elif ref != alt:
                Qual = round(-10.0 * math.log(p_GD[which(genotype, ref)], 10.0), 2)
                if Qual < quality:
                    continue
                ALT = alt
                Qual = round(-10.0 * math.log(p_GD[which(genotype, ref)], 10.0), 2)
                INFO = 'AC=1;AF=1.00;AN=1;DP=' + str(count)
                PL = str(int(round(-10.0*math.log(p_GD[which(genotype, ref)], 10.0)))) + ',' + str(int(round(-10.0*math.log(p_GD[which(genotype, alt)], 10.0))))
                Sample1 = '1:' + str(count_list[which(genotype, ref)]) + ',' + str(sum(count_list)-count_list[which(genotype, ref)]) + ':' + str(count) + ':' +str(Qual)+ ':' + str(PL)
                f2.write(chromosome + '\t' + position + '\t' + '.' + '\t' + ref + '\t' + ALT + '\t' + str(Qual) + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + Sample1 + '\n')

        f2.close()


time_end = time.time()
print("running time of AntCaller:" + str(time_end - time_start))