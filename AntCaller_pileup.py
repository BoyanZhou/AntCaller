#!/usr/bin/env python
import sys

# define the cigar_mode
import re
cigar_mode = re.compile("([0-9]*)([A-Z])")

def cut(x):
    return str(x).lstrip("[").rstrip("]")

#   define the function of pileup a read
def pileup(inread, inref, inquals, inMquals, inpos, inend5, inend3):
    for e in range(len(inread)):
        try:
            if inpos[e] in poslist:
                location = poslist.index(inpos[e])
                countlist[location] += 1
                Ntlist[location] += inread[e]
                qualist[location] += inquals[e]
                Mqualist[location].append(inMquals)
                end5list[location].append(inend5[e])
                end3list[location].append(inend3[e])
            else:
                poslist.append(inpos[e])
                reflist.append(ref[e])
                countlist.append(1)
                Ntlist.append(inread[e])
                qualist.append(inquals[e])
                Mqualist.append([inMquals])
                end5list.append([inend5[e]])
                end3list.append([inend3[e]])
        except IndexError:
            print line
            print inread
            print inref
            print inquals
            print inMquals
            print inpos
            print inend5
            print inend3
            exit()
# end the function of pileup a read


poslist = []
reflist = []
countlist = []
Ntlist = []
qualist = []
Mqualist = []
end5list = []
end3list = []
chromosome = []

#   read each line and process
for line in sys.stdin:
    if line[0] == '@':
        print line.rstrip('\n')
        continue

    col = line.split('\t')
    readname = col[0]
    flag = col[1]
    if col[2] != chromosome and chromosome:
        while poslist:
            print (chromosome + '\t' + str(poslist.pop(0)) + '\t' + reflist.pop(0) + '\t' + str(countlist.pop(0)) + '\t' + Ntlist.pop(0) + '\t' + qualist.pop(0) + '\t' + cut(end5list.pop(0)) + '\t' + cut(end3list.pop(0)) + '\t' + cut(Mqualist.pop(0)))
    chromosome = col[2]
    position = int(col[3])
    MAPQ = int(col[4])
    if MAPQ < 20:
        continue

    cigar = col[5]
    read = col[9]
    quals = col[10]
    noMD = 0

    try:
        MD = line.split('MD:Z:')[1].split()[0].rstrip('\n')
    except IndexError:
        noMD += 1
        continue
    mdlist = re.findall('(\d+|\D+)', MD)

#   determine if the read is reverse
    reverse = False
    if flag.isdigit():
        if int(flag) & 16:
            reverse = True
    else:
        if 'r' in flag:
            reverse = True
#   end determine if the read is reverse


#   output the pileup
    while poslist and poslist[0] < position:
        print (chromosome + '\t' + str(poslist.pop(0)) + '\t' + reflist.pop(0) + '\t' + str(countlist.pop(0)) + '\t' + Ntlist.pop(0) + '\t' + qualist.pop(0) + '\t' + cut(end5list.pop(0)) + '\t' + cut(end3list.pop(0)) + '\t' + cut(Mqualist.pop(0)))
#   end output the pileup


#   create the pileup input
    if 'I' in cigar or 'S' in cigar or 'H' in cigar:
        end5 = []
        end3 = []
        newcigar = ''
        newread = ''
        newquals = ''
        cigarcount = 0
        cigarcomponents = cigar_mode.findall(cigar)
        for p in cigarcomponents:
            if 'I' in p[1] or 'S' in p[1] or 'H' in p[1]:
                cigarcount += int(p[0])

            elif 'D' in p[1]:
                continue
            else:
                newcigar += ''.join(list(p))
                for i in range(int(p[0])):
                    end5.append(cigarcount + i + 1)
                    end3.append(len(read) - cigarcount - i)

                cigarcount2 = cigarcount + int(p[0])
                newread += read[cigarcount:cigarcount2]
                newquals += quals[cigarcount:cigarcount2]
                cigarcount = cigarcount2

    else:
        newread = read[:]
        newquals = quals[:]
        newcigar = cigar[:]
        end5 = []
        end3 = []
        for n in range(len(newread)):
            end5.append(n + 1)
            end3.append(len(read) - n)

#   create ref and position in reference
    pos_inref = []
    ref = ''
    mdcount = 0
    refcount = position + 0
    for m in mdlist:
        if m.isdigit():
            pos_inref.extend(range(refcount, refcount + int(m)))
            refcount += int(m)
            ref += newread[mdcount: mdcount + int(m)]
            mdcount += int(m)

        elif '^'in m:
            refcount += len(m) -1
            continue

        elif m.isalpha():
            ref += m
            mdcount += 1
            pos_inref.append(refcount)
            refcount += 1

    if reverse:
        newread = newread.lower()
        end_midlle = end3[:]
        end3 = end5[:]
        end5 = end_midlle[:]
    else:
        newread = newread.upper()

    pileup(newread, ref, newquals, MAPQ, pos_inref, end5, end3)

while poslist:
    print (chromosome + '\t' + str(poslist.pop(0)) + '\t' + reflist.pop(0) + '\t' + str(countlist.pop(0)) + '\t'+ Ntlist.pop(0) + '\t' + qualist.pop(0) + '\t'+ cut(end5list.pop(0)) + '\t' + cut(end3list.pop(0)) + '\t' + cut(Mqualist.pop(0)))