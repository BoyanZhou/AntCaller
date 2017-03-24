#!/usr/bin/env python
fvcf = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/AltaiNea.hg19_1000g.21.mod.vcf','r')
f2 = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/AltaiNea.hg19_1000g.21.indel_removed.vcf','w')
f3 = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/AltaiNea.hg19_1000g.21.indel.vcf','w')
indel = 0
for line in fvcf:
    if line[0] == '#':
        f2.write(line)
        continue

    col = line.split('\t')
    chromosome = col[0]
    pos = int(col[1])
    ref = col[3]
    alt = col[4]
    if len(ref) != 1 or len(alt) != 1:
        if not (',' in alt):
            indel += 1
            f3.write(line)
            continue

    f2.write(line)

print(indel)
fvcf.close()
f2.close()
