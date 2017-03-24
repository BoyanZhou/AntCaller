f = open(r'/home/Lab/homo/Neandertal/vindija/human_g1k_v37_chr21.fasta','r')
fvcf = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/AltaiNea.hg19_1000g.21.indel_removed.vcf','r')
fout = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/human_g1k_v37_chr21.modify.fasta','w')
fout2 = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/human_g1k_v37_chr21.modify2.fasta','w')
indel = open(r'/home/Lab/homo/Neandertal/vindija/simulationchr21/indel.vcf','w')

unmatch = 0
pos_start = 1
line_vcf = '#'

for line in f:
    if line[0] == '>':
        fout.write(line)
        fout2.write(line)
        continue

    line = line.strip('\n')
    len_line = len(line)
    line1 = list(line)
    line2 = list(line)

    if not line_vcf:
        fout.write(''.join(line1) + '\n')
        fout2.write(''.join(line2) + '\n')
        pos_start += len_line
        continue

    while line_vcf[0] == '#':
        line_vcf = fvcf.readline()

    col = line_vcf.split('\t')
    pos = int(col[1])

    if pos >= (pos_start + len_line):
        fout.write(''.join(line1) + '\n')
        fout2.write(''.join(line2) + '\n')
        pos_start += len_line
        continue

    while pos < (pos_start + len_line):
        chromosome = col[0]
        ID = col[2]
        ref = col[3]
        alt = col[4]
        qual = col[5]
        geno = col[9].split(':')[0]
        geno1 = geno.split('/')[0]
        geno2 = geno.split('/')[1]
        if not line[pos - pos_start] == ref:
            unmatch += 1
            indel.write(line_vcf + '\n')

        if geno1 == '1':
            line1[pos - pos_start] = alt[0]
        if geno2 == '1':
            line2[pos - pos_start] = alt[0]

        line_vcf = fvcf.readline()
        if not line_vcf:
            break
        col = line_vcf.split('\t')
        pos = int(col[1])

    fout.write(''.join(line1) + '\n')
    fout2.write(''.join(line2) + '\n')
    pos_start += len_line

print(unmatch)
f.close()
fvcf.close()
fout.close()
fout2.close()
indel.close()


