#!/bin/sh
for((i=1;i<=1000;i++))
do
art_illumina -i /home/Lab/homo/Neandertal/vindija/simulationMT/Tianyuan_KC417443.fas -l 60 -f 10 -ss HS20 -o simu_Tianyuan_MT
python make-fq-damage-3.0.py -f /home/Lab/homo/Neandertal/vindija/simulationMT/repeat/5fold/5high/simu_Tianyuan_MT.fq -o simu_Tianyuan_MT.damage.fq -r 0.40
bwa aln -n 0.02 -o 2 -l 16500 /home/Lab/download/pipeline/GATK_bundle/MT/modern/forGATK/rCRS-NC_012920.1.fa simu_Tianyuan_MT.damage.fq > simu_Tianyuan_MT.damage.sai
bwa samse /home/Lab/download/pipeline/GATK_bundle/MT/modern/forGATK/rCRS-NC_012920.1.fa simu_Tianyuan_MT.damage.sai simu_Tianyuan_MT.damage.fq > simu_Tianyuan_MT.damage.sam
samtools view -Sb simu_Tianyuan_MT.damage.sam > simu_Tianyuan_MT.damage.bam
samtools sort simu_Tianyuan_MT.damage.bam simu_Tianyuan_MT.damage.sort
samtools rmdup -S simu_Tianyuan_MT.damage.sort.bam simu_Tianyuan_MT.damage.sort.rmdup.bam
samtools index simu_Tianyuan_MT.damage.sort.rmdup.bam

samtools mpileup -q 20 -Q 0 -u -I -f /home/Lab/download/pipeline/GATK_bundle/MT/modern/forGATK/rCRS-NC_012920.1.fa simu_Tianyuan_MT.damage.sort.rmdup.bam > simu_Tianyuan_MT.damage.sort.rmdup.samtools.pileup
bcftools view -s sample.txt -vcg simu_Tianyuan_MT.damage.sort.rmdup.samtools.pileup > simu_Tianyuan_MT.damage.sort.rmdup.samtools.vcf

java -Xmx10g -jar /home/Lab/src/picard-tools-1.137/picard.jar AddOrReplaceReadGroups \
    INPUT=simu_Tianyuan_MT.damage.sort.rmdup.bam \
    OUTPUT=simu_Tianyuan_MT.damage.sort.rmdup.group.bam \
    RGID=group RGLB=lib RGPL=illumina RGPU=unit RGSM=test

samtools index simu_Tianyuan_MT.damage.sort.rmdup.group.bam

java -Xmx10g -jar /home/Lab/src/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R /home/Lab/download/pipeline/GATK_bundle/MT/modern/forGATK/rCRS-NC_012920.1.fa \
   -ploidy 1 \
   -stand_call_conf 0 \
   -I simu_Tianyuan_MT.damage.sort.rmdup.group.bam \
   -o simu_Tianyuan_MT.damage.sort.rmdup.group.gatk.vcf
   
samtools view -h simu_Tianyuan_MT.damage.sort.rmdup.bam | python AntCaller-1.1.py --pileup -o simu_Tianyuan_MT.ant
samtools view -h simu_Tianyuan_MT.damage.sort.rmdup.bam | python AntCaller-1.1.py --extract -o simu_Tianyuan_MT.ant
python AntCaller-1.1.py --snpcalling -o simu_Tianyuan_MT.ant -d simu_Tianyuan_MT.ant.damageinfo -f simu_Tianyuan_MT.ant.AntCaller.pileup -p 1 -q 0

python snp_compare.py -n 16
done
