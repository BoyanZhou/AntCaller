# AntCaller
SNP calling and genotyping for damaged ancient DNA

# Dependency
1. Python 2.7
2. samtools

# Usage
1. This program is built on Python2.7 and can only be used for SAM/BAM files with MD field under Linux system.
2. Copy three .py files to the directory where SAM/BAM files are stored.
3. Print "python AntCaller-1.1.py -help" for detailed usage information.
4. SNP calling can only be performed after extracting damage information and pileuping reads.

# Option
1. Default length of the longest reads in BAM file is set to 300. If there are reads that are longer than 300, "-l" option should be used in extracting damage information. For example, you can use "samtools view -h example.bam | python AntCaller-1.1.py --extract -l 600 -o example", if the longest reads are longer than 300bp but shorter than 600bp.

# Examples
1. for BAM file (Require SAMtools installed)

  1) pileup reads, output example.AntCaller.pileup
  
  samtools view -h example.bam | python AntCaller-1.1.py --pileup -o example
  
  2) extract damage information, output example.damageinfo
  
  samtools view -h example.bam | python AntCaller-1.1.py --extract -o example
  
  3) SNP calling, output vcf file
  
  python AntCaller-1.1.py --snpcalling -o example -d example.damageinfo -f example.AntCaller.pileup

2. for SAM file

  1) pileup reads, output example.AntCaller.pileup
  
  cat example.sam | python AntCaller-1.1.py --pileup -o example
  
  2) extract damage information, output example.damageinfo
  
  cat example.sam | python AntCaller-1.1.py --extract -o example
  
  3) SNP calling, output vcf file
  
  python AntCaller-1.1.py --snpcalling -o example -d example.damageinfo -f example.AntCaller.pileup

# Contact
For any problems, please contact boyanzhou1992@gmail.com
