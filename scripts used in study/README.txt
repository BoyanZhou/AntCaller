##Simulation

"Simulation" folder contains scripts and files used for the simulation study in article.

Python script "make-fq-damage.py" was used to generate damage pattern for fastq files.

#Diploid genome simulation

1.The original vcf file of chr21 from the altai Neandertal was downloaded from http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/

2.InDels were removed from the original vcf file using python script "indel-remove.py"

3.The fasta file of chr21(GRCh37) was modified as two chr21 of the
Neandertal according to the InDels-removed vcf file using python script "chr21-modify.py"

#Haploid genome simulation

1.Tianyuan_KC417443.fas: haploid genome reference

2.The process of generating data and calling SNPs was detailed in bash script "haploid genome simulation.sh"

3.Python script "snp_count.py" was used to count called SNPs


##Real Data

1.Python script "autosome-select": Filter autosomal variants from called variants

2.Python script "compare-vcf.py": Compare called variants with the database

