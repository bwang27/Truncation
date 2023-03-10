1. Convert SAM to BAM for sorting
 samtools view -S -b my.sam > my.bam

2. Sort BAM for SNP calling
 samtools sort my.bam my-sorted

3. Index the genome assembly
 samtools faidx my.fasta

4. Run 'mpileup' to generate VCF format
 samtools mpileup -g -f my.fasta my-sorted1.bam my-sorted-2.bam my-sorted-n.bam > myraw.bcf

5. Call SNPs...
 bcftools view -bvcg my-raw.bcf > my-var.bcf

6. Filter SNPs
 bcftools view my.var.bcf | vcfutils.pl varFilter - > my.var-final.vcf

#/drive2/bwang/software/samtools-0.1.19/bcftools/bcftools view -bvcg my-raw_ID43.bcf >my-var_ID43.bcf
