#!/usr/bin/perl -w
##usage: perl step4_get_plus_minus_mapping_reads.pl  ssAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments >plus.minus.step4.ssAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments
open A,"$ARGV[0]"||die $!;
while(<A>){
   chomp;
   if(($_=~/\t\+\t/)&&($_=~/\t\-\t/)){
   print "$_\n";}
}

close A;


