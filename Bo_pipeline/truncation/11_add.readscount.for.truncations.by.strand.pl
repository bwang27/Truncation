#!/usr/bin/perl -w
##  perl add.readscount.for.truncations.by.strand.pl  ../../../../bedgraph_by_strand_0_based/ssAAV-SCA3-3.bedgraph ../uniq.ssAAV-SCA3-3.truncation.vs.total.reads.plusandminus >uniq.ssAAV-SCA3-3.truncation.vs.total.reads.plusandminus.refine
open A,"$ARGV[0]"||die $!;
$bed=join'',<A>;
@bedgraph=split(/\n/,$bed);
close A;

open B,"$ARGV[1]"||die $!;
while(<B>){
    chomp;
    print "$_\t";
     @aa=split(/\t/,$_);
      foreach $line(@bedgraph){
          @bb=split(/\t/,$line);
         #if(($aa[1]==$bb[1])&&($aa[5]=~/\$bb[4]/)){
         if($aa[1]==$bb[1]){
      print "$bb[3]\t$bb[4]\t";}
     }
 
   print "truncation\n";}

close B;

