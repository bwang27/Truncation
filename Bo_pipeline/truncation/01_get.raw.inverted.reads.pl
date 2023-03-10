#!/usr/bin/perl -w
##perl get.inverted.reads.pl ssAAV-SCA3-3.step1.same.genome ssAAV-SCA3-3.step2.candidate.reads.ids
#
open A,"$ARGV[0]"||die $!;
$all=join'',<A>;
@bed=split(/\n/,$all);
close A;

open B,"$ARGV[1]"||die $!;
open O,">$ARGV[0].step3.sameread.all.alignments"||die $!;
while(<B>){
   chomp;
   print O "$_\t";
   foreach $line(@bed){
     @aa=split(/\t/,$line);
     if($aa[3] eq $_){
   print O "$aa[1]\t$aa[2]\t$aa[5]\t";}
}
   print O "\n";
}

close B;

