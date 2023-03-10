#!/usr/bin/perl -w
##perl merge.rtuncation.reads.and.total.reads.pl number_of_truncation_reads_per_point ssAAV-SCA3-3.truncation.point.bed;
#
open A,"$ARGV[0]"||die $!;
$count=join'',<A>;
@truncation_reads=split(/\n/,$count);
close A;

open B,"$ARGV[1]"||die $!;
open O,">$ARGV[2]"||die $!; #ssAAV-SCA3-3.truncation.vs.total.reads"||die $!;
while(<B>){
    chomp;
    @aa=split(/\t/,$_);
      foreach $line(@truncation_reads){
           @bb=split(/\s+/,$line);
  # print "$bb[1]\n";
           if($aa[1]==$bb[2]){
      print O "$aa[0]\t$aa[1]\t$aa[2]\t$bb[1]\t$aa[3]\t$aa[4]\t$aa[5]\t$aa[6]\n";}
     }
}

close B;

