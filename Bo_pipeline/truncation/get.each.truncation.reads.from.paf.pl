#!/usr/bin/perl -w
##usage: perl get.each.truncation.reads.from..paf.pl ssAAV-SCA3-3.coordinates 1073.ssAAV-SCA3-3.uniq.sense.reads
## perl get.each.truncation.reads.from.paf.pl ID23_scAAV-SCA3-4E10x2_AAV.reads.coordinates uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.sense.truncation_point_step2 uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.antisense.truncation_point_step2
#
open A,"$ARGV[0]"||die $!;
$all=join'',<A>;
@paf=split(/\n/,$all);
close A;

open B,"$ARGV[1]"||die $!;
open SENSE,">IR.sense.reads.position"||die $!;
while(<B>){
    chomp;
    @bb=split(/\t/,$_);
   print SENSE "$bb[0]\t";
    foreach $line(@paf){
       @aa=split(/\t/,$line);
       if($aa[0] eq $bb[0]){
     print SENSE "$aa[2]/$aa[7]\t$aa[3]/$aa[8]\t$aa[4]\t";}
}
     print SENSE "\n";
}

close B;
close SENSE;

open C,"$ARGV[2]"||die $!;
open ANTISENSE,">IR.antisense.reads.position"||die $!;
while(<C>){
    chomp;
    @bb=split(/\t/,$_);
   print ANTISENSE "$bb[0]\t";
    foreach $line(@paf){
       @aa=split(/\t/,$line);
       if($aa[0] eq $bb[0]){
     print ANTISENSE "$aa[2]/$aa[7]\t$aa[3]/$aa[8]\t$aa[4]\t";}
}
     print ANTISENSE "\n";
}

close C;
close ANTISENSE;

