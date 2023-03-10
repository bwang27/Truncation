#!/usr/bin/perl -w
## perl revised_09_merge.plus.minus.table.pl ID23_scAAV-SCA3-4E10x2_AAV.truncation.vs.total.reads.sense ID23_scAAV-SCA3-4E10x2_AAV.truncation.vs.total.reads.antisense 
open A,"$ARGV[0]"|| die $!;
open B,"$ARGV[1]"|| die $!;
open AO,">$ARGV[0].plus"||die $!;
open BO,">$ARGV[1].mius"||die $!;
while(<A>){
     chomp;
     @aa=split(/\t/,$_);
     $percent=($aa[3]/$aa[4])*100;
     $percentage= sprintf("%.3f",$percent);
 #    print "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t\+\n";}
     print AO "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t$percentage%\t\+\n";}
while(<B>){
     chomp;
     @aa=split(/\t/,$_);
     $percent=($aa[3]/$aa[4])*100;
     $percentage= sprintf("%.3f",$percent);
 #    print "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t\+\n";}
     print BO "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t$percentage%\t\-\n";}
 #
close A;
close B;
close AO;
close BO;

