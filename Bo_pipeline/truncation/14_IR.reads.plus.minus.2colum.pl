#!/usr/bin/perl -w
## perl IR.reads.plus.minus.2colum.pl final.summary.txt
open A,"$ARGV[0]"||die $!;
print "Ref	Start_position	Number_of_total_reads	Percentage_of_truncation	IR_read(+)	IR_read(-)	Truncation_or_not	Read_start(+)	Read_start(-)	Read_end(+)	Read_end(-)\n";
while(<A>){
    chomp;
    @aa=split(/\t/,$_);
    $aa[1]=$aa[1]+1;
    if($aa[5]=~/\+/){
       print "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t$aa[6]\t$aa[8]\t$aa[9]\t$aa[10]\t$aa[13]\t$aa[12]\t$aa[11]\n";}
   elsif($aa[5]=~/\-/){
       print "$aa[0]\t$aa[1]\t$aa[3]\t$aa[4]\t$aa[8]\t$aa[6]\t$aa[9]\t$aa[10]\t$aa[13]\t$aa[12]\t$aa[11]\n";}
}

close A;

