#!/usr/bin/perl -w
##perl get.truncation.position.pl uniq.plus.minus.step4.ssAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments.sense.truncation.point_step2 ../ssAAV-SCA3-3.sorted.bedgraph
#
open A,"$ARGV[0]"||die $!;
open B,"$ARGV[1]"||die $!;
open O,">$ARGV[2]"||die $!; #ssAAV-SCA3-3.truncation.point.bed"||die $!;
$bed=join'',<B>;
@bedgraph=split(/\n/,$bed);
while(<A>){
    chomp;
    @aa=split(/\s+/,$_);
    if($aa[1]==$aa[4]){
     $end=$aa[1]+1;
      foreach $line(@bedgraph){
        if($line=~/ID23_scAAV-SCA3-4E10x2/){
        @bb=split(/\t/,$line);
        if(($aa[1]==$bb[1])&&($bb[2]==$end)){
    print O "$line\t$_\n";}}}}
#          print "$aa[0]\t$aa[2]\t$end\n";}
  #  print "ssAAV-SCA3-3_shifted\t$aa[2]\t$end\n";}
 
   else{
        if($aa[1]<$aa[4]){
        $truncation=$aa[1]; #$aa[1]+int(abs($aa[4]-$aa[1])/2);
        $end=$truncation+1;}
        if($aa[1]>$aa[4]){
        $truncation=$aa[4]; #$aa[4]+int(abs($aa[4]-$aa[1])/2);   
        $end=$truncation+1;}
     foreach $line(@bedgraph){
        if($line=~/ID23_scAAV-SCA3-4E10x2/){
        @bb=split(/\t/,$line);
        if(($truncation==$bb[1])&&($bb[2]==$end)){
    print O  "$line\t$_\n";}}}}

#       print "$aa[0]\t$truncation\t$end\n";}
  #  print "ssAAV-SCA3-3_shifted\t$truncation\t$end\n";}
}

close A;

