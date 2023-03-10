#!/usr/bin/perl -w

open A,"$ARGV[0]"||die $!;
#open B,"$ARGV[1]"||die $!;
#open O,">ssAAV-SCA3-3.truncation.point.bed"||die $!;
#$bed=join'',<B>;
#@bedgraph=split(/\n/,$bed);
while(<A>){
    chomp;
    @aa=split(/\s+/,$_);
    if($aa[2]==$aa[5]){
     $end=$aa[2]+1;
 #     foreach $line(@bedgraph){
  #      @bb=split(/\t/,$line);
  #       if(($aa[2]==$bb[1])&&($bb[2]==$end)){
  #  print O  "$line\n;"}}}
#          print "$aa[0]\t$aa[2]\t$end\n";}
    print "ID23_scAAV-SCA3-4E10x2\t$aa[2]\t$end\n";}
 
   else{
       #$truncation=$aa[2]+int(abs($aa[5]-$aa[2])/2);
       if($aa[2]>$aa[5]){
       $truncation=$aa[2];
       $end=$truncation+1;}
       if($aa[2]<$aa[5]){
       $truncation=$aa[5];
       $end=$truncation+1;}
   
  #   foreach $line(@bedgraph){
  #      @bb=split(/\t/,$line);
  #       if(($truncation==$bb[1])&&($bb[2]==$end)){
  #  print O  "$line\n;"}}}

#       print "$aa[0]\t$truncation\t$end\n";}
    print "ID23_scAAV-SCA3-4E10x2\t$truncation\t$end\n";}
}

close A;

