#!/usr/bin/perl -w
## perl change.bedgraph.to.continous.position.pl ssAAV-SCA3-3.sorted.bedgraph | grep 'ssAAV-SCA3-3_shifted' | wc -l
open A,"$ARGV[0]"||die $!;
while(<A>){
   chomp;
  if($_=~/ID23_scAAV-SCA3-4E10x2/){
   @aa=split(/\t/,$_);
   if($aa[2]-$aa[1]==1){
      $start=$aa[1];
      $end=$aa[2];
      print "$aa[0]\t$start\t$end\t$aa[3]\n";}
   elsif($aa[2]-$aa[1]>1){
       $minus=$aa[2]-$aa[1];
       for($i=0;$i<$minus;$i++){
           $start=$aa[1]+$i;
           $end=$start;
           print "$aa[0]\t$start\t$end\t$aa[3]\n";}
   }
}

}

close A;

