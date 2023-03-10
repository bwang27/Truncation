#!/usr/bin/perl -w
##perl get.ref.start.end.counts.pl ssAAV-SCA3-3.bed sorted.merged.truncation.normal.reads.coverage
#
open A,"$ARGV[0]"||die $!;
$all=join'',<A>;
@reads=split(/\n/,$all);
close A;

open B,"$ARGV[1]"||die $!;
open O,">plus.sorted.merged.truncation.normal.reads.coverage.withref.startend"||die $!;
while(<B>){
    chomp;
    @aa=split(/\t/,$_);
    $start_plus=0;
    $start_minus=0;
    $end_plus=0;
    $end_minus=0;
        foreach $read(@reads){
            @bb=split(/\t/,$read);
            if(($aa[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($bb[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($aa[1]==$bb[1])&&($bb[5]=~/\+/)){
         $start_plus++;}
            if(($aa[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($bb[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($aa[1]==$bb[1])&&($bb[5]=~/\-/)){
         $start_minus++;}
            if(($aa[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($bb[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($aa[1]==$bb[2])&&($bb[5]=~/\+/)){
         $end_plus++;} 
            if(($aa[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($bb[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($aa[1]==$bb[2])&&($bb[5]=~/\-/)){
         $end_minus++;}

      }
  
print O "$_\t$start_plus\t$start_minus\t$end_plus\t$end_minus\n";}



close B;
close O;
