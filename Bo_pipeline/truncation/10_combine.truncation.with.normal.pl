#!/usr/bin/perl -w
## perl combine.truncation.with.normal.pl uniq.ssAAV-SCA3-3.truncation.vs.total.reads.plusandminus ../../../ssAAV-SCA3-3.sorted.bedgraph
open A, "$ARGV[0]"||die $!;
$truncation=join'',<A>;
@truncation_reads=split(/\n/,$truncation);
close A;

open B,"$ARGV[1]"||die $!;
open O,">newtmp.coverage"||die $!;
#open O,">normal.reads.without.truncations.coverage"||die $!;
while(<B>){
     chomp;
     $num=scalar@truncation_reads;
     @aa=split(/\t/,$_);
     $count=0;
     foreach $line(@truncation_reads){
          @bb=split(/\t/,$line);
        #  if(($aa[1]==$bb[1])&&($aa[0]=~/ssAAV-SCA3-3_shifted/)&&($bb[0]=~/ssAAV-SCA3-3_shifted/)){
#   	      next;
     #         print "$_\n";}
    # last;
     if(($aa[1]!=$bb[1])&&($aa[0]=~/ID23_scAAV-SCA3-4E10x2/)&&($bb[0]=~/ID23_scAAV-SCA3-4E10x2/)){
           $count++;
        if($count==$num){
         print O "$aa[0]\t$aa[1]\t0\t$aa[3]\t0\t\.\tnormal\n";} 
      }
    } 
# { print O "$aa[0]\t$aa[1]\t0\t$aa[3]\t0\t\.\n";}
}
#}
#}

system "cat newtmp.coverage | uniq >newnormal.reads.without.truncations.coverage";
