#!/usr/bin/perl -w
## perl get.reads.position.supported.IR.pl uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.sense.truncation_point_step2  final.uniq.IR.sense.reads.position.10.together.continuous  uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.antisense.truncation_point_step2  final.uniq.IR.antisense.reads.position.10.together.continuous
open A,"$ARGV[0]"||die $!;
$IR=join'',<A>;
@irreads=split(/\n/,$IR);
close A;

open B,"$ARGV[1]"||die $!;
open BO,">sense.IR.supported.by.reads.position"||die $!;
open BOO,">sense.IR.supported.by.reads.position.for.downstream.analysis"||die $!;
while(<B>){
     chomp;
     @aa=split(/\t/,$_);
      foreach $line(@irreads){
          @bb=split(/\t/,$line);
      if($bb[0] eq $aa[0]){
    print BO "$_\t$line\t";
    print BOO "$line\n";}
}

print BO "\n";}

close B;
close BO;
close DOO;

open C,"$ARGV[2]"||die $!;
$IRanti=join'',<C>;
@irreadsanti=split(/\n/,$IRanti);
close C;

open D,"$ARGV[3]"||die $!;
open DO,">antisense.IR.supported.by.reads.position"||die $!;
open DOO,">antisense.IR.supported.by.reads.position.for.downstream.analysis"||die $!;
while(<D>){
     chomp;
     @aa=split(/\t/,$_);
      foreach $line(@irreadsanti){
          @bb=split(/\t/,$line);
      if($bb[0] eq $aa[0]){
    print DO "$_\t$line\t";
    print DOO "$line\n";}
}

print DO "\n";}

close D;
close DO;
close DOO;


