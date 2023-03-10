#!/usr/bin/perl -w
#usage: perl judge.the.read.inverted.or.not.pl plus.minus.step4.scAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments;
#perl modify.judge.the.read.inverted.or.not.pl plus.minus.step4.ssAAV-SCA3-4.step1.same.genome.step3.sameread.all.alignments 10;
#
open A,"$ARGV[0]"||die $!;
open O,">$ARGV[0].chunk"||die $!;
while(<A>){
    chomp;
    chomp;
    @aa = split(/\t/,$_);
    $read_id = shift(@aa);
    print O "$read_id\t";
    my @ret; 
    my @new = array_chunk(3, @aa);
     foreach $line(@new){
        my @array = @$line;
          print O "@array\t";}
 print O "\n";
   }

#close A;

open B,"$ARGV[0].chunk"||die $!;
open SENSE,">$ARGV[0].sense"||die $!;
open ANTISENSE,">$ARGV[0].antisense"||die $!;
while(<B>){
   chomp;
   @nd = split(/\t/,$_);
   $read=$nd[0];
 #  $id = $nd[0];
   shift(@nd);
 #  $first = shift(@nd);
 #  @seperate_first = split(/\s/,$first);
     for ($i=0;$i<@nd;$i++){
        $qry=$nd[$i];
        @split_qry=split(/\s/,$qry);
        foreach $element(@nd){
        if($element ne $qry){
          #  @split_qry=split(/\s/,$qry);
        @seperate_element = split(/\s/,$element);
#        if(($seperate_first[2]!~/\$seperate_element[2]/)&&($seperate_first[1]==$seperate_element[1])){
        $left=$seperate_element[1]-"$ARGV[1]"; #30; #20;
        $right=$seperate_element[1]+"$ARGV[1]"; #30; #20;
        $left_start=$seperate_element[0]-"$ARGV[1]"; #30;
        $right_start=$seperate_element[0]+"$ARGV[1]"; #30; 
        if(($split_qry[2] ne $seperate_element[2])&&($split_qry[1]>$left)&&($split_qry[1]<$right)){
        print SENSE "$read\t$split_qry[0]\t$split_qry[1]\t$split_qry[2]\t$seperate_element[0]\t$seperate_element[1]\t$seperate_element[2]\n";}
        elsif(($split_qry[2] ne $seperate_element[2])&&($split_qry[0]>$left_start)&&($split_qry[0]<$right_start)){
        print ANTISENSE "$read\t$split_qry[0]\t$split_qry[1]\t$split_qry[2]\t$seperate_element[0]\t$seperate_element[1]\t$seperate_element[2]\n";}

          
  		}
 	}	
  }
}

close A;
close B;

system "cat $ARGV[0].sense | sort | uniq >uniq.$ARGV[0].sense";
system "cat $ARGV[0].antisense | sort | uniq >uniq.$ARGV[0].antisense";
system "cat $ARGV[0].antisense $ARGV[0].sense | sort | uniq >uniq.$ARGV[0].total";
system "grep -Fxf uniq.$ARGV[0].sense uniq.$ARGV[0].antisense >uniq.$ARGV[0].sense.antisense.overlap"; 
system "wc -l uniq.$ARGV[0].*";


sub array_chunk {
   my($num,@arr)=@_;
   my @ret;
   while(@arr){
     push(@ret,[splice@arr,0,$num]);
   }
   return @ret;
}

