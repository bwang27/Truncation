#!/usr/bin/perl -w
#usage: perl judge.the.read.inverted.or.not.pl plus.minus.step4.scAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments;
#perl modify.judge.the.read.inverted.or.not.pl plus.minus.step4.ssAAV-SCA3-4.step1.same.genome.step3.sameread.all.alignments 10;
#perl judge.reads.position.continous.or.not.pl IR.sense.reads.position 10
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
open SENSE,">$ARGV[0].$ARGV[1].both.confirmed.continuous"||die $!;
#open ANTISENSE,">$ARGV[0].antisense"||die $!;
while(<B>){
   chomp;
   @nd = split(/\t/,$_);
   $read = $nd[0];
   shift(@nd);
 #  $first = shift(@nd);
 #  @seperate_first = split(/\s/,$first);
     for ($i=0;$i<@nd;$i++){
        $qry=$nd[$i];
        @split_qry=split(/\s/,$qry);
        @qrystart=split(/\//,$split_qry[0]);
        @qryend=split(/\//,$split_qry[1]);
        foreach $element(@nd){
         if($element ne $qry){
          #  @split_qry=split(/\s/,$qry);
            @seperate_element = split(/\s/,$element);
            @targetstart=split(/\//,$seperate_element[0]);
            @targetend=split(/\//,$seperate_element[1]);
#        if(($seperate_first[2]!~/\$seperate_element[2]/)&&($seperate_first[1]==$seperate_element[1])){
             $leftread=$targetend[0]-"$ARGV[1]"; #30; #20;
             $rightread=$targetend[0]+"$ARGV[1]"; #30; #20;
             $left_startread=$targetstart[0]-"$ARGV[1]"; #30;
             $right_startread=$targetstart[0]+"$ARGV[1]"; #30;
             $leftread_ir=$targetend[1]-"$ARGV[1]";
             $rightread_ir=$targetend[1]+"$ARGV[1]";
             $left_startread_ir=$targetstart[1]-"$ARGV[1]";
             $right_startread_ir=$targetstart[1]+"$ARGV[1]";
          # if(($split_qry[2] ne $seperate_element[2])&&((($qryend[0]>$left_startread)&&($qryend[0]<$right_startread))||(($qrystart[0]>$leftread)&&($qrystart[0]<$rightread)))&&((($qryend[1]>$leftread_ir)&&($qryend[1]<$rightread_ir))||(($qrystart[1]>$left_startread_ir)&&($qrystart[1]<$right_startread_ir)))){
   #   print SENSE "$_\n";}
       if(($split_qry[2] ne $seperate_element[2])&&((($qryend[0]>$left_startread)&&($qryend[0]<$right_startread))||(($qrystart[0]>$leftread)&&($qrystart[0]<$rightread)))&&((($qryend[1]>$leftread_ir)&&($qryend[1]<$rightread_ir)))){
  
     print SENSE "$read\t$split_qry[0]\t$split_qry[1]\t$split_qry[2]\t$seperate_element[0]\t$seperate_element[1]\t$seperate_element[2]\n";}    
      # elsif(($split_qry[2]!~/\$seperate_element[2]/)&&($split_qry[0]>$left_start)&&($split_qry[0]<$right_start)){
    #  print ANTISENSE "$_\n";}

          
  		}
 	}	
  }
}

close A;
close B;

#system "cat $ARGV[0].sense | sort | uniq >uniq.$ARGV[0].sense";
system "cat $ARGV[0].$ARGV[1].both.confirmed.continuous | sort | uniq >uniq.$ARGV[0].$ARGV[1].both.confirmed.continuous";
system "awk '\$4~/+/ {print}' uniq.$ARGV[0].$ARGV[1].both.confirmed.continuous >final.uniq.$ARGV[0].$ARGV[1].both.confirmed.continuous";

#system "cat $ARGV[0].antisense $ARGV[0].sense | sort | uniq >uniq.$ARGV[0].total";
#system "grep -Fxf uniq.$ARGV[0].sense uniq.$ARGV[0].antisense >uniq.$ARGV[0].sense.antisense.overlap"; 
#system "wc -l uniq.$ARGV[0].*";
open FINAL,"final.uniq.$ARGV[0].$ARGV[1].both.confirmed.continuous"||die $!;
open OUT,">$ARGV[0].$ARGV[1].both.confirmed.continuousfor.downstream.analysis"||die $!;
while(<FINAL>){
    chomp;
    @aa=split(/\t/,$_);
    @bb=split(/\//,$aa[1]);
    @cc=split(/\//,$aa[2]);
    @dd=split(/\//,$aa[4]);
    @ee=split(/\//,$aa[5]);
  print OUT "$aa[0]\t$bb[1]\t$cc[1]\t$aa[3]\t$dd[1]\t$ee[1]\t$aa[6]\n";}

close FINAL;
close OUT;



sub array_chunk {
   my($num,@arr)=@_;
   my @ret;
   while(@arr){
     push(@ret,[splice@arr,0,$num]);
   }
   return @ret;
}

