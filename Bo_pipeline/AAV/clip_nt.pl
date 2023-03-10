#!/usr/bin/perl -w

open A,"$ARGV[0]"||die $!;
my @clip;
while(<A>){
	chomp;
	@aa=split(/\t/,$_);
	if($aa[4]>30){
            if($aa[5]!~/(\d+)H/){
                if($aa[5]=~/^(\d+)S.*\D(\d+)S$/){
                $clipped_nt=$1+$2;
 }#           print "$clipped_nt\t";}
                elsif(($aa[5]=~/^(\d+)S/)&&($aa[5]!~/(\d+)S$/)){
                $clipped_nt=$1;
}#print "$clipped_nt\t";}
                elsif(($aa[5]!~/^(\d+)S/)&&($aa[5]=~/(\d+)S$/)){
                 $clipped_nt=$1;
}#print "$clipped_nt\t";}

      @bb=split(/\D+/,$aa[5]);
      $sum=0;
      foreach $num(@bb){
         $sum+=$num;}
    #  print "$sum\t";
      $percentage_clip=$clipped_nt/$sum;
   #   print "$percentage_clip\n";
      if($percentage_clip<0.4){
#      print "$_\t$clipped_nt\t$percentage_clip\n";
      print "$_\n";
}
}
}
}


close A;


