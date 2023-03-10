#!/usr/bin/perl -w
#
open A,"$ARGV[0]"||die $!;
while(<A>){
      chomp;
      @aa=split(/\t/,$_);
       if(($aa[5]=~/\-/)&&($aa[5]!~/\+/)){
        print "$aa[0]\t$aa[1]\t$aa[2]\t$aa[3]\t$aa[4]\t\-\t$aa[2]\t\+\t0\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\t$aa[10]\n";}
      elsif(($aa[5]=~/\+/)&&($aa[5]!~/\-/)){
        print "$aa[0]\t$aa[1]\t$aa[2]\t$aa[3]\t$aa[4]\t\+\t$aa[2]\t\-\t0\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\t$aa[10]\n";}
      elsif(($aa[5]=~/\+/)&&($aa[5]=~/\-/)){
        print "$aa[0]\t$aa[1]\t$aa[2]\t$aa[3]\t$aa[4]\t\+\t$aa[2]\t\-\t$aa[2]\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\t$aa[10]\n";}
      elsif($aa[5]=~/\./){
        print "$aa[0]\t$aa[1]\t$aa[2]\t$aa[3]\t$aa[4]\t\+\t0\t\-\t0\t$aa[6]\t$aa[7]\t$aa[8]\t$aa[9]\t$aa[10]\n";}
}

close A;


