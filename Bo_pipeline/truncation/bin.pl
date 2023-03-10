#!/usr/bin/perl -w

#$len=4824;
open A,"$ARGV[0]"||die $!;
$all=join'',<A>;
@each=split(/\n/,$all);
my @pos;
my @fouth;
my @fifth;
my @eleventh;
my @twelevth;

foreach $line(@each){
      @aa=split(/\t/,$line);
      push (@pos,$aa[1]);
      push (@fouth,$aa[4]);
      push (@fifth,$aa[5]);
      push (@eleventh,$aa[11]);
      push (@twelevth,$aa[12]);      
}


#print "@st";
# print "@nd";   
#my $count=@st;
#print "$count\n";

for ($i=0;$i<4808;$i++){
      $position=($pos[$i]+$pos[$i+4])/2;
      $IR_plus=($fouth[$i]+$fouth[$i+1]+$fouth[$i+2]+$fouth[$i+3]+$fouth[$i+4]);
      $IR_minus=($fifth[$i]+$fifth[$i+1]+$fifth[$i+2]+$fifth[$i+3]+$fifth[$i+4]);
      $terminal_plus=($eleventh[$i]+$eleventh[$i+1]+$eleventh[$i+2]+$eleventh[$i+3]+$eleventh[$i+4]);
      $terminal_minus=($twelevth[$i]+$twelevth[$i+1]+$twelevth[$i+2]+$twelevth[$i+3]+$twelevth[$i+4]);
      $average_IRplus=$IR_plus/5;
      $average_IRminus=$IR_minus/5;
      $average_terminalplus=$terminal_plus/5;
      $average_terminalminus=$terminal_minus/5;     

       print "$position\t$average_IRplus\t$average_IRminus\t$average_terminalplus\t$average_terminalminus\n";}

close A;

