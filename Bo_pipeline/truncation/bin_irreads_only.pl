#!/usr/bin/perl -w
##usage:  perl bin_irreads_only.pl tmp;
#$len=4824;
open A,"$ARGV[0]"||die $!;
$all=join'',<A>;
@each=split(/\n/,$all);
my @pos;
my @rd;
my @fouth;
my @fifth;
my @eleventh;
my @twelevth;

print "Start_position\tNumber_of_total_reads\tIR_read(+)\tIR_read(-)\t%read_terminal_+\t%read_terminal_-\n";
#shift(@each);
foreach $line(@each){
      @aa=split(/\t/,$line);
      push (@pos,$aa[1]);
      push (@rd,$aa[2]);
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
      #$coverage=($rd[$i]+$rd[$i+1]+$rd[$i+2]+$rd[$i+3]+$rd[$i+4]);
      $coverage=$rd[$i+2];
      $IR_plus=($fouth[$i]+$fouth[$i+1]+$fouth[$i+2]+$fouth[$i+3]+$fouth[$i+4]);
      $IR_minus=($fifth[$i]+$fifth[$i+1]+$fifth[$i+2]+$fifth[$i+3]+$fifth[$i+4]);
      #$terminal_plus=($eleventh[$i]+$eleventh[$i+1]+$eleventh[$i+2]+$eleventh[$i+3]+$eleventh[$i+4]);
      $terminal_plus=$eleventh[$i+2];
      #$terminal_minus=($twelevth[$i]+$twelevth[$i+1]+$twelevth[$i+2]+$twelevth[$i+3]+$twelevth[$i+4]);
      $terminal_minus=$twelevth[$i+2];
      #$average_coverage=$coverage/5;
      $average_IRplus=$IR_plus/5;
      $average_IRminus=$IR_minus/5;
      #$average_terminalplus=$terminal_plus/5;
      #$average_terminalminus=$terminal_minus/5;     

      # print "$position\t$average_coverage\t$average_IRplus\t$average_IRminus\t$average_terminalplus\t$average_terminalminus\n";}
       print "$position\t$coverage\t$average_IRplus\t$average_IRminus\t$terminal_plus\t$terminal_minus\n";}
close A;

