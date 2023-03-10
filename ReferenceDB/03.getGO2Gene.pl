# from jiyeon's script
#from NCBI gid2go file and go anotation file pasing from .obo file in http://www.geneontology.org/,
#create gene IDs to all GO IDs file

#!/usr/bin/perl -w

require "../../SharedCodes/perl.fun.inc.pl";
my $db = "gunzip -c 01.download/gene2go.gz |"; 
my $input = "02.GO_assn/go2go.all.out";
my $out_folder="03.go2gene/";

use Getopt::Std; 
getopt("st",\%args);

$sp = "hs";  $tax = "9606";
$sp = "mm";	$tax = "10090";
$sp = "rn";  $tax = "10116";
$sp = "chicken";  $tax = "9031";

$sp=$args{s} if $args{s};
$tax=$args{t} if $args{t};

my $output = "$out_folder/$sp.go2gene.tbl";
my $output2 = "$out_folder/$sp.go2gene.num.tbl";

########################################################
create_dir_ifNotExist($output);
print "$sp, $tax\n";

my (%cnt, %h, %h2);
open (DB, $db) || die $!;
while (<DB>) {
	chomp;
	next if (/^\#/);
	$cnt{"DB: $db"}++;

	my @items = split /\t/;
	my $tid = &trim($items[0]);
	next if ($tid ne $tax);  #restrict with species

	my $gid = &trim($items[1]);
	my $go_id = &trim($items[2]);

	$h{$go_id}{$gid} = 1;
}
close DB;

open (OUT, "> $output") || die $!;
open (OUT2, "> $output2") || die $!;
print OUT "GO_ID\tGene_ID\n";
print OUT2 "Category\tGO_ID\tGO_term\tNum_Child\tNum_Parent\tNum_Genes\n";

open (IN, $input) || die $!;
while (<IN>) {
	chomp;
	next if (/^\#/);

	$cnt{"INPUT: $input"}++;
	my @items = split /\t/;
	if (scalar @items != 6) { print "error1: $_\n"; }
	my ($cat, $go_id, $go_name, $num_pa, $num_ch, $assn_chs)=@items;

	if ($go_id !~ /^GO:/){ print "error2: $_\n"; }

	my @cgo = split /\|/, $assn_chs;
	if ($cgo[0] eq "-") { shift @cgo; }
	if ($num_ch != scalar @cgo) { print "error3: $_\n"; }

	my @gs=keys %{$h{$go_id}};
	foreach my $cgo_id (sort @cgo) {
		if ($cgo_id !~ /^GO:/){ print "error4: $_\n"; }

		push (@gs, keys %{$h{$cgo_id}});  #add genes of all the children
	}
	my %ug = map { $_, 1 } @gs; #remove duplicates

	foreach my $gid (sort {$a <=> $b} keys %ug){
		$cnt{"OUTPUT1: $output"}++;
		print OUT "$go_id\t$gid\n";
	}

	$cnt{"OUTPUT2: $output2"}++;
	print OUT2 "$cat\t$go_id\t$go_name\t$num_ch\t$num_pa\t", scalar keys %ug, "\n";
}
close IN; close OUT; close OUT2;


foreach my $tmp (sort keys %cnt) {
	print "$tmp\t$cnt{$tmp}\n";
}
sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
