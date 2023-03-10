#get all GO childs and parents for any GO term  
#2018/05/05 edited to not use "GO::Parser"

#!/usr/bin/perl -w

require "../../SharedCodes/perl.fun.inc.pl";
#use GO::Parser;
$goVersion="";
$go_defi_f="01.download/go-basic.obo";
#my $parser = new GO::Parser({handler=>'obj'}); # create parser object
#$parser->parse("$go_defi_f"); # parse file -> objects

my $output = "02.GO_assn/$goVersion/go2go.all.out";

#########RUN


create_dir_ifNotExist($output);
my %cnt;

##1, read go_defi_f and load hashs %go2term_h, %go2cat_h, %go2parent_h, %go2child_h, %all_goids_h
open (GOOBO,$go_defi_f) || die "error open $go_defi_f\n";
my $current_go;
while(<GOOBO>){
	s/\s+$//;
	if(s/^id: //){
		$current_go=$_;
		next if $current_go !~/^GO:/;
		$all_goids_h{$current_go}++;
	}
	if(s/^name: //){ $go2term_h{$current_go}=$_; }
	if(s/^namespace: //){ $go2cat_h{$current_go}=$_; }
	if(s/^is_a: |^relationship: part_of //){
		s/ ! .*//;
		my $parent_GO=$_;
		next if $parent_GO!~/^GO:/;
		$go2parent_h{$current_go}{$parent_GO}++;
		$go2child_h{$parent_GO}{$current_go}++;
		$all_goids_h{$parent_GO}++;
	}
}
close GOOBO;

##2, search all parents and all childs through the GO tree
foreach my $goid (sort keys %go2parent_h){
	print "search parents for $goid\n";
	search_all_parents($goid,$goid);
}
#

sub search_all_parents{ #
	my (@goids)=@_;
	my $child_id=shift @goids;
	my %old_parent_h=();
	map($old_parent_h{$_}++, @goids);
	my %new_parents_h=();
	foreach my $goid(@goids){
		my @parents_gos=keys %{$go2parent_h{$goid}};
		map($new_parents_h{$_}++, @parents_gos);
	}
	
	foreach my $parent_id(keys %new_parents_h){
		next if $old_parent_h{$parent_id};
		$go2parent_h{$child_id}{$parent_id}++;
		$go2child_h{$parent_id}{$child_id}++;
		search_all_parents($child_id,$parent_id);
	}
}


##3, output:
open (OUT, "> $output");	
print OUT "#Category\tGO_ID\tGO_term\tNum_Parent\tNum_Child\tALL_Children\n";

#my $graph = $parser->handler->graph;  # get L<GO::Model::Graph> object
#my $all_terms = $graph->get_all_nodes();
my @all_goids=sort keys %all_goids_h;

foreach my $goid (@all_goids) {

	# $cnt{'all'}++;
	# my $nsp = $term->namespace;
	# if ($term->is_obsolete) {
	# 	$cnt{'obs'}++;
	# 	next;
	# }
	# my $parent_terms = $graph->get_recursive_parent_terms($term->acc);
	# my $child_terms = $graph->get_recursive_child_terms($term->acc);
	# my @a = &proc($parent_terms);   #remove duplicate terms
	# my @b = &proc($child_terms);	#remove duplicate terms
	# my $child_list = join ("|",@b);
	# if ($child_list eq "") { $child_list = "-"; }
	# my $line = sprintf "%s\t%s\t%s\t%s\t%s", $term->acc, $term->name, scalar @a, scalar @b, $child_list;
	# printf OUT "$nsp\t$line\n";
	# $cnt{"OUTPUT: $output"}++;
	my $num_child=scalar keys %{$go2child_h{$goid}}; 
	$num_child=0 if !$num_child;
	my $child_gos=join("|", sort keys %{$go2child_h{$goid}});
	$child_gos="-" if !$child_gos;
	print OUT $go2cat_h{$goid}."	$goid	".$go2term_h{$goid}."	".(scalar keys %{$go2parent_h{$goid}})."	$num_child	$child_gos\n";
}
close OUT;
