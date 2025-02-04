#perl functions
#require "/Home/cocochen/analyze/perl.fun.inc.pl";
use POSIX;

sub change_f_name_if_exist{
	my($file)=@_;
	my $counter=0;
	my $file_out=$file;
	while( -e $file_out){
	  $counter++;
	  $file_out = "$file.$counter";
	}
	$file_out;
}

sub create_dir_ifNotExist{ #input xxxx/xxxxx/  or xxx/xxx/xxx.xxx or /xxx/xxx/xxx
	$newdir=shift;
	$newdir=~s/[^\/]*$//;
	print "newdir=".$newdir."\n";
	if(!(-d $newdir)){
		@folder_arr=split(/\//,$newdir);
		my $i;
		for($i=0;$i<@folder_arr;$i++){
			$dir_i=join("/",@folder_arr[0..$i]);
			next if (!$dir_i); ##root
			if(! (-d $dir_i)){
				print "mkdir $dir_i, ";
				print(mkdir($dir_i,0774)."\n");
			}
		}
	}
}


sub operate_str_array{
	# $operator= "mi" or "ad"
	my ($str1,$str2,$seperator,$operator,$join_with_str)=@_;
	if(! $str1){return $str2;}else{
		my @str1_arr=split(/$seperator/,$str1);
		my @str2_arr=split(/$seperator/,$str2);
		my $outstr="";
		for(my $i=0;$i<@str1_arr;$i++){
			if($i>0){$outstr.=$seperator;}
			$str1_arr[$i]=0 if (!$str1_arr[$i]);
			$str2_arr[$i]=0 if (!$str2_arr[$i]);
			if($operator eq "ad" || $operator eq "+"){ #add
				$outstr.=$str1_arr[$i]+$str2_arr[$i];
			}elsif($operator eq "mi" || $operator eq "-"){ #minus
				$outstr.=$str1_arr[$i]-$str2_arr[$i];
			}elsif($operator eq "div"  || $operator eq "/"){
				$outstr.=$str1_arr[$i]/$str2_arr[$i];
			}elsif($operator eq "mu"  || $operator eq "*"){
				$outstr.=$str1_arr[$i]*$str2_arr[$i];
			}elsif($operator eq "jo"  || $operator eq "."){ #join
				$outstr.=$str1_arr[$i].$join_with_str.$str2_arr[$i];
			}
			
		}
		return $outstr;
	}
}

sub operate_multi_str_array{ #more than 2 string operation
	my ($seperator,$operator,$join_with_str,@str_arr)=@_;
	my $outstr=$str_arr[0];
	for(my $i=1;$i<@str_arr;$i++){
		$outstr=operate_str_array($outstr,$str_arr[$i], $seperator,$operator,$join_with_str);
	}
	return($outstr);
}

sub extract_seq_fr_contig{ 
	my ($contig_file_dir,$contig_name,$strand,$regions_str,$contig_file)=@_;
	#if AACCCTTT + 3:6, return CCCT
	#if AACCCTTT - 3:6, return AGGG
	#$contig_name is the contig name following ">" in the contig file
	#when $search_reverse eq "T" and the strand is minus, program will search the reversed strand. the output will be reversed complementary NTs
	
	if(!$contig_seq{$contig_name} && !(scalar keys %contig_seq>0 && $contig_file) ){ #not loaded the seq
		#if (keys %contig_seq >6){%contig_seq={};}; #to clear memory, use this if the contig is very large
		$contig_file="$contig_file_dir/$contig_name.fa" if !$contig_file; 
		open (CONTIG_FILE_IN, $contig_file) || die ("error $contig_file\n"); #read the contig file and get the seq
		print "$contig_file opened...\n";
		my $contig_name_infile='';
		while (<CONTIG_FILE_IN>){
			s/\s+$//;
			if(substr($_,0,1)eq">"){
				$_=~s/>//;
				$contig_name_infile=$_;
				print " find sequence $contig_name_infile;\n";
				next;
			}
			$contig_seq{$contig_name_infile}.=$_;
		}
		close CONTIG_FILE_IN;
	}
	
	#extract seq from $contig_seq{$contig_name}
	my @regions=split(/,|\|/,$regions_str); #for multi-region (like junction)
	my $reg_id=0,$out_seq="";
	if($contig_seq{$contig_name}){
		while ($reg_id<@regions){
			my ($from,$to)=sort {$a<=>$b} split(/:|\-/,$regions[$reg_id]);
			my $seq_length=$to-$from+1;
			$out_seq.=substr($contig_seq{$contig_name},$from-1,$seq_length);
			$reg_id++;
		}
		if($strand=~/\-/){ # to reverse and complement
			$out_seq=reverse_complement($out_seq,"F");
		}
	}
	return uc($out_seq);
}

sub reverse_complement{ #for DNA or RNA sequence only
	my ($input_seq,$is_RNA)=@_;
	$input_seq=~tr/actguACTGU/tgacaTGACA/;
	if($is_RNA eq "T"){$input_seq=~tr/tT/uU/;}
	my $out_seq=reverse_str($input_seq);
	#for($i=length($input_seq)-1;$i>=0;$i--){
	#	$out_seq.=substr($input_seq,$i,1);
	#}
	return $out_seq;
}

sub reverse_str{
	my $input=shift;
	join("", reverse(split(//,$input)) )
}

sub flip_region_str{ #convert 12-34|56-78 to 78-56|34-12
	my $in_region_str=shift;
	$in_region_str=~s/([\|\-:])/	$1	/g;
	join("", reverse(split(/	/,$in_region_str)));	
}

sub get_region_bundary{
	my $in_region_str=shift;
	my @coor_arr=sort {$a<=>$b} split(/\||\-/,$in_region_str);
	return(@coor_arr[0,-1]);
}

sub return_time{
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	return (($year+1900)."/".($mon+1)."/".$mday." ".$hour.":".$min.":".$sec);
}



sub judge_2region_relation{ 
	#given 2 regions a and b (4 coordinates), output position of b relative to a (return =, 0, 1, -1, inc, ov1, ov-1)
	#if $ifout_overlapLen=1, output overlap info (overlap vode, overlap length, overlap coordinates and unoverlapped coordinates)
	my ($a1,$a2,$b1,$b2,$strand_a,$ifout_overlapLen)=@_;
	$strand_a=1 if ($strand_a eq "+");
	$strand_a=-1 if ($strand_a eq "-");
	$strand_a=1 if $strand_a ne '-1'; #$strand_a= 1 or -1 (represent + or - strand)
	if($a1>$a2){($a1,$a2)=($a2,$a1);  }
	($b1,$b2)=sort {$a<=>$b} ($b1,$b2);
	if(($a1 eq $b1 && $a2 eq $b2) || ($a1 eq $b2 && $a2 eq $b1) ){
		return $ifout_overlapLen ? ("=", $a2-$a1+1, "$a1:$a2", "") : "=";
	}else{
		$b1pos=judge_1pos_rel_1reg($a1,$a2,$b1);
		$b2pos=judge_1pos_rel_1reg($a1,$a2,$b2);
		if($b1pos==$b2pos){
			if($ifout_overlapLen){
				if($b1pos eq "0"){
					my @non_overlapped_a_arr=();
					if($b1>$a1){push @non_overlapped_a_arr, "$a1:".($b1-1);}
					if($b2<$a2){push @non_overlapped_a_arr, ($b2+1).":$a2";}
					return ($b1pos, abs($b2-$b1)+1,"$b1:$b2", join("|",@non_overlapped_a_arr) );
				}else{
					return ($b1pos*$strand_a, 0, "", "$a1:$a2");
				}
				
			}else{
				return $b1pos*$strand_a; #0 1 or -1
			}
			
		}else{
			if($b1pos*$b2pos==-1 || ($b1pos==-1 && $b2==$a2) ||  ($b2pos==1 && $b1==$a1) ){
				return $ifout_overlapLen ? ("inc", $a2-$a1+1, "$a1:$a2", "") : "inc"; #b include a
			}else{ #overlap
				my $posStr_rel=$b1pos+$b2pos; #if positive strand, relation: 1 or -1
				if($ifout_overlapLen){
					my @arr4=sort {$a<=>$b}($a1,$a2,$b1,$b2);
					my $overlap_len=$arr4[2]-$arr4[1]+1;
					return ("ov".($posStr_rel*$strand_a) , $overlap_len,  $arr4[1].":".$arr4[2],  $posStr_rel==1? "$a1:".($b1-1) : ($b2+1).":$a2"   );
				}else{
					return ( "ov".($posStr_rel*$strand_a) ); #ov1 or ov-1
				}
			}
		}		
	}
}


%match_syb2val_hash=(
	"="=>5,
	"0"=>4,
	"inc"=>3,
	"ov1"=>2,
	"ov-1"=>2,
	"1"=>1,
	"-1"=>1,
);

sub judge_1pos_rel_1reg{ #given 1 region a (a1 to a2) and another position x, judge position x relative to a (return -1, 0 or 1)
	my ($a1,$a2,$x,$strand_a)=@_; #a1 should<a2
	($a1,$a2)=sort{$a<=>$b}($a1,$a2);
	$strand_a=1 if ($strand_a eq "+");
	$strand_a=-1 if ($strand_a eq "-");
	$strand_a=1 if !$strand_a;
	if( ($a1-$x)*($a2-$x)<=0 ){
		return 0; #inside
	}else{
		if($x<$a1){
			return -1*$strand_a; #upstream if strand = 1 (+)
		}else{
			return $strand_a; #downstream if strand = 1 (+)
		}
	}
}

sub create_chr_reg_info_hash{ #create standard indexed chromosome region info hash
	my ($pos1,$pos2,$key1,$val_str,$chr_block_size,$key2)=@_; 
	#$key1 usally contains chromosome, [strand]
	#$key2 usally is a region (it should be unique for all values)
	my $posfrom=($pos1<$pos2)?$pos1:$pos2;
	my $posto=($pos1<$pos2)?$pos2:$pos1;
	my $posfrom_block_num=int($posfrom/$chr_block_size);
	my $posto_block_num=int($posto/$chr_block_size);
	$key2="$posfrom:$posto" if (!$key2);
	foreach my $pos_block_num(($posfrom_block_num..$posto_block_num)){
		$chro2reginfo_hash{"$key1:$pos_block_num"}{$key2}=$val_str;
	}
	($key1,$posfrom_block_num,$posto_block_num,$key2);
}


###############.sam parsing functions

sub cal_coverage_fr_cigar{ #calculate reads coverage region and junction/intron info from cigar string
	my ($from_pos,$ciger_str)=@_;
	return "" if $ciger_str!~/M/i;
	my @exonR_arr=(); #exon region
	my @intron_arr=();
	my $match_len;
	$ciger_str=~s/\d+[^\dMIDNSHP]//g;
	while($ciger_str){
		if($ciger_str=~s/^(\d+)M//){ #aligned region
			$match_len=$1;
			if(scalar @exonR_arr>0 && $exonR_arr[@exonR_arr-1] ne "|"){
				$exonR_arr[@exonR_arr-1]+=$match_len-1;
			}else{
				push (@exonR_arr, "$from_pos:", ($from_pos+$match_len-1) );
			}
			$from_pos+=$match_len-1;
		}
		if($ciger_str=~s/^(\d+)N//){ #gap, intron
			$match_len=$1;
			push (@intron_arr, ($exonR_arr[@exonR_arr-1]+1).":". ($exonR_arr[@exonR_arr-1]+$match_len) );
			push (@exonR_arr, "|") if @exonR_arr>0;
			$from_pos+=$match_len+1;
		}
		if($ciger_str=~s/^(\d+)[SHPI]//){
			#do nothing
		}
		if($ciger_str=~s/^(\d+)D//){ #deletion in reads
			$match_len=$1;
			$exonR_arr[@exonR_arr-1]+=$match_len+1 if @exonR_arr>0;
			$from_pos+=$match_len+1;
		}
	}
	return (join("",@exonR_arr), join("|",@intron_arr));
}


sub from_flag2strand { ###now only for tophat single end outfile
	my $flag_str=shift;
	($flag_str & 16)?"-":"+";
}

sub from_flag2info{
	my $flag_str=shift;
	#print "$flag_str	";
	my $if_PE=($flag_str & 1)?"PE":"SE";
	my $ifProperPair=($flag_str & 2)?1:0;
	my $ifUnmapped=($flag_str & 4)?1:0;
	my $Rstrand=($flag_str & 16)?"-":"+";
	my $mRstrand=$if_PE eq "PE"? (($flag_str & 0x0020)?"-":"+") : "";
	my $ifNotPrimary=($flag_str & 0x0100)?1:0;
	#print "$if_PE	$ifProperPair	$ifUnmapped	$Rstrand	$mRstrand	$ifNotPrimary\n";
	return($if_PE,$ifProperPair,$ifUnmapped,$Rstrand,$mRstrand,$ifNotPrimary);
}

sub judge_reads_quality{
	my ($OutSamFrom,$mismatch_cut,$clip_cut,$cigar,@tags_arr)=@_;
	##1, if clip length > clip_cut (bad quality)
	if($clip_cut ne 'none'){ #eg. set $clip_cut ='none' to invalid  this parameter
		my $clip_len=0;
		my $sum=0;
		my $clip_percent=0;
		while($cigar!~/(\d+)H/){
                if($cigar=~/^(\d+)S.*\D(\d+)S$/){
                $clip_len=$1+$2;
        }#           print "$clipped_nt\t";}
                elsif(($cigar=~/^(\d+)S/)&&($cigar!~/(\d+)S$/)){
                $clip_len=$1;
        }#print "$clipped_nt\t";}
                elsif(($cigar!~/^(\d+)S/)&&($cigar=~/(\d+)S$/)){
                 $clip_len=$1;}
        @bb=split(/\D+/,$cigar);
        foreach $num(@bb){
        $sum+=$num;}
        $clip_percent=$clip_len/$sum;
        if($clip_percent>$clip_cut){
            return 0;
        } else {return 1;}
}}

	##2, if mismatch > $mismatch_cut (bad quality)
	if($mismatch_cut ne 'none'){
		my $tag_str=join("	",@tags_arr);
		my $mismatch=0;
		if($OutSamFrom eq "bwasw"){ #$cigar=1S22M  AS:i:22
			$cigar=~/(\d+)M/;
			my $match_len=$1;
			if($tag_str=~/AS:i:(\d+)/){
				$mismatch=($match_len-$1)/4 ; #1 mismatch score =-3, match score =1 for bwasw
			}
		}else{
			if($tag_str=~/NM:i:(\d+)/i){$mismatch=$1;} #NM:i:1
		}
		if($mismatch>$mismatch_cut){return 0;}else{return 1;}
	}
	return 1;
}

sub cigar2cover_len{
	my $cigar2=shift;
	$cigar2=~s/\d+[HSIP]//g;
	my $total_cover_len=0; #length in reference covered by reads aligned region
	foreach $len1(split(/[MNDX=]/,$cigar2)){	$total_cover_len+=$len1;} #
	$total_cover_len;
}
sub cigar2readMapCorr{ #from cigar, get te information of which region of a read was mapped to the reference
	my $cigar2=shift;
	$cigar2=~s/\d+[NDP]//g;
	$cigar2=~s/\d+H$//; #remove right-hand hard clip
	$cigar2=~s/\d+S$//; #remove right-hand soft clip
	my $readCovLen=0; #length of a read in aligned region
	my $leftClipLen=0; #length of read left clips
	while($cigar2=~/(\d+)[SH]/g){ $leftClipLen+=$1;}
	while($cigar2=~/(\d+)[MIX=]/g){ $readCovLen+=$1;}
	return($leftClipLen, $leftClipLen+$readCovLen);
}

sub cigar2read_len{
	my $cigar2=shift;
	$cigar2=~s/\d+[NDP]//g;
	my $readLen=0; #length of a read
	foreach $len1(split(/[HSMIX=]/,$cigar2)){	$readLen+=$len1;} #
	$readLen;
}

###################
sub sum{
	my @invars=@_;
	my $sumval=0;
	map($sumval+=$_,@invars);
	$sumval;
}
sub average{
	my @invars=@_;
	return sum(@invars)/(scalar @invars) ;
}
sub stdev{
	my(@data) = @_;
	if(@data == 1){
	        return 0;
	}
	my $average = &average(@data);
	my $sqtotal = 0;
	foreach(@data) {
	        $sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@data-1)) ** 0.5;
	return $std;
}

sub sem{ #standard error of mean
	my(@data) = @_;
	return stdev(@data)/ (scalar @data ** 0.5) ;
}

##remove array elements
sub rm_arr_elements{
	my ($arr_addr,@rm_ids)=@_;	
	my @new_arr=();
	my %rm_id_h=();
	map($rm_id_h{$_}=1,@rm_ids);
	for(my $i=0;$i<@{$arr_addr};$i++){
		if(!$rm_id_h{$i}){
			push @new_arr,$$arr_addr[$i];
		}
	}
	@new_arr;
}


###################polyA related
sub find_pA_signal{ 
	my ($in_seq)=@_;
	$in_seq=~tr/uU/tT/;
	my %occurance_num_hash=(); 
	my %occurance_detail_hash=();
	my ($match,$signal_pos);
	
	my @pas_motif_types=("AATAAA","ATTAAA","OtherPAS","Arich");
	my %motif_search_str_hash=(
		"AATAAA"=>"AATAAA",
		"ATTAAA"=>"ATTAAA",
		"OtherPAS"=>"TATAAA|AGTAAA|AATATA|AATACA|CATAAA|GATAAA|AATGAA|TTTAAA|ACTAAA|AATAGA", #AAGAAA take out
		"Arich"=>"AAAAA|A[CTG]AAAA|AA[CG]AAA|AAA[TCG]AA|AAAA[CTG]A",
	);
	
	foreach $motif_type( keys %motif_search_str_hash){
		$search_str=$motif_search_str_hash{$motif_type};
		while($match= $in_seq=~/$search_str/ig){
			$signal_pos=length($in_seq)-pos($in_seq);
			$occurance_num_hash{$motif_type}++;
			$occurance_detail_hash{$motif_type}.="$signal_pos,";
		}
	}
	
	$pas_type="NoPAS";
	for($i=@pas_motif_types-1;$i>=0;$i--){
		if($occurance_num_hash{$pas_motif_types[$i]}){
			$pas_type=$pas_motif_types[$i];
		}
	}
	
	return (join("	",map($occurance_num_hash{$_}."	".$occurance_detail_hash{$_},@pas_motif_types)), $pas_type);
}


%readsTypes2sco_h=( #higher score represent higher priority (applicable to RNA-seq); + - represent strand relative to reference gene only useful for directional reads 
	"CDS+"=>18,
	"CDS"=>17,
	"3UTR+"=>16,
	"3UTR"=>15,
	"5UTR+"=>14,
	"5UTR"=>13,
	"UTR3e+"=>12,
	"UTR3e"=>11,
	"intron+"=>10,
	"intron"=>9,
	"intron-"=>8,
	"UTR5e-"=>7,
	"UTR5e"=>6,
	"UTR5e+"=>5,
	"UTR3e-"=>4,
	"3UTR-"=>3,
	"CDS-"=>2,
	"5UTR-"=>1,
	"intergenic"=>0,
);

sub get_genomic_region{
	my ($ref_pos, $direction, $rel_pos1, $rel_pos2)=@_; #$direction = 1 or -1
	my $pos1=$ref_pos+$direction*$rel_pos1;
	my $pos2=$ref_pos+$direction*$rel_pos2;
	$pos1=1 if $pos1<=0;
	$pos2=1 if $pos2<=0;
	if($pos1>1 || $pos2>1){
		join(":", sort {$a<=>$b}($pos1,$pos2));
	}else{
		"";
	}
}


sub convert_regName2_coordinates{ #like: p101p2100, m100m60, convert to (101,2100) (-100,-60)
	my $region_name=shift; 
	$region_name=~s/([mp])/	$1/g;
	$region_name=~s/m/\-/g;
	$region_name=~s/p//g;
	$region_name=~s/^	//;
	split(/	/,$region_name);
}

%Gmodl2header_h=(
	"refFlat"=>"gene_symbol|name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends",
	"knownGene"=>"name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends|proteinID|alignID",
	"ensGene"=>"bin|name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends|score|name2|cdsStartStat|cdsEndStat|exonFrames", #for mm9 hg19 hg18
	"ensGene2"=>"name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends", #for mm8

);

sub bed2bigbed{
	my ($geno,$bedf,$group_name,$bed_name)=@_;
	my $bbf_name="$bedf.bb";
 	print "\n####convert $bedf to bigbed file: $bbf_name\n";
	system("sed \'1d\' $bedf | sort -k1,1 -k2,2n  > $bedf.sorted");
	system("bedToBigBed $bedf.sorted /Home/cocochen/soft/biosoft/genome/UCSC/$geno.chrom.sizes $bbf_name");
	system("rm $bedf.sorted");
	$bbf_name=~s/^.*\///;
	print  "####add track to UCSC:####\ntrack type=bigBed itemRgb=On visibility=3 colorByStrand=\'255,0,0 0,0,255\' group=$group_name priority=1 name=\"$bed_name\" description=\"$bed_name\" bigDataUrl=http://rna:rna\@rna.umdnj.edu/cocochen/bigwig/$group_name/$bbf_name\n\n";
}

sub load_taxo_hash { #will load %speName2taxID_h 
	#use load_taxo_hash("/Home/cocochen/data/ncbi/taxonomy/ucsc_genos.taxid.txt"), will load a small hash with ucsc genome version map to taxids
	($taxo_f, @species)=@_;
	if (!$taxo_f){$taxo_f="/Home/cocochen/data/ncbi/taxonomy/names.dmp";}
	@headername_arr=("taxid","tmp1","name","tmp2","unique_name","tmp3","name_class");
	open(TAX_F, $taxo_f) || die "error open $taxo_f\n";
	print " open $taxo_f\n";
	while(<TAX_F>){
		s/\s*$//;
		if(/^name/i){ #name	taxid
			@headername_arr=split(/\t/);
			print "header: ".join(" ",@headername_arr)."\n";
		}else{
			@temp_arr = split(/\t/);
			for($i=0;$i<@headername_arr;$i++){
				${$headername_arr[$i]}=$temp_arr[$i];
			}
			$speName2taxID_h{"$name"}=$taxid;
		}
	}
	close TAX_F;
	if(scalar @species>0){
		return map($speName2taxID_h{$_}, @species);
	}
}


sub cal_relative_pos2cdsEnd{
	my($utr3_str,$strand,$pos1,$if_consider_ups_or_intronic)=@_; #utr3_str like 123-456|600-800
	my $ifInUTRorExt=0; #whether in UTR
	my @utr3_poss=split(/\||\-/,$utr3_str);
	my $strand_s=1;
	if($strand eq "-" || $strand eq "-1"){
		@utr3_poss=reverse(@utr3_poss);
		$strand_s=-1;
	}
	$rel_exonic_pos=0;
	if($pos1*$strand_s<$utr3_poss[0]*$strand_s){ #$pos1 upstream of UTR
		if($if_consider_ups_or_intronic){
			$rel_exonic_pos=($pos1-$utr3_poss[0])*$strand_s; 
			return ($rel_exonic_pos);
		}else{
			return 0;
		}
	}
	for(my $i=0;$i<=@utr3_poss-2;$i+=2){
		my $rel_pos=judge_1pos_rel_1reg($utr3_poss[$i],$utr3_poss[$i+1], $pos1);
		if($rel_pos eq "0"){
			$rel_exonic_pos+=abs($pos1-$utr3_poss[$i])+1;
			$ifInUTRorExt=1;
			last;
		}elsif($pos1*$strand_s<$utr3_poss[$i]*$strand_s){ #in upstream intron
			if($if_consider_ups_or_intronic){
				#$rel_exonic_pos+=($pos1-$utr3_poss[$i])*$strand_s;
				last;
			}else{
				return 0;
			}
		}elsif(($i==@utr3_poss-2) && $pos1*$strand_s>=$utr3_poss[$i]*$strand_s){ #terminal exon extention region
			$rel_exonic_pos+=abs($pos1-$utr3_poss[$i])+1;
			$ifInUTRorExt=1;
		}else{
			$rel_exonic_pos+=abs($utr3_poss[$i+1]-$utr3_poss[$i])+1; #not terminal exon
		}
	}
	if(!$ifInUTRorExt){$rel_exonic_pos=0;}
	return ($rel_exonic_pos);
}

sub percentile {
    my ($p,$aref) = @_;
    my $percentile1 = floor($p/100 * $#{$aref});
    my $percentile2 = ceil($p/100 * $#{$aref});
    my @arr_sorted=sort {$a<=>$b} @$aref;
    return ($arr_sorted[$percentile1] + $arr_sorted[$percentile2])/2;
}


sub cal_refseq_genelen{ #input are from ucsc genome browser refflat file (or other similar format)
	my ($cds_start,$cds_end,$exon_starts,$exon_ends)=@_;
	my @exon_start_arr=split(/,/,$exon_starts);
	my @exon_end_arr=split(/,/,$exon_ends);	
	my $cds_len=0; my $transcript_len=0;
	my $if_cds=0;
	for(my $i=0;$i<@exon_start_arr;$i++){ 
		my $exon_start=$exon_start_arr[$i];
		my $exon_end=$exon_end_arr[$i];
		$transcript_len+= abs($exon_end-$exon_start);
		
		my $cds_start2exon= judge_1pos_rel_1reg($exon_start,$exon_end, $cds_start);
		my $cds_end2exon= judge_1pos_rel_1reg($exon_start,$exon_end, $cds_end);
		if($cds_start2exon eq "0"){
			$exon_start=$cds_start;
			$if_cds=1;
		}
		if($cds_end2exon eq "0"){
			$exon_end=$cds_end;
			$cds_len += abs($exon_end-$exon_start);
			$if_cds=0;
		}
		
		if($if_cds){
			$cds_len+= abs($exon_end-$exon_start);
		}
	}	
	return ($cds_len,$transcript_len);
}

sub translate_one_DNA {
	my $seq=@_[0];
	$seq=~s/\w{3}/$& /g;
	$seq=~s/ \w{1,2}$//g;
	$seq=~s/TTT|TTC/F/gi;
	$seq=~s/TAT|TAC/Y/gi;
	$seq=~s/CAT|CAC/H/gi;
	$seq=~s/CAA|CAG/Q/gi;
	$seq=~s/AAT|AAC/N/gi;
	$seq=~s/AAA|AAG/K/gi;
	$seq=~s/GAT|GAC/D/gi;
	$seq=~s/GAA|GAG/E/gi;
	$seq=~s/TGT|TGC/C/gi;
	$seq=~s/TTA|TTG|CTT|CTC|CTA|CTG/L/gi;
	$seq=~s/AGT|AGC|TCT|TCC|TCA|TCG/S/gi;
	$seq=~s/AGA|AGG|CGT|CGC|CGA|CGG/R/gi;
	$seq=~s/ATT|ATC|ATA/I/gi;
	$seq=~s/TGG/W/gi;
	$seq=~s/ATG/M/gi;
	$seq=~s/GTT|GTC|GTA|GTG/V/gi;
	$seq=~s/CCT|CCA|CCC|CCG/P/gi;
	$seq=~s/ACT|ACA|ACC|ACG/T/gi;
	$seq=~s/GCT|GCA|GCC|GCG/A/gi;
	$seq=~s/GGT|GGA|GGC|GGG/G/gi;
	$seq=~s/TAA|TAG|TGA/-/gi; # stop codon to -
	$seq=~s/\w{3}/?/g; #other unconverted
	$seq=~s/ //g;
	return $seq;
}

@all_codons=split(/;/,"TTT;TTC;TAT;TAC;CAT;CAC;CAA;CAG;AAT;AAC;AAA;AAG;GAT;GAC;GAA;GAG;TGT;TGC;TTA;TTG;CTT;CTC;CTA;CTG;AGT;AGC;TCT;TCC;TCA;TCG;AGA;AGG;CGT;CGC;CGA;CGG;ATT;ATC;ATA;TGG;ATG;GTT;GTC;GTA;GTG;CCT;CCA;CCC;CCG;ACT;ACA;ACC;ACG;GCT;GCA;GCC;GCG;GGT;GGA;GGC;GGG;TAA;TAG;TGA;");
sub count_codon_freq{
	my $seq=@_[0];
	$seq=~s/\w{3}/$& /g;
	$seq=uc($seq);
	my @codons=split(/ /, $seq);
	my %codon_fre=();
	map($codon_fre{$_}++, @codons);
	return( map(($codon_fre{$_}?$codon_fre{$_}:0),  @all_codons) );
}


sub cal_2_kmer_mismatch{
	my($seq1, $seq2)=@_;
	my $xor12=$seq1 ^ $seq2;
	$xor12=~s/./ord $& ? "^" : ""/ge;
	return length($xor12);
}


1;


