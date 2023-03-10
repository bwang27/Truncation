#2019/12/2: from .bam file, count read mapped to different chromosome and the desired region of the vector (the region to be packed into the virus)
#2019/12/2: count read number and accumulated length
#2019/12/2: this script is a modified version of 01.Cal_geno_junc_readsnum_frSAM.pl
#2019/12/19: add code to calculate the distribution of read coverage in different regions or chromosomes
#2020/07/02: add function to cluster reads mapped to host genomes
#2021/02/08: add code to also able to deal with paired-end reads
#2022/4: add codes to output read terminals
#2022/06/29-2022/07/14: added codes to output chimeric reads (including IR reads) information using supplementary alignment
#2022/11/01: added code to set clip_cut to a percentage relative to the mapped length
#2022/11/01: not considering mapq, and percent mismatch, clip etc for chimeric read analysis; can also output chimeric read sam file (-C $out_chimeric_read_samf)

use Getopt::Std; 
require "../../../Pipeline/SharedCodes/perl.fun.inc.pl";
require "../../../Pipeline/RNAseq/01.sam.parse.inc.pl";

getopt("nsbpouBtmcvSgC",\%args);

$study_name=$args{n}?$args{n}:"2019-11-11_AAV9-CMV-CDKL5-03";
$sample=$args{s}?$args{s}:"Vg50ulLam25ng"; 
$bam_f=$args{b}?$args{b}:"";
$outPrefix=$args{p}?$args{p}:"$study_name/02_Bam2ReadCounts/";
$out_readnum_f=$args{o}?$args{o}:"$outPrefix$sample.read.counts.txt";
$out_readCov_distr_f="$outPrefix$sample.read.Cov.distr.txt";
$out_readMapCoor_f="$outPrefix$sample.read.map.coor.txt";
$out_readterminal_f="$outPrefix$sample.read.map.terminal.txt";
$out_chimericMapping_f="$outPrefix$sample.read.map.chimeric.txt";
$log_file="$out_readnum_f.log";
create_dir_ifNotExist($log_file);
open (LOGF,">$log_file") || die "error write $log_file\n";

@read_terminal_types=("Plus_start","Plus_end","Minus_start","Minus_end");

$uniqueness_mapq_minSco=$args{u} ne "" ? $args{u} : 10;
$if_best_score=$args{B} ne ""? $args{B} : 1;
$OutSamFrom=$args{t}? $args{t} : "minimap2";
$min_mismatch_in100bp=$args{m} ne ""? $args{m}: 30 ; #minimal number of mismatch for both mates
$clip_cut=$args{c} ne ""? $args{c}: "none"; #maximal number of soft clip
$pairend_gapmax="6000"; # 

$virusGenome_info_f=$args{v}? $args{v} : "target_index/hg19_CDKL5/GOI_info.txt";
$out_highQuality_read_samf=$args{S} ne "" ? $args{S} : ""; #output high quality reads mapping to a sam file
$out_chimeric_read_samf=$out_highQuality_read_samf;
$out_chimeric_read_samf=~s/highQuality/chimeric/i;
$out_chimeric_read_samf=$args{C} if  $args{C} ne "" ;

$debug=0;
my %counts;


##2, read virusGenome_info_f, load  %geno_info_h
open (VIRU_GENO_INFO, $virusGenome_info_f) || die "error open $virusGenome_info_f\n";
print "#opened $virusGenome_info_f\n";
print LOGF "#opened $virusGenome_info_f\n";
while(<VIRU_GENO_INFO>){
	s/\s+$//;
	if(/^vector|^chromosome/){ #vector	start	end	name #chromosome	from	to	name
		@headername_arr=split(/\t/);
		print $_."\n";
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		$geno_info_h{$chromosome}{"$from	$to"}=$name;
		print "\$geno_info_h{$chromosome}{$from	$to}=$name;\n";
	}

}
close VIRU_GENO_INFO;

##3a, read bam file and store all read IDs with supplementary alighments ($flag & 0x800)
my %supp_alignNum_h=();
my $openSamCmd;
	
if($bam_f=~/.bam$/){
	$openSamCmd="samtools view $bam_f | ";
}else{
	$openSamCmd=" $bam_f "
}
print "#open \$bam_f=$openSamCmd\n";
open(SAM_F,$openSamCmd) || die "error open $openSamCmd\n";
while(<SAM_F>){
	next if /^@/; #comment line: @HD     VN:1.0  SO:sorted
	s/\s+$//;
	($qname_raw,$flag,$rname,$posi,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,@tags_arr)=split(/\t/);
	$qname=$qname_raw;
	$qname=~s/\.[12]$//; 
	if($flag & 0x800){
		$supp_alignNum_h{$qname}++;
	}
}
close SAM_F;



##3, read bam file, filter reads, count reads mapped to the virus genome

print "#open \$bam_f=$openSamCmd\n";
print LOGF "#open \$bam_f=$openSamCmd\n";
open(SAM_F,$openSamCmd) || die "error open $openSamCmd\n";
open (RNAME_RNUM_OUT,">$out_readnum_f") || die "error write $out_readnum_f\n";
if($out_highQuality_read_samf){
	create_dir_ifNotExist($out_highQuality_read_samf);
	open(OUT_SAM,">$out_highQuality_read_samf") || die "error write $out_highQuality_read_samf\n";
}
if($out_chimeric_read_samf){
	create_dir_ifNotExist($out_chimeric_read_samf);
	open(OUT_SAM2,">$out_chimeric_read_samf") || die "error write $out_chimeric_read_samf\n";
}

my %rname2rnum_h=();
my %ReadCov2region2Count_h=();
my %readMapCoor_h=(); #store read mapping coordinates (reference genomic location)
my %readTerminal_h=(); #store read terminal read counts information
my %chimeric_reads_info_h=(); #store chimeric reads (including IR reads) information

$rowi=0;
while(<SAM_F>){
	next if /^@/; #comment line: @HD     VN:1.0  SO:sorted
	s/\s+$//;
	($qname_raw,$flag,$rname,$posi,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,@tags_arr)=split(/\t/);
	$qname=$qname_raw;
	$qname=~s/\.[12]$//; 
	$rowi++;
	if($rowi % 1000000==0){
		print  return_time()."|read SAM file $sample line $rowi\n";
	}
	if(eof){
		print  return_time()."|read SAM file $sample line $rowi\n";
		print LOGF "#".return_time()."|read SAM file $sample line $rowi\n";
	}
	
	last if($rowi>100000 && $debug); 
	$counts{$sample}{"0.reads.row_inSam"}++;
	if($flag & 0x0004){next;} #not mapped
	
	my $read_map_primary_orNot;
	if($flag & 0x0100){
		$read_map_primary_orNot="secondaryAli";
	}else{
		$read_map_primary_orNot="primaryAli";
	}
	next if ($read_map_primary_orNot ne "primaryAli");
	$counts{$sample}{"1.reads.$read_map_primary_orNot"}++;

	my $reads_strand=from_flag2strand($flag);
	if($flag & 0x0080){ #reverse the strand for mate2
		$reads_strand=$reads_strand=~/-/?"+":"-";
	}
	my $read_cover_len=cigar2cover_len($cigar); #for paired-end reads, this value will be the coverage of one read instead of one pair of reads
	my $readPair_Cov_len=$read_cover_len; #for single end reads, $readPair_Cov_len=$read_cover_len (for code compatibility with paired-end reads)
	my $posto=$posi+$read_cover_len-1;

	my $mate1or2=($flag & 0x0080)?"mate2":"mate1";
	my $tag_str=join("	",@tags_arr);
	my $aliSco="";
	if($tag_str=~/AS:i:(\d+)/i){$aliSco=$1;}

	my $mismatch=0; #number of mismatch
	if($tag_str=~/NM:i:(\d+)/i ){$mismatch=$1;} #only consider mismatch for single end reads

	##save chimeric reads information (not considering mapq, and percent mismatch, clip etc for chimeric read analysis)
	if($supp_alignNum_h{$qname}>0){ 
		my $readLen=cigar2read_len($cigar);
		my ($leftClipLen,$readEndCoor)=cigar2readMapCorr($cigar);
		if($reads_strand eq "-"){ #if a read map to - strand, the cigar is based on the reverse complementary the original read sequence. In this case, convert the coordinates to original strand
			($leftClipLen,$readEndCoor)=($readLen-$readEndCoor, $readLen-$leftClipLen);
		}
		$chimeric_reads_info_h{$qname}{$leftClipLen}= ($leftClipLen+1).":$readEndCoor:".($readEndCoor-$leftClipLen).":$reads_strand:$rname:$posi:$posto:$mapq:$read_cover_len:$aliSco:$mismatch";
		print OUT_SAM2 "$_\n" if $out_chimeric_read_samf;
	}


	#uniqueness
	if($uniqueness_mapq_minSco && $mapq < $uniqueness_mapq_minSco){ #judge uniqueness by $mapq
		next;
	}
	$counts{$sample}{"2.reads.$read_map_primary_orNot.MAPQ>=$uniqueness_mapq_minSco"}++;
	
	if($if_best_score){
		if($tag_str=~/XS:i:(\d+)/i){
			$xs=$1;
			next if $aliSco<$xs;
		}
		$counts{$sample}{"4.reads.$read_map_primary_orNot.BestAliSco"}++;
	}


	my $mismatch_cut=$min_mismatch_in100bp*($read_cover_len/100);
	my $clipcut_nt=$clip_cut=~/%$/ ? $read_cover_len*$clip_cut : $clip_cut;
	if(! judge_reads_quality($OutSamFrom,$mismatch_cut,$clipcut_nt,$cigar,@tags_arr)){
		next;
	}
	$counts{$sample}{"5.reads.$read_map_primary_orNot.mismatch<=$min_mismatch_in100bp/100bp.clip<=$clip_cut.$mate1or2"}++;
	
	my $simple_cigar=modify_cigar_rm_I_D($cigar);
	my $readmap_info="$rname $posi $simple_cigar $posto";
	my $if_duplicate_readmap=$allReadMap_info_h{$readmap_info} ? 1 : 0; #whether the same read mapping (different UMI is considered different mapping) has happened for a previous read
	$allReadMap_info_h{$readmap_info}=1;


	my ($reads_cigarcode,$readmap_type, $readPair_start, $readPair_end);

	if($flag & 0x1){ #paired-end
		if(!($flag & 0x0002) || ($flag & 0x0008)){next;} #if not mapped in a proper pair
		$counts{$sample}{"6.reads.$read_map_primary_orNot.properMapped"}++;
		$readPair_Cov_len=abs($isize);
		my $read1_strand=($flag & 0x10) ? "-" : "+";
		if($read1_strand eq "+" && $isize>0){
			($readPair_start, $readPair_end)=($posi, $posi+$isize-1);
			if ($mate1or2 eq "mate1"){
				$readMapCoor_h{$rname}{"+	$readPair_start	$readPair_end"}++ ;
				$readTerminal_h{$rname}{$readPair_start}{"Plus_start"}++;
				$readTerminal_h{$rname}{$readPair_end}{"Plus_end"}++;
			}
			}elsif($read1_strand eq "-" && $isize<0){
			($readPair_start, $readPair_end)=($mpos, $mpos-$isize-1); 
			if ($mate1or2 eq "mate1"){
				$readMapCoor_h{$rname}{"-	$readPair_start	$readPair_end"}++ ;
				$readTerminal_h{$rname}{$readPair_end}{"Minus_start"}++;
				$readTerminal_h{$rname}{$readPair_start}{"Minus_end"}++;
			}
			
		}elsif($read1_strand eq "+" && $isize<0){ ##rare situation, read 1 map to + strand at the end of a circular plasmid
			($readPair_start, $readPair_end)=($posi, $posto);
			$readPair_Cov_len=""; #read-pair coverage not easy to calculate 
		}elsif($read1_strand eq "-" && $isize>0){  ##rare situation, read 1 map to - strand at the start of a circular plasmid
			($readPair_start, $readPair_end)=($posi, $posto);
			$readPair_Cov_len=""; #read-pair coverage not easy to calculate 
		}
		$readmap_type="exon";

	}else{#single end
		($reads_cigarcode,$readmap_type)=judge_readstype_fr_ciger($cigar);
		($readPair_start, $readPair_end)=($posi, $posto);
		$readMapCoor_h{$rname}{"$reads_strand	$posi	$posto"}++;
		if($reads_strand eq "+"){
			$readTerminal_h{$rname}{$posi}{"Plus_start"}++;
			$readTerminal_h{$rname}{$posto}{"Plus_end"}++;
		}else{
			$readTerminal_h{$rname}{$posto}{"Minus_start"}++;
			$readTerminal_h{$rname}{$posi}{"Minus_end"}++;
		}

	}

	print OUT_SAM "$_\n" if $out_highQuality_read_samf;
	$counts{$sample}{"7.reads.Used.$mate1or2.$readmap_type"}++;
	$counts{$sample}{"8.reads.collapsed.$mate1or2.$readmap_type"}++ if (!$if_duplicate_readmap);
	

	if($readmap_type eq "exon"){
		# $rname	$reads_strand	$posi	$posto
		my @annotated_regions=keys %{$geno_info_h{$rname}};
		if(scalar @annotated_regions >0){ # at least one defined region in this vector to look at
			my $not_on_target_align_len=$read_cover_len;
			foreach my $region_str(@annotated_regions){
				my $target_region_name=$geno_info_h{$rname}{$region_str};
				my ($reg_start, $reg_end)=split(/\t/, $region_str);
				my ($read_rel2_anoReg, $overlap_len)=judge_2region_relation($reg_start, $reg_end, $posi,$posto,"+", 1); # read relative position to annotated genomic region
				if($read_rel2_anoReg=~/0|=/ || ($read_rel2_anoReg=~/0|=|ov|inc/ && $target_region_name!~/GOI/ ) ){ #read inside region or read inside or overlap for non-GOI gene features
					$rname2rnum_h{"$rname	$target_region_name"}{"counts"}++;
					$rname2rnum_h{"$rname	$target_region_name"}{"length"}+=$read_cover_len;
					$rname2rnum_h{"$rname	$target_region_name"}{"mismatch"}+=$mismatch;
					$ReadCov2region2Count_h{$readPair_Cov_len}{"$rname	$target_region_name"}++ if $readPair_Cov_len;
					$not_on_target_align_len=0;
				}elsif($read_rel2_anoReg=~/ov/ ){ #ov1 or ov-1
					# $rname2rnum_h{"$rname	$target_region_name"}{"length"}+=$overlap_len; #only add length, not add count
					$not_on_target_align_len-=$overlap_len;
					my $ovlap_key="$rname	$target_region_name".($read_rel2_anoReg eq "ov-1"?".OvLeft":".OvRight");
					$rname2rnum_h{ $ovlap_key }{"counts"}++;
					$rname2rnum_h{ $ovlap_key }{"length"}+=$read_cover_len;
					$rname2rnum_h{ $ovlap_key }{"mismatch"}+=$mismatch;
					$ReadCov2region2Count_h{$readPair_Cov_len}{ $ovlap_key }++ if $readPair_Cov_len;
				}elsif($read_rel2_anoReg=~/inc/){ #read include the whole region 
					#$rname2rnum_h{"$rname	$target_region_name"}{"length"}+=$overlap_len;
					$not_on_target_align_len-=$overlap_len;
					foreach my $OvType ("OvLeft", "OvRight"){
						$rname2rnum_h{"$rname	$target_region_name.$OvType"}{"counts"}+=1;
						$rname2rnum_h{"$rname	$target_region_name.$OvType"}{"length"}+=$read_cover_len;
						$rname2rnum_h{"$rname	$target_region_name.$OvType"}{"mismatch"}+=$mismatch;
						$ReadCov2region2Count_h{$readPair_Cov_len}{ "$rname	$target_region_name.$OvType" }++ if $readPair_Cov_len;
					}

				}elsif($read_rel2_anoReg=~/1/){ #1 or -1
					##do nothing
				}
			}
			if($not_on_target_align_len>0){
				#$rname2rnum_h{"$rname	Others"}{"length"}+=$not_on_target_align_len;
				if($not_on_target_align_len==$read_cover_len){
					$rname2rnum_h{"$rname	Others"}{"counts"}++;
					$rname2rnum_h{"$rname	Others"}{"length"}+=$read_cover_len;
					$rname2rnum_h{"$rname	Others"}{"mismatch"}+=$mismatch;
					$ReadCov2region2Count_h{$readPair_Cov_len}{ "$rname	Others" }++ if $readPair_Cov_len;
				}
			}
		}
		foreach my $rname1 ($rname, "All"){
			$rname2rnum_h{"$rname1	All"}{"counts"}++;
			$rname2rnum_h{"$rname1	All"}{"length"}+=$read_cover_len;
			$rname2rnum_h{"$rname1	All"}{"mismatch"}+=$mismatch;
			$ReadCov2region2Count_h{$readPair_Cov_len}{ "$rname1	All" }++ if $readPair_Cov_len;
		}
	}


}
close SAM_F;
close OUT_SAM if $out_highQuality_read_samf;
close OUT_SAM2 if $out_chimeric_read_samf;

##write RNAME_RNUM_OUT
print RNAME_RNUM_OUT "reference	region	Num_$sample	Len_$sample	Mismatch_$sample\n";
foreach my $rname_reg (sort keys %rname2rnum_h){
	print RNAME_RNUM_OUT "$rname_reg	". 
		join("	", map($rname2rnum_h{$rname_reg}{$_}, "counts","length","mismatch" ) ) ."\n";
}
close RNAME_RNUM_OUT;

##output $out_readMapCoor_f

open (READ_MAP_COOR,">$out_readMapCoor_f") || die "error write $out_readMapCoor_f\n";
print READ_MAP_COOR "rname	reads_strand	posfr	posto	readsnum\n";
foreach my $ref(sort keys %readMapCoor_h){
	print READ_MAP_COOR join("", map("$ref	$_	$readMapCoor_h{$ref}{$_}\n", sort keys %{$readMapCoor_h{$ref}}));
}
close READ_MAP_COOR;


##output %readTerminal_h
open (READ_MAP_TERMINAL,">$out_readterminal_f") || die "error write $out_readterminal_f\n";
print READ_MAP_TERMINAL "rname	position	".join("	",map("Num_$_",@read_terminal_types))."\n";
foreach my $ref(sort keys %readTerminal_h){
	foreach my $coor ( sort {$a<=>$b}  keys %{$readTerminal_h{$ref}} ){
		print READ_MAP_TERMINAL "$ref	$coor	".join("	",map($readTerminal_h{$ref}{$coor}{$_} ? $readTerminal_h{$ref}{$coor}{$_} :0, @read_terminal_types))."\n";
	}
}
close READ_MAP_TERMINAL;

#output chimeric_reads_info_h #$chimeric_reads_info_h{$qname}{$leftClipLen}= ($leftClipLen+1).":$readEndCoor:".($readEndCoor-$leftClipLen).":$reads_strand:$rname:$posi:$posto:$mapq:$read_cover_len:$aliSco:$mismatch";
open (READ_CHIMERIC,">$out_chimericMapping_f") || die "error write $out_chimericMapping_f\n";
print READ_CHIMERIC "ReadID	ChimericMappingInfo\n"; #ReadStart:ReadEnd:ReadStrand:TargetName:TargetStart:TargetEnd
foreach my $qname(sort keys %chimeric_reads_info_h){
	if(scalar keys %{$chimeric_reads_info_h{$qname}}>1){
		print READ_CHIMERIC "$qname	".join("|", map($chimeric_reads_info_h{$qname}{$_}, sort {$a<=>$b} keys %{$chimeric_reads_info_h{$qname}}))."\n";
		$counts{$sample}{"9a.chimeric.mapping"}+=scalar keys %{$chimeric_reads_info_h{$qname}};
		$counts{$sample}{"9b.reads.withChimericMapping"}++;
	}
}
close READ_CHIMERIC;



####output numbers:
my %all_readGrps_h=();
foreach my $sample1( keys %counts){
	map($all_readGrps_h{$_}=1, keys %{$counts{$sample1}});
}
my @all_readGrps=sort keys %all_readGrps_h;
$print_txt="Sample	".join("	",@all_readGrps)."\n";
foreach my $sample1( keys %counts){
	$print_txt.="$sample1	". join("	",map($counts{$sample1}{$_}, @all_readGrps))."\n";
}

print $print_txt;
print LOGF $print_txt;


##output read coverage distribution matrix
my @all_readCovLengths=sort {$a<=>$b} keys %ReadCov2region2Count_h;
open (READCOV_MATRIX_OUT,">$out_readCov_distr_f") || die "error write $out_readCov_distr_f\n";
my @all_ref_names=sort keys %rname2rnum_h;
print READCOV_MATRIX_OUT "CovLength";
foreach my $region1(@all_ref_names){
	my $region1_txt=$region1;
	$region1_txt=~s/\t/./;
	print READCOV_MATRIX_OUT "\t$region1_txt";
}
print READCOV_MATRIX_OUT "\n";

print join(" ",@all_ref_names)."\n";
foreach my $covLen (@all_readCovLengths){
	print READCOV_MATRIX_OUT "$covLen	".join("\t", map($ReadCov2region2Count_h{$covLen}{$_}, @all_ref_names))."\n";
}
close READCOV_MATRIX_OUT;

close LOGF;
