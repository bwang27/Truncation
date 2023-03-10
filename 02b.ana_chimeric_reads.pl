##read the output of 02.count_AAV_read_fromBam.pl (.read.map.chimeric.txt files) and do futher classification of chimeric reads to several groups: IR reads, chimeric reads, reads mapped to plasmid bundary (not real chimeric), etc
##2022/08/09 added function to output bed format files for different types of chimeric reads
##2022/08/11 added function to output BEDPE format files for visualization
##2022/08/17, deleted/commented codes to output bed format
##2022/09/21, added function to extract the DNA template sequences near the chimeric events fusion point (need extract_seq_fr_contig function from perl.fun.inc.pl)
##2022/11/01, changed some default parameters to call a chimeric reads: $read_block_distance_max=50; $min_mapq_cutoff=0; 

require "../../../Pipeline/SharedCodes/perl.fun.inc.pl";
use Getopt::Std; 
getopt("nspioSgQDIMF",\%args);

$study_name=$args{n}?$args{n}:"2019-11-11_AAV9-CMV-CDKL5-03";
$sample=$args{s}?$args{s}:"Vg50ulLam25ng";
$in_chimeric_read_f=$args{i}?$args{i}:"$study_name/02_Bam2ReadCounts/$sample.read.map.chimeric.txt";
$outPrefix=$args{p}?$args{p}:"$study_name/02b_ChimericRead/";
$out_chimericRead_f=$args{o}?$args{o}:"$outPrefix$sample.chimericRead.detail.txt"; #output each individual chimeric read information
$out_chimericRead_stats_f="$out_chimericRead_f.stats.txt"; #overall stats of counts for different categories
$out_chimericRead_count_f="$outPrefix$sample.chimericReadCount.txt"; #bin high quality chimeric reads based on joint point coordinates and output their counts
$log_file="$out_chimericRead_f.log";

$bed_thin_size=50; #number of NTs shows as a thin line in the bed format visualization
$read_block_distance_max=$args{D}?$args{D}:50; #maximal distance between the read end coordinates for block1 and read start coordinates for block 2
$short_indel_size_max=$args{I}?$args{I}:5; #the distance of the two alignment blockes above this value will be considered as a chimeric reads (not small indel)
$plasmid_bundary_tolerance_nt=8; #the maximal distance between the mapping position and the plasmid bundary to be deemd as "mapping to bundary"

$min_mapq_cutoff=$args{Q} ne ""?$args{Q}:50; #Minimal MAPQ score of the two alignments
$max_percMM_cutoff=$args{M} ne ""?$args{M}:0.15; #maximal percent mismatch of the two alignments

create_dir_ifNotExist($log_file);
open (LOGF,">$log_file") || die "error write $log_file\n";

my %strand2rgb=("+"=>"255,0,0", "-"=>"0,0,255");


$genome_fa_f=$args{F}? $args{F} : ""; #$index_folder/$index_name/$index_name.fa
$ref_sizes_f=$args{S}? $args{S} : "$genome_fa_f.sizes"; #$index_folder/$index_name/$index_name.fa.sizes
$GOI_name=$args{g}? $args{g} : "";
$extr_seq_inside_size=30; #sequence length extracted from inside of each fusion point
$extr_seq_outside_size=19; #sequence length extracted from outside of each fusion point, the actual extract sequence will be 1 nt longer than this setting

sub build_extract_seq_regions {
	my ($refName1,$refName2,$strand1,$strand2,$refEnd1,$refStart2, $extr_seq_inside_size, $extr_seq_outside_size)=@_;
	my $strand_sign1=$strand1 eq "-"?1:-1; #extract the antisense strand sequence of the read strand
	my $strand_sign2=$strand2 eq "-"?1:-1;
	my @extract_regions=();
	if($refName1 ne $refName2){
		@extract_regions=(
			$refName2."|". ($strand2 eq "-"?"+":"-")."|".  ($refStart2-$extr_seq_outside_size*$strand_sign2).":".($refStart2),
			$refName2."|". ($strand2 eq "-"?"+":"-")."|".  ($refStart2+$strand_sign2).":".($refStart2+$extr_seq_inside_size*$strand_sign2),
			$refName1."|". ($strand1 eq "-"?"+":"-")."|".  ($refEnd1-$extr_seq_inside_size*$strand_sign1).":".($refEnd1-$strand_sign1) , 
			$refName1."|". ($strand1 eq "-"?"+":"-")."|".  ($refEnd1).":".($refEnd1+$extr_seq_outside_size*$strand_sign1)  
		);
	}else{ #same reference
		my ($leftcorr,$rightcorr)=sort {$a<=>$b} ($refEnd1, $refStart2);
		if(abs($refEnd1-$refStart2)+1<=2*$extr_seq_inside_size){
			@extract_regions=(
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($leftcorr-$extr_seq_outside_size).":".($leftcorr), 
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($leftcorr+1).":".($rightcorr-1),
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($rightcorr).":".($rightcorr+$extr_seq_outside_size)  
			); #middle region has small size, merged
		}else{
			@extract_regions=(
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($leftcorr-$extr_seq_outside_size).":".($leftcorr),
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($leftcorr+1).":".($leftcorr+$extr_seq_inside_size),
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($rightcorr-$extr_seq_inside_size).":".($rightcorr-1),
				$refName1 ."|". ($strand1 eq "-"?"+":"-") ."|". ($rightcorr).":".($rightcorr+$extr_seq_outside_size)
			);
		}
		if($strand1 eq "+"){
			@extract_regions=reverse(@extract_regions);
		}
	}
	return(@extract_regions);
}

##1, read ref_sizes_f, load  %ref_sizes_h
die "Error: ref_sizes_f not defined!" if !$ref_sizes_f;
open (REF_SIZES, $ref_sizes_f) || die "error open $ref_sizes_f\n";
print "#opened $ref_sizes_f\n";
print LOGF "#opened $ref_sizes_f\n";
while(<REF_SIZES>){
	s/\s+$//;
	@temp_arr = split(/\t/);
	$ref_sizes_h{$temp_arr[0]}=$temp_arr[1];
	print "\$ref_sizes_h{".$temp_arr[0]."}=".$temp_arr[1]."\n";
}
close REF_SIZES;
print $ref_sizes_h{"ID23_scAAV-SCA3-4E10x2"}."\n";

##2, read $in_chimeric_read_f and process 
die "Error: ref_sizes_f not defined!" if !$in_chimeric_read_f;
open (IN_CHIMERIC, $in_chimeric_read_f) || die "error open $in_chimeric_read_f\n";
print LOGF "#opened $in_chimeric_read_f\n";
print "#opened $in_chimeric_read_f\n";
open (OUT_CHIMERIC,">$out_chimericRead_f") || die "error write $out_chimericRead_f\n";
print OUT_CHIMERIC "ReadID	chimeric_Ref_type	chimeric_type	refName1	refName2	ReadMapBlock1	ReadMapBlock2	min_mapq	min_refCovLen	min_ali_sco	max_percMM	twoBlockDist	ifHighQuality\n";
print LOGF "#write $out_chimericRead_f\n";
print "#write $out_chimericRead_f\n";
my %chimeric_stats_h=();
my %bed_lines_h=();
my %chimericInfo2counts_h=();
while(<IN_CHIMERIC>){
	s/\s+$//;
	next if !$_;
	if(/^ReadID|ChimericMappingInfo/){ #ReadID	ChimericMappingInfo
		@headername_arr=split(/\t/);
		print $_."\n";
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		my @mappingInfo_arr=split(/\|/, $ChimericMappingInfo);
		my $mapping_block_num=scalar @mappingInfo_arr;
		next if $mapping_block_num<2;
		foreach my $i( 0..($mapping_block_num-2)){
			foreach my $j( ($i+1)..($mapping_block_num-1) ){
				my ($readStart1,$readEnd1,$readCovLen1,$strand1,$refName1,$refStart1,$refEnd1,$mapq1,$refCovLen1,$aliSco1,$mismatch1)=split(/:/, $mappingInfo_arr[$i]); # ($leftClipLen+1).":$readEndCoor:".($readEndCoor-$leftClipLen).":$reads_strand:$rname:$posi:$posto:$mapq:$read_cover_len:$aliSco:$mismatch
				my ($readStart2,$readEnd2,$readCovLen2,$strand2,$refName2,$refStart2,$refEnd2,$mapq2,$refCovLen2,$aliSco2,$mismatch2)=split(/:/, $mappingInfo_arr[$j]);
				next if abs($readStart2-$readEnd1)>$read_block_distance_max; #the coordinates in read didn't connect
				if($strand1 eq "-"){($refStart1,$refEnd1)=($refEnd1,$refStart1)} #if - strand, start position > end position
				if($strand2 eq "-"){($refStart2,$refEnd2)=($refEnd2,$refStart2)}
				next if ($refName1 eq $refName2) && ($strand1 eq $strand2) && abs($refStart2-$refEnd1)<=$short_indel_size_max; ##may be short indels
				my $twoBlockDist=-1;
				my $refType1=$refName1 eq $GOI_name ? "GOI":"others";
				my $refType2=$refName2 eq $GOI_name ? "GOI":"others";
				my $chimeric_Ref_type=$refName1 eq $refName2 ? "intra-":"inter-";
				$chimeric_Ref_type=$chimeric_Ref_type.join("-", sort($refType1,$refType2) );
				my $chimeric_type="fusion"; #fusion between different reference sequences 
				if($refName1 eq $refName2){ #block 1 and 2 in same reference, define details of chimeric_type
					$twoBlockDist=abs($refStart2-$refEnd1); #the distance between block1_end and block2_start when they map to the same reference
					if($strand1 ne $strand2){
						$chimeric_type="Inversion$strand1";
					}else{#same strand
						my $strand_sign=$strand1 eq "+"?1:-1;
						if($refStart2*$strand_sign>$refEnd1*$strand_sign){ #block 2 downstream of block 1 #   ----1---> ----2---> 
							$chimeric_type="Deletion$strand1";
						}elsif($refStart2*$strand_sign>=$refStart1*$strand_sign || $refEnd2*$strand_sign>$refStart1*$strand_sign ){ #block 2 start in the middle of block 1 ( -----1-=====>---2----> ) or block 2 end in the middle of block 1 ( -----2-=====>---1----> )
							$chimeric_type="Insertion$strand1";
						}elsif($refStart2*$strand_sign<$refStart1*$strand_sign){ #block 2 start upstream of block 1:   ----2--->   ----1--->
							if($strand1 eq "+" && $refStart2<=$plasmid_bundary_tolerance_nt && abs($ref_sizes_h{$refName1}-$refEnd1)<=$plasmid_bundary_tolerance_nt){ # |----2--->      ----1--->|
								$chimeric_type="Boundary+";
							}elsif($strand1 eq "-" && $refEnd1<=$plasmid_bundary_tolerance_nt && abs($ref_sizes_h{$refName1}-$refStart2)<=$plasmid_bundary_tolerance_nt){ # |<---1---    <----2----|
								$chimeric_type="Boundary-";
							}else{
								$chimeric_type="Reversion$strand1";
							}
						}
					}
				}else{
					#fusion point close to bundary of plasmid, potential mapping artifact
					if($refEnd1<=$plasmid_bundary_tolerance_nt || $refStart2<=$plasmid_bundary_tolerance_nt || abs($ref_sizes_h{$refName2}-$refStart2)<=$plasmid_bundary_tolerance_nt ||  abs($ref_sizes_h{$refName1}-$refEnd1)<=$plasmid_bundary_tolerance_nt ){
						$chimeric_type="FuseToBound";
					}
				}
				#output a pair of alignment blocks
				my $min_mapq=$mapq1<$mapq2?$mapq1:$mapq2;
				my $min_refCovLen=$refCovLen1<$refCovLen2?$refCovLen1:$refCovLen2;
				my $min_ali_sco=$aliSco1<$aliSco2?$aliSco1:$aliSco2;
				my $max_percMM=($mismatch1/$refCovLen1>$mismatch2/$refCovLen2) ? ($mismatch1/$refCovLen1) : ($mismatch2/$refCovLen2);
				my $ifHighQuality=($min_mapq>=$min_mapq_cutoff && $max_percMM<=$max_percMM_cutoff)? "Y":"N";
				$chimeric_stats_h{"$chimeric_Ref_type	$chimeric_type	$ifHighQuality"}++;
				$chimericInfo2counts_h{"$chimeric_Ref_type	$chimeric_type	$refName1	$refName2	$strand1	$strand2	$refEnd1	$refStart2"}++ if $ifHighQuality eq "Y";
				print OUT_CHIMERIC "$ReadID	$chimeric_Ref_type	$chimeric_type	$refName1	$refName2	$readStart1:$readEnd1:$readCovLen1:$strand1:$refStart1:$refEnd1:$refCovLen1	".
					"$readStart2:$readEnd2:$readCovLen2:$strand2:$refStart2:$refEnd2:$refCovLen2	$min_mapq	$min_refCovLen	$min_ali_sco	$max_percMM	$twoBlockDist	$ifHighQuality\n";
			}
		}
	}
}
close IN_CHIMERIC;
close OUT_CHIMERIC;

##output %chimeric_stats_h to out_chimericRead_stats_f
open (OUT_CHIMERIC_STATS,">$out_chimericRead_stats_f") || die "error write $out_chimericRead_stats_f\n";
print LOGF "#write stats: $out_chimericRead_stats_f\n";
print "#write stats: $out_chimericRead_stats_f\n";
print OUT_CHIMERIC_STATS "chimeric_Ref_type	chimeric_type	ifHighQuality_MAPQ$min_mapq_cutoff.MM$max_percMM_cutoff	NumberOfMapping\n";
foreach my $chimeric_type (sort keys %chimeric_stats_h){
	print OUT_CHIMERIC_STATS "$chimeric_type	$chimeric_stats_h{$chimeric_type}\n";
}
close OUT_CHIMERIC_STATS;

##output %chimericInfo2counts_h to out_chimericRead_count_f; also extract the sequence near the fusion point, create  %outBEDPE_row_h, 
my %outbedrows_h=();
my %outBEDPE_row_h=();
open (OUT_CHIMERIC_COUNTS,">$out_chimericRead_count_f") || die "error write $out_chimericRead_count_f\n";
print LOGF "#write stats: $out_chimericRead_count_f\n";
print "#write stats: $out_chimericRead_count_f\n";
print OUT_CHIMERIC_COUNTS "chimeric_Ref_type	chimeric_type	refName1	refName2	strand1	strand2	refEnd1	refStart2	count	template_seq\n";
foreach my $chimeric_key (sort { $chimericInfo2counts_h{$b} <=> $chimericInfo2counts_h{$a} } keys %chimericInfo2counts_h){
	my $read_count=$chimericInfo2counts_h{$chimeric_key};
	my ($chimeric_Ref_type,$chimeric_type,$refName1,$refName2,$strand1,$strand2,$refEnd1,$refStart2)=split(/\t/,$chimeric_key);
	my $template_seq="";
	if($chimeric_type!~/Bound/){
		my @out_regions=build_extract_seq_regions($refName1,$refName2,$strand1,$strand2,$refEnd1,$refStart2, $extr_seq_inside_size, $extr_seq_outside_size);
		my @out_seq_arr=(); #four or three sub-regions
		foreach my $i(0..(scalar @out_regions-1) ){
			my $out_region1=$out_regions[$i];
			my ($chro,$strand,$region1)=split(/\|/,$out_region1);
			my ($region1Start,$region1End)=sort {$a<=>$b} split(/:/,$region1);
			$region1Start=1 if $region1Start<=0;
			my $out_seq1=extract_seq_fr_contig("", $chro, $strand, "$region1Start:$region1End", $genome_fa_f);
			if($i==0 || $i==scalar @out_regions-1){
				$out_seq1=lc($out_seq1);
			}
			push (@out_seq_arr, $out_seq1);
		}
		$template_seq=join("|",@out_seq_arr);
	}
	print OUT_CHIMERIC_COUNTS "$chimeric_key	$read_count	$template_seq\n";
	my $feature_name="$refName1:$strand1:${refEnd1}~".($refName1 eq $refName2?"":"$refName2:$strand2:")."$refStart2:$read_count";
	# my $bed_row=convert_coor_to_bedRow($chimeric_type, $refName1,$refName2,$strand1,$strand2, $refEnd1, $refStart2, $bed_thin_size, $read_count);
	# my $score=$read_count*($strand1 eq "+"?1:-1);
	# $chimeric_type2=$chimeric_type;
	# $chimeric_type2=~s/[+-]//;
	# if($bed_row){
	# 	$outbedrows_h{$chimeric_type2}.="$bed_row\n";
	# }
	my $bedpe_row="$refName1	".($refEnd1-1)."	$refEnd1	$refName2	".($refStart2-1)."	$refStart2	$feature_name	$read_count	$strand1	$strand2	$chimeric_Ref_type	$chimeric_type";
	if($chimeric_type=~/Inversion|Insertion|Deletion|Reversion|fusion/){
		$outBEDPE_row_h{$chimeric_type}.="$bedpe_row\n";
	}
}
close OUT_CHIMERIC_COUNTS;

#output bed format files
# foreach my $chimeric_type(keys %outbedrows_h){
# 	my $out_bed_f="${outPrefix}/visu/$sample.$chimeric_type.bed";
# 	create_dir_ifNotExist($out_bed_f);
# 	open (OUTBED, ">$out_bed_f") || die "error write $out_bed_f\n";
# 	print OUTBED $outbedrows_h{$chimeric_type};
# 	close OUTBED;
# }
foreach my $chimeric_type(keys %outBEDPE_row_h){
	my $out_bedpe_f="${outPrefix}/visu/$sample.$chimeric_type.bedpe";
	create_dir_ifNotExist($out_bedpe_f);
	open (OUTBEDPE, ">$out_bedpe_f") || die "error write $out_bedpe_f\n";
	print OUTBEDPE $outBEDPE_row_h{$chimeric_type};
	close OUTBEDPE;
}
close LOGF;


sub convert_coor_to_bedRow { #convert chimeric read information to a bed format row
	my ($chimeric_type, $refName1,$refName2,$strand1,$strand2, $refEnd1, $refStart2, $bed_thin_size, $read_count)=@_;
	my $bed_row="";
	my $bed_start; my $bed_end;
	if($chimeric_type=~/Inversion/){
		my @wide2coors=sort {$a <=> $b} ($refEnd1, $refStart2);
		if($strand1 eq "+"){
			$bed_start=$wide2coors[0]-$bed_thin_size-1; #0-based 
			$bed_start=0 if $bed_start<0;
			$bed_row="$refName1	$bed_start	$wide2coors[1]	$refName1:$strand1:$wide2coors[0]-$wide2coors[1]:$read_count	$read_count	$strand1	".($wide2coors[0]-1)."	$wide2coors[1]	$strand2rgb{$strand1}";
		}else{ # - strand
			$bed_end=$wide2coors[1]+$bed_thin_size;
			$bed_row="$refName1	".($wide2coors[0]-1)."	$bed_end	$refName1:$strand1:$wide2coors[0]-$wide2coors[1]:$read_count	$read_count	$strand1	".($wide2coors[0]-1)."	$wide2coors[1]	$strand2rgb{$strand1}";
		}
	}elsif($chimeric_type=~/Insertion/){
		my @wide2coors=sort {$a <=> $b} ($refEnd1, $refStart2);
		my @all_4coor=($wide2coors[0]-$bed_thin_size, @wide2coors, $wide2coors[1]+$bed_thin_size); #all 1-based coordinates
		$all_4coor[0]=1 if $all_4coor[0]<1;
		$bed_row="$refName1	".($all_4coor[0]-1)."	$all_4coor[3]	$refName1:$strand1:$wide2coors[0]-$wide2coors[1]:$read_count	$read_count	$strand1	".($wide2coors[0]-1)."	$wide2coors[1]	$strand2rgb{$strand1}";
	}elsif($chimeric_type=~/Deletion/){
		my @wide2coors=sort {$a <=> $b} ($refEnd1, $refStart2);
		$bed_row="$refName1	".($wide2coors[0]-1)."	$wide2coors[1]	$refName1:$strand1:$wide2coors[0]-$wide2coors[1]:$read_count	$read_count	$strand1	".($wide2coors[0]-1)."	$wide2coors[1]	$strand2rgb{$strand1}	2	1,1	0,".($wide2coors[1]-$wide2coors[0]);
	}elsif($chimeric_type=~/Reversion/){
		my @wide2coors=sort {$a <=> $b} ($refEnd1, $refStart2);
		my @all_4coor=($wide2coors[0]-$bed_thin_size, @wide2coors, $wide2coors[1]+$bed_thin_size); #all 1-based coordinates
		$all_4coor[0]=1 if $all_4coor[0]<1;
		$bed_row="$refName1	".($all_4coor[0]-1)."	$all_4coor[3]	$refName1:$strand1:$wide2coors[0]-$wide2coors[1]:$read_count	$read_count	$strand1	".($wide2coors[0]-1)."	$wide2coors[1]	$strand2rgb{$strand1}	2	".($all_4coor[1]-$all_4coor[0]+2).",".($bed_thin_size+1)."	0,".($all_4coor[3]-$all_4coor[0]);
	}elsif($chimeric_type=~/fusion/){
		my $bedName="$refName1:$strand1:${refEnd1}_$refName2:$strand2:$refStart2:$read_count";
		$bed_row=convert_fusion_coor_to_1bedRow($refName1,$strand1, $refEnd1, $bed_thin_size, $bedName, $read_count)."\n".convert_fusion_coor_to_1bedRow($refName2,$strand2, $refStart2, $bed_thin_size, $bedName, $read_count);
	}
	return($bed_row);
}

sub convert_fusion_coor_to_1bedRow{
	my ($refName,$strand, $fusion_pos, $bed_thin_size, $bedName, $read_count)=@_;
		if($strand eq "+"){
			my $bed_start=$fusion_pos-$bed_thin_size-1; #0-based 
			$bed_start=0 if $bed_start<0;
			$bed_row="$refName	$bed_start	$fusion_pos	$bedName	$read_count	$strand	".($fusion_pos-1)."	$fusion_pos	$strand2rgb{$strand}";
		}else{ # - strand
			$bed_end=$fusion_pos+$bed_thin_size;
			$bed_row="$refName	".($fusion_pos-1)."	$bed_end	$bedName	$read_count	$strand	".($fusion_pos-1)."	$fusion_pos	$strand2rgb{$strand1}";
		}
}