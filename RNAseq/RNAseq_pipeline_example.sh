
##this pipeline start with fastq files, do mapping (using STAR), gene expression quantification and differential gene expression analysis (using DESeq2)
#programs needed to install in linux:
# perl, R
# STAR
# samtools
# bedtools
# R library: DESeq2, ggplot2, edgeR, BiocParallel, fmsb, gplots, 


root_dir=/HPCTMP_NOBKUP/wl314/
study_name=Feng_p14p16_RD_RNAseq #project name
pipeline_root=$root_dir/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
project_root=$root_dir/analyze/Projects/TJ/$study_name/ #project folder (all the input and output files related to this project)
mapping_data_root=$root_dir/analyze/Projects/TJ/$study_name/ #use this directory to save large files (eg. fastq files and mapped bam files). expect a fastq folder with all fastq files under this directory

source $pipeline_root/setGlobalVars.sh
if_PE=1; #whether this is pair end sequencing ()
if_stranded=1; #whether the RNA library is stranded (contain gene strand information)
if_read1_antiSense=1; #whether the first read is antisense to gene. This is usually true for TruSeq Stranded mRNA Library (second read is sense strand)
if_do_splicing_analysis=1; #thether do splicing analysis

####R.1, mapping
##using STAR to map
study_name=Feng_p14p16_RD_RNAseq
geno=hg19
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"

cd $mapping_data_root
mkdir STAR_map
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 
genomeDir=$pipeline_root/ReferenceDB/map_index/STAR/$geno #the directory of mapping index file


if [ $if_PE -eq 1 ]; then
   echo "#treat pair-end read..."
   ## unzip (in some server, STAR will report error if not using unzipped files)
   # p=0
   # for sample1 in $samples; do
   #   fq1=$fq_dir/${sample1}.1.fq.gz
   #   fq2=$fq_dir/${sample1}.2.fq.gz
   #   gunzip $fq1 &
   #   gunzip $fq2 &
   #   p=$(($p+2)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
   # done
   # wait
   
   samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
   for sample1 in $samples; do
     fq1=$fq_dir/${sample1}.1.fq.gz
     fq2=$fq_dir/${sample1}.2.fq.gz
     out_dir=STAR_map/$sample1/
     mkdir $out_dir
     STAR --runMode alignReads -f --genomeDir $genomeDir --readFilesIn $fq1 $fq2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate \
       --genomeLoad LoadAndKeep --limitBAMsortRAM  10000000000 \
       --readFilesCommand "gunzip -c" \
       --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD jM jI XS --outFileNamePrefix $out_dir
   done

else #single end reads
   echo "#treat single-end read..."
   # p=0
   # for sample1 in $samples; do
   #   fq1=$fq_dir/${sample1}.fq.gz
   #   gunzip $fq1 &
   #   p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
   # done
   # wait
   
   samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
   for sample1 in $samples; do
     fq1=$fq_dir/${sample1}.fq.gz
     out_dir=STAR_map/$sample1/
     mkdir $out_dir
     STAR --runMode alignReads -f --genomeDir $genomeDir --readFilesIn $fq1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate \
       --genomeLoad LoadAndKeep --limitBAMsortRAM  10000000000 \
       --readFilesCommand "gunzip -c" \
       --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD jM jI XS --outFileNamePrefix $out_dir
   done
fi


STAR  --genomeDir $genomeDir --genomeLoad Remove  #release memory
rm -fr _STARtmp
rm Log.out
rm Log.progress.out
rm Aligned.out.sam


cd STAR_map
out_readnum_f=Input.read.num.log.out
echo -e "Sample\tRaw read" >$out_readnum_f
grep "Number of input reads" */Log.final.out  >>$out_readnum_f
sed -i -e 's/\/Log.final.out:\s*Number of input reads |//g' $out_readnum_f


##create index for IGV visualization
cd $mapping_data_root/STAR_map
p=0
for sample1 in $samples; do
  ln -f $sample1/Aligned.sortedByCoord.out.bam $sample1.bam
  echo $sample1
  samtools index $sample1.bam &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done
wait





####R.2, analyze RNA-seq data
####R2.1  find uniquely mapped reads with good quality from sam or bam file and output table format files for downstream analysis
study_name=Feng_p14p16_RD_RNAseq
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
cd $pipeline_root/RNAseq

p=0
for sample in $samples; do
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "$sample"  -d "$mapping_data_root/STAR_map/" -o "$project_root/01.ReadTable/" \
    -e ".bam" -t STAR -q 10 -m 5 -u 1  &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done
wait


cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable



if [ $if_PE -eq 1 ]; then
  ####R2.2  extract pair-end reads mapping info
  cd $pipeline_root/RNAseq
  samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
  study_name=Feng_p14p16_RD_RNAseq
  
  p=0
  for sample in $samples; do
    perl 02.Cal_PEreadsCov_frSam.pl -s $study_name -S "$sample"  -i "$mapping_data_root/STAR_map/" -e ".bam" \
      -q 10  -d 0 -t STAR -o "$project_root/02.PE.ReadTable/" &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
  done
  wait
  
  #combine PE read stats
  cd $pipeline_root/RNAseq
  Rscript 02a.combinePE.readstats.R -in_dir "$project_root/02.PE.ReadTable/"
  
  
  ####R2.3a calculate gene expression using 14cal_gene_cds_readsnum.pl
  study_name=Feng_p14p16_RD_RNAseq
  geno=hg19
  cd $pipeline_root/RNAseq
  samples=(siA1 siA2 siB1 siB2 siCtl1 siCtl2)
  p=0
  for i in `seq 1 ${#samples[*]}`; do
    indx=$(($i-1))
    perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i refseqid -p refseqcds_PE. \
      -e CDS -d "$if_stranded" -r "$if_read1_antiSense" -a "${samples[indx]}" \
      -D $project_root/02.PE.ReadTable/ -E reads.cov.tbl -o $project_root/04.GeneReadNum/ &
    p=$(($p+1));  if [ "$p" -ge 6 ]; then p=0; wait;  fi
  done
  wait
  
  cd $pipeline_root/RNAseq/
  perl 04a.cal_readMapStats.pl -o $project_root/04.GeneReadNum/ -p refseqcds_PE
fi



####R2.3b treat all reads as single end (SE) read, calculate junction read map to gene info (as well as gene expression based on CDS region of refseq transcripts)
##the output of this step is not used for PE reads gene expression analysis. bjut it will be used for splcing analysis for both PE and SE reads
geno=hg19
study_name=Feng_p14p16_RD_RNAseq
cd $pipeline_root/RNAseq
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
p=0
for sample in $samples; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i refseqid -p refseqcds_SE. -e CDS -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 1 -E "JunReadNum GenoMapNum" -o $project_root/04.GeneReadNum/ &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $pipeline_root/RNAseq/
perl 04a.cal_readMapStats.pl -o $project_root/04.GeneReadNum/ -p refseqcds_SE



####R2.4 calculate gene expression regulation 
study_name=Feng_p14p16_RD_RNAseq
geno=hg19
cd $pipeline_root/RNAseq
if [ $if_PE -eq 1 ]; then
  gene_rnum_f=$project_root/04.GeneReadNum/refseqcds_PE.ReadNum.tbl
else
  gene_rnum_f=$project_root/04.GeneReadNum/refseqcds_SE.ReadNum.tbl
fi

sample_repl_l="siA:siA1 siA2;siB:siB1 siB2;siCtl:siCtl1 siCtl2;"
test_samples="siA siB"
ref_samples="siCtl siCtl"
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
out_prefix=""
rnum_idname=refseqid

Rscript 05.cal_gexch.template.R -p_cut 0.05 -pval_fun "DESeq2" -foldchange_cut 1.5 -IfCal_RPM 0 \
  -gene_model refseqcds  -geno $geno -study_name $study_name -gene_rnum_f $gene_rnum_f \
  -output_root  $project_root/05.Gene_DE/ \
  -nCores 8 \
  -all_sample_name "$samples" \
  -gene_len_var cds_len -test_samples "$test_samples" -ref_samples "$ref_samples" \
  -g_ano_f "../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt" -rnum_idname "$rnum_idname" \
  -comb_sample_l "$sample_repl_l" \
  -out_prefix "$out_prefix" \
  -what2do  all

  # -what2do  all
  # -what2do  formatGeneTb
  # -what2do  clust_rawSample_exp




####R2.5, create bigwig files for visualization:
study_name=Feng_p14p16_RD_RNAseq
geno=hg19
PROG_DIR=$pipeline_root/Software/UCSC
genofasta=$pipeline_root/ReferenceDB/ucsc/genomes/$geno/$geno.fa
fext='.bam' #string after sample name in input sam file name

function sam2bw {
 totalReadNum=`samtools view -bc -q 10 $sample$fext`
 #totalReadNum=`wc -l $sample$fext | sed s/[[:blank:]].*//`
 echo "step2, ____ $sample for file $sample$fext, TotalReadNum=$totalReadNum ____"
 #echo step3, ____ $sample, convert to bam ____
 #samtools view -S -b -T $genofasta -o $sample$fext.bam  $sample$fext 
 echo "step3, ____ $sample filter original bam file (keep uniquely mapped reads only) ____"
 samtools view -b -q 10 -o $sample.filtered.bam $sample$fext

 echo step4,  ____ $sample, sort bam file ____
 #samtools sort $sample$fext.bam  $sample$fext.bam.sorted ##!!!no .bam in output name
 echo step5, ____ $sample, convert to bedgraph ____
 genomeCoverageBed -bg -split -ibam $sample.filtered.bam -g $PROG_DIR/$geno.chrom.sizes > $sample.bedgraph
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand + -g $PROG_DIR/$geno.chrom.sizes > $sample.p.bedgraph #-split 
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand '-' -g $PROG_DIR/$geno.chrom.sizes > $sample.m.bedgraph #-split 
 echo step6, ____ $sample, normolize bedgraph counts ____
 norm_bedgraph.pl -t $totalReadNum -i "$sample.bedgraph" 
 echo step8, ____ $sample, convert to bw ____
 $PROG_DIR/bedGraphToBigWig $sample.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.p.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.p.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.m.bedgraph.normolized $PROG_DIR/$geno.chrom.sizes $sample.m.bw 
}


samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
p=0
cd $mapping_data_root/STAR_map/
for sample in  $samples; do 
  sam2bw  &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $mapping_data_root/STAR_map/
rm *.bedgraph
rm *.bedgraph.normolized
rm *.bam.unique
rm *.filtered.bam

<<create_UCSC_track_strings
  ####create track strings
  study_name=Feng_p14p16_RD_RNAseq
  colors=(000,000,125 000,000,255 125,000,150 125,000,020 255,000,000)
  samples=(siCtl1 siCtl2 siA1 siA2 siB1 siB2) 
  col_indexs=(0 0 2 2 3 3); priority_start=41
  
   for samp_i in `seq 1 ${#samples[*]}`; do
    indx=$(($samp_i-1))
    priority=$(( $priority_start+$indx ))
    echo track type=bigWig visibility=2 alwaysZero=on color=${colors[ ${col_indexs[$indx]} ]}   graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=$study_name \
     priority=$priority name=\"${samples[$indx]}\" description=\"RNA-seq ${samples[$indx]}\" \
     bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/$study_name/${samples[$indx]}.bw
   done
  
  track type=bigWig visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=41 name="siCtl1" description="RNA-seq siCtl1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siCtl1.bw
  track type=bigWig visibility=2 alwaysZero=on color=000,000,125 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=42 name="siCtl2" description="RNA-seq siCtl2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siCtl2.bw
  track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=43 name="siA1" description="RNA-seq siA1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siA1.bw
  track type=bigWig visibility=2 alwaysZero=on color=125,000,150 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=44 name="siA2" description="RNA-seq siA2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siA2.bw
  track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=45 name="siB1" description="RNA-seq siB1" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siB1.bw
  track type=bigWig visibility=2 alwaysZero=on color=125,000,020 graphType=bar maxHeightPixels=20:30:50 itemRgb=On group=Feng_p14p16_RD_RNAseq priority=46 name="siB2" description="RNA-seq siB2" bigDataUrl=http://www.bioinformatics.org/blogo/bigwig/Feng_p14p16_RD_RNAseq/siB2.bw
  #autoScale=off yLineOnOff=on viewLimits=0:6
create_UCSC_track_strings




if [ $if_do_splicing_analysis != 1 ]; then 
  echo "No splicing analysis was requested. Exit!"
  exit 0; 
fi
echo "Start doing splicing analysis..."


#####Splicing analysis pipeline
####S.1 combine junction read mapping info 
geno=hg19
study_name=Feng_p14p16_RD_RNAseq
cd $pipeline_root/RNAseq_Splicing
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
Rscript 01.comb_junc_map_info.R -study_name $study_name -geno $geno -samples "$samples" -indir $project_root/04.GeneReadNum/refseqcds_SE. \
  -out_cbf $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -id_header refseqid



####S.2, cut genes to blocks (exon and intron part) and do expression analysis
geno=hg19
cd $pipeline_root/RNAseq_Splicing
perl 02.cut_gene3blocks_basedOnSplicing.pl  -g $geno -R 4 -S 2 -F 0.05 -i refseqid \
  -j $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -o $project_root/RNAseq_Splicing/02.geneBlocks

##create bed and bigbed file
perl 02b.Gblock2bigBed.pl -g $geno -a $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl.ano


####S.3, extract splice site flanking region sequence
cd $pipeline_root/SequenceAnalysis
geno=hg19

#extract 100 bp upstream and downstream of splice site (this input file has more junctions, some of them are from gene model only, not from reads)
perl 01.extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s $project_root/RNAseq_Splicing/02.geneBlocks/AllJunc.tbl \
  -i "contig strand juncpos5" -e juncpos5 -n ss5_PM100 -c m99p100 \
  -o $project_root/RNAseq_Splicing/SS_flank_seq/ss5_PM100.fa &
 
perl 01.extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s $project_root/RNAseq_Splicing/02.geneBlocks/AllJunc.tbl \
  -i "contig strand juncpos3" -e juncpos3 -n ss3_PM100 -c m100p99 \
  -o $project_root/RNAseq_Splicing/SS_flank_seq/ss3_PM100.fa &
wait


####S.4, calculate gene block expression (exonic and intronic parts) 
geno=hg19
study_name=Feng_p14p16_RD_RNAseq
cd $pipeline_root/RNAseq
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
p=0
for sample in $samples; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i Gblock_id -p Gblock. -e transcript -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 0 -E "JunReadNum GenoMapNum" -o $project_root/RNAseq_Splicing/geneBlockReadNum/ \
    -A $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl  &
    p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait
cd $pipeline_root/RNAseq/
perl 04a.cal_readMapStats.pl -o $project_root/RNAseq_Splicing/geneBlockReadNum/ -p Gblock



####S.5 do splicing analysis using DEXSeq 
cd $pipeline_root/RNAseq_Splicing

study_name=Feng_p14p16_RD_RNAseq
geno=hg19


sample_repl_l="siA:siA1 siA2;siB:siB1 siB2;siCtl:siCtl1 siCtl2;"
test_samples="siA siB"
ref_samples="siCtl siCtl"
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"
ana_unit=geneBlock 
#note: need relatively large memory, do not use too many workers if memory is small
Rscript 05.ana_splicing_DEXseq.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -workers 2 \
   -gene_id_header refseqid \
   -in_cbf $project_root/RNAseq_Splicing/geneBlockReadNum/Gblock.ReadNum.tbl \
   -juncAno_f  $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl.ano \
   -as_ano_root $project_root/RNAseq_Splicing/02.geneBlocks/ASano. \
   -geneAno_f "../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt" \
   -ana_unit $ana_unit \
   -out_root $project_root/RNAseq_Splicing/05.DEXSeq \
   -what2do all -use_p_type Padj -P_cut 0.001 -Log2ratio_cut 1
 


####S.6, study all exon splicing using PSI
cd $pipeline_root/RNAseq_Splicing
study_name=Feng_p14p16_RD_RNAseq
geno=hg19


sample_repl_l="siA:siA1 siA2;siB:siB1 siB2;siCtl:siCtl1 siCtl2;"
test_samples="siA siB"
ref_samples="siCtl siCtl"
samples="siA1 siA2 siB1 siB2 siCtl1 siCtl2"

Rscript 06.ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -basal_PSI_samples "DMSO" \
   -gblock_f "$project_root/RNAseq_Splicing/05.DEXSeq/geneBlock/combine_d.tbl" \
   -read_num_f $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
   -out_root $project_root/RNAseq_Splicing/06.exon_PSI/ \
   -workers 8 \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 5 -delta_PSI_cut 10 \
   -what2do all


##study all annotated exons
Rscript 06.ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "$test_samples" -ref_samples "$ref_samples" -all_samples "$samples" \
   -basal_PSI_samples "DMSO" \
   -read_num_f $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
   -out_root $project_root/RNAseq_Splicing/06.exon_PSI/ \
   -workers 6 \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 5 -delta_PSI_cut 10 \
   -study_exon_type  annotated_exon -what2do all


