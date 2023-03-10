##2022/06/24 WL: changed way of demultiplexing, use both qcat and the results from the sequencer (MinKNOW)
##2022/07/07 WL: added step 8b, do futher classification of chimeric reads

##settings of dataset (make changes based on the dataset)
study_name=2022-04-19_AAV_SEQ
seq_dir=/home/Research_Archive/RawData/Gene_Therapy/Nanopore/$study_name/
seq_sum_f=$seq_dir/sequencing_summary_FAS00538_c41b1205.txt
src_dir=/home/Research_Archive/RawData/Gene_Therapy/Nanopore/$study_name/fastq_pass
ont_kit_name="NBD104/NBD114" #used for qcat for demultiplexing step #see 'qcat -h': Auto,DUAL,NBD103/NBD104,NBD104/NBD114,NBD114,PBC001,PBC096,PBK004/LWB001,RAB204,RAB204/RAB214,RAB214,RBK001,RBK004,RPB004/RLB001,VMK001
nCPU=8

##main settings:
samples=(ID23_scAAV-SCA3-4E10x2_AAV SC-AAV9-SCA3-E10x2-CAGDel7-ID43-01-02)
helpers=(pAdDeltaF6 pHelper-KanR)
RepCaps=(pAAV-Rep2Cap9-KanR pAAV-Rep2Cap9-KanR)
GOIs=(ID23_scAAV-SCA3-4E10x2 scAAV-SCA3-E10x2_CAGdel7)
GOIfolders=(SCA3 SCA3)
colors=(orange green)
barcodes=(NBD01 NBD02)


##common settings (usually no need to change)

MinIONQC_program=/drive2/wli/software2/biosoft/seq_ana/minion_qc/MinIONQC.R

GOI_fa_dir=/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences
index_folder=/home/Research_Archive/ProcessedData/Nanopore_seq/target_index

root_dir=/drive2/wli/
pipeline_root=/drive2/wli/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
AAVseqPipeline_root=/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq
mapping_data_root=/home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/
project_root=$mapping_data_root #project folder (all the input and output files related to this project)
sample_info_f=$project_root/$study_name-Samples.xlsx
MinIONQC_out_dir=$mapping_data_root/MinIONQC/$study_name/


aav_sample_ids=(`seq 1 ${#samples[*]}`)

all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
done
echo ${all_AAV_samples[*]}


##1, QC
out_dir=$project_root/MinIONQC/
mkdir -p $out_dir
cd $out_dir

/bin/Rscript $MinIONQC_program  -i $seq_sum_f -o $out_dir -s TRUE


##2, demultiplexing:
mkdir -p $mapping_data_root/fastq/demultiplex
cd $mapping_data_root/fastq/
#cat $src_dir/*/*.fastq | qcat -b ./ --trim  -k NBD104/NBD114 >demultiplex.stats.txt 2>&1  & 
p=0
for i in `seq 1 ${#barcodes[*]}`; do
  indx=$(($i-1));
  barcodei=${barcodes[indx]}
  barcodei_num=`echo $barcodei | sed -e 's/NBD//' -e 's/BC//'`
  echo $barcodei_num
  mkdir -p $mapping_data_root/fastq/demultiplex/barcode$barcodei_num/
  cd $mapping_data_root/fastq/demultiplex/barcode$barcodei_num/
  cat $src_dir/barcode$barcodei_num/*.fastq | qcat -b ./ --trim  -k $ont_kit_name  >demultiplex.stats.txt 2>&1  & 
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait


##rename samples
cd $mapping_data_root/fastq/
mkdir renamed
cd renamed
for i in `seq 1 ${#barcodes[*]}`; do
  indx=$(($i-1));
  barcodei=${barcodes[indx]}
  sample1=${samples[indx]}
  barcodei_num=`echo $barcodei | sed -e 's/NBD//' -e 's/BC//'`
  mv $mapping_data_root/fastq/demultiplex/barcode$barcodei_num/barcode$barcodei_num.fastq $sample1.fastq
done


###use this "only" if qcat doesn't work
#  cd $mapping_data_root/fastq/
#  mkdir renamed
#  cd renamed
#  for i in `seq 1 ${#barcodes[*]}`; do
#    indx=$(($i-1));
#    barcodei=${barcodes[indx]}
#    sample1=${samples[indx]}
#    barcodei_num=`echo $barcodei | sed -e 's/NBD//' -e 's/BC//'`
#    cat $src_dir/barcode$barcodei_num/*fastq >$sample1.fastq
#  done



#run convert_snapgeneFeatures_2_bed_gtf_GOI_etc.R
#run Build_minimap2_index.R

# minimap2 folder: export PATH="$PATH:"`pwd`
# samtools folder: export PATH="$PATH:"`pwd`
# Use Build_minimap2_index.R(/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq/) to build an index folder


##3, calculate total read number, and gzip fastq files
cd /home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/fastq/renamed

out_readNum_f=readnum.txt
echo -e "Sample\tRaw read" >$out_readNum_f
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample1=${samples[indx]}
  rowNum=`wc -l $sample1.fastq | sed s/[[:blank:]].*//`
  rawReadNum=$(( $rowNum / 4 ))
  echo -e "$sample1\t$rawReadNum" >>$out_readNum_f
  gzip -f $sample1.fastq &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait


##4, mapping

##4.1 map to GOI-vector first
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/renamed/ #the directory of the fastq files 

for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.fastq.gz
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.map2GOI.sam
  minimap2_index_f=$GOI_fa_dir/${GOIfolders[indx]}/${GOIs[indx]}.fa
  mkdir $out_dir
  minimap2 -ax map-ont -t $nCPU $minimap2_index_f $fq1 >$out_sam_f
done


##4.2, find unmapped and mapped reads (to GOI plasmid) and output fastq format 
#2022/06/30 modified the code for output mapped reads, current code have redundant output for a same read
#2022/07/08 need to modify the code to output the original fastq, current code will output the reverse complimentary sequence of mapped reads if it map to the - strand (for minimap2 mapping results)
cd $pipeline_root/RNAseq

p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "${samples[indx]}"  -d "$mapping_data_root/minimap2/" -o "$project_root/01.ReadTable/map2GOI/" \
    -e ".map2GOI.sam" -t minimap2 -q 0 -m 30 -u 0 -U $mapping_data_root/fastq/${samples[indx]}.unmapped2GOI.fq -M $mapping_data_root/fastq/${samples[indx]}.mapped2GOI.fq &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait

cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable/map2GOI/


# ##5.1 map to all other GOI sequences (performed in same sequencing batch), 5.1 and 5.2 can be skipped if the GOI is the only GOI in the sequencing batch
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 
 
 
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.unmapped2GOI.fq
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.map2otherGOI.sam
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  minimap2_index_f=$project_root/othGOI_fa/$index_name.othGOI.fa
  mkdir $out_dir
  minimap2 -ax map-ont -t $nCPU $minimap2_index_f $fq1 >$out_sam_f
done

# ##5.2, find unmapped  reads (to other GOI plasmid) and output in fastq format
cd $pipeline_root/RNAseq
 
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "${samples[indx]}"  -d "$mapping_data_root/minimap2/" -o "$project_root/01.ReadTable/map2othGOI/" \
    -e ".map2otherGOI.sam" -t minimap2 -q 0 -m 30 -u 0 -U $mapping_data_root/fastq/${samples[indx]}.unmapped2GOI.fq  &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait
# 
cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable/map2othGOI/


##6, combine the reads mapped to GOI and not mapped to all other GOIs and map to the combined genome index with GOI+all other sequences (with sequences homologoud to GOI masked)
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 

for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.mapped2GOI.fq
  fq2=$fq_dir/${samples[indx]}.unmapped2GOI.fq
  fq3=$fq_dir/${samples[indx]}.filtered.fq #this is the final read file used for mapping
  cat $fq1 $fq2 >$fq3 #combine fastq files
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.sam
  echo out_sam_f
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  minimap2_index_f=$index_folder/$index_name/$index_name.mmi #make sure this index is GOI + masked all other references
  mkdir $out_dir
  minimap2 -ax map-ont -t $nCPU $minimap2_index_f $fq3 >$out_sam_f
done

##7, convert sam to bam
cd $mapping_data_root
out_dir=minimap2
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  samtools view -b -o $out_dir/${samples[indx]}.bam $out_dir/${samples[indx]}.sam &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait


##8, calculate read mapping number and total coverage length to different plasmids or chromosomes; also output read terminals, IR reads, chimeric reads information
cd $AAVseqPipeline_root
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  GOI_info_f=$index_folder/$index_name/$index_name.feature.txt

  perl 02.count_AAV_read_fromBam.pl -n $study_name -s "$sample"  -b "$mapping_data_root/minimap2/$sample.bam" \
    -p "$project_root/02_Bam2ReadCounts/"  \
    -t minimap2 -u 10 -m 30 -c 0.4 -v "$GOI_info_f" \
    -S "$mapping_data_root/minimap2/$sample.highQuality.sam" &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait

##combine mapping summary
cd $pipeline_root/RNAseq/
Rscript 02a.combinePE.readstats.R -in_dir "$project_root/02_Bam2ReadCounts/" -out_f $project_root/02_Bam2ReadCounts/All.read.counts.txt \
  -inf_pattern .read.counts.txt.log 

##8b, do futher classification of chimeric reads
cd $AAVseqPipeline_root
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  GOI_info_f=$index_folder/$index_name/$index_name.feature.txt

  perl 02b.ana_chimeric_reads.pl -n $study_name -s "$sample" -i "$project_root/02_Bam2ReadCounts/$sample.read.map.chimeric.txt" \
    -p "$project_root/02b_ChimericRead/"  \
    -S $ref_sizes_f -g "${GOIs[indx]}" &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait

##8c, further analyze the output of 02b.ana_chimeric_reads.pl, find the best threshold to define chimeric-reads (optional)
cd $AAVseqPipeline_root
Rscript 02c.ana_chimericReads.R -study_name $study_name -samples "${samples[*]}" \
    -read_in_dir $project_root/02_Bam2ReadCounts/  -out_dir $project_root/03_ana_readMap/$sample. \
    -ref_sizes_f $ref_sizes_f -GOI_info_f $GOI_info_f


##9, create bam index (and bigwig) files for visualization:
geno=hg19
PROG_DIR=$pipeline_root/Software/UCSC
fext='.highQuality.sam' #string after sample name in input sam file name

function sam2bw {
 #totalReadNum=`samtools view -c -q 10 $sample$fext`
 #echo "# $sample for file $sample$fext, TotalReadNum=$totalReadNum "
 echo "# $sample, convert to bam "
 samtools view -S -b -T $genofasta -o $sample.highQuality.bam  $sample$fext 
 #echo "# $sample filter original bam file (keep uniquely mapped reads only) "
 #samtools view -b -q 10 -o $sample.filtered.bam $sample.bam
 
 echo "# $sample, sort bam file"
 samtools sort -o $sample.sorted.bam $sample.highQuality.bam
 
 echo "# $sample, create index file for bam"
 samtools index $sample.sorted.bam
 #rm $sample.bam
 
 echo "# $sample, convert to bedgraph"
 genomeCoverageBed -bg -split -ibam $sample.sorted.bam -g $genofasta.sizes > $sample.bedgraph
 echo "# $sample, convert to bw"
 $PROG_DIR/bedGraphToBigWig $sample.bedgraph $genofasta.sizes $sample.bw 
}


p=0
cd $mapping_data_root/minimap2/
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  genofasta=$index_folder/$index_name/$index_name.fa
  sam2bw  &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPU ]; then  p=0 ;  wait; fi
done
wait

cd $mapping_data_root/minimap2/
# rm *.bedgraph
# rm *.sam
# rm *.map2GOI.bam
# rm *.map2oths.bam
# rm *.highQuality.bam


##10, analyze read mapping results
cd $AAVseqPipeline_root

for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  GOI_info_f=$index_folder/$index_name/$index_name.feature.txt

  Rscript 03.ana_readMapping.R -study_name $study_name -samples "$sample" \
    -read_in_dir $project_root/02_Bam2ReadCounts/  -out_dir $project_root/03_ana_readMap/$sample. \
    -ref_sizes_f $ref_sizes_f -GOI_info_f $GOI_info_f
done


# PATH=/home/wli/soft/pandoc/pandoc-2.2.1/bin:$PATH  
##11, draw read coverage in different chromosomes:
cd $AAVseqPipeline_root

Rscript 04GeneReadCovPlot.R -study_name $study_name -sample_info_f $sample_info_f \
   -mapping_data_root $mapping_data_root -index_folder $index_folder


##12, call peaks of human genome mapping

##or 12.b use a custom R script to cluster reads
cd $AAVseqPipeline_root
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  Rscript 07clusterReads.R -study_name $study_name -samples "$sample" \
    -read_in_dir $project_root/02_Bam2ReadCounts/  -out_dir $project_root/07_readCluster/$sample. \
    -ref_sizes_f $ref_sizes_f 
done

all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
done
echo ${all_AAV_samples[*]}

##cluster all samples:
  i=${aav_sample_ids[0]}
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  Rscript 07clusterReads.R -study_name $study_name -samples "${all_AAV_samples[*]}" \
    -read_in_dir $project_root/02_Bam2ReadCounts/  -out_dir $project_root/07_readCluster/AllSamples. \
    -ref_sizes_f $ref_sizes_f 

##12.c annotate peaks or clusters (bed format file)
cd $AAVseqPipeline_root

Rscript 08annotate_cluster.R -study_name $study_name -samples "AllSamples ${all_AAV_samples[*]}"

#add: PATH=/home/wli/soft/pandoc/pandoc-2.2.1/bin:$PATH ; to correct the graph size

##13.perform mutation analysis
#cd $AAVseqPipeline_root
cd $project_root
all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
done
echo ${all_AAV_samples[*]}

 i=${aav_sample_ids[0]}
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  genome=$index_folder/$index_name/$index_name.fa

out_dir=mutation
mkdir $out_dir
mkdir $sample

source /drive2/bwang/software/2022miniconda3/miniconda3/bin/activate nanocaller_env

/drive2/bwang/software/NanoCaller/NanoCaller --bam $mapping_data_root/minimap2/$sample.sorted.bam --ref $genome --cpu 3 --mode all --output $out_dir/$sample --prefix $sample

source /drive2/bwang/software/2022miniconda3/miniconda3/bin/deactivate nanocaller_env

##14, generate sequencing report

cd $AAVseqPipeline_root
/bin/Rscript 05prepare_report.R -study_name $study_name -samples "${samples[*]}" -barcodes "${barcodes[*]}" -project_root $project_root \
  -mapping_data_root $mapping_data_root -MinIONQC_out_dir "$MinIONQC_out_dir" -sample_info_f $sample_info_f \
  -AAV_samples "AllSamples ${all_AAV_samples[*]}"

cd $project_root/05.Report
for sample in All ${samples[*]}; do
  echo $sample
  /bin/Rscript -e "out_sample='$sample'; rmarkdown::render('NanoporeSeq_report.Rmd', output_format='html_document', output_file = 'NanoporeSeq_report.$sample.html') " 
  /bin/Rscript -e "out_sample='$sample'; rmarkdown::render('NanoporeSeq_report.Rmd', output_format='word_document', output_file = 'NanoporeSeq_report.$sample.docx') " 
done
