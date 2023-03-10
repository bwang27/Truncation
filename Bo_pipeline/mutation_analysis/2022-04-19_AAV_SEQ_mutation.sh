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


##13.perform mutation analysis
#cd $AAVseqPipeline_root
cd $project_root
out_dir=09_mutation #mutation_analysis
#mkdir $out_dir
cd $out_dir
all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  genome=$index_folder/$index_name/$index_name.fa
 # mkdir $sample
  source /drive2/bwang/software/2022miniconda3/miniconda3/bin/activate nanocaller_env
  /drive2/bwang/software/NanoCaller/NanoCaller --bam $mapping_data_root/minimap2/$sample.sorted.bam --ref $genome --cpu 3 --mode all --output $out_dir/$sample --prefix $sample
  cd $project_root/$out_dir/$sample
  bcftools view -e 'QUAL<=10 || DP<114 || FQ <0.5' $sample.snps.vcf.gz >$sample.filtered.snps.vcf.gz
done

source /drive2/bwang/software/2022miniconda3/miniconda3/bin/deactivate nanocaller_env

#echo ${all_AAV_samples[*]}
