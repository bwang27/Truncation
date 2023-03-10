##1.1, bowtie index:
bowtie-build [options]* <reference_in> <ebwt_base>

##genome sequence and rDNA sequences (combined)
bowtie-build $root_dir/data/ucsc/genomes/mm8/mm_mm8.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_mm8/mm8 &
bowtie-build $root_dir/data/ucsc/genomes/hg18/hs_hg18.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_hg18/hg18 &

ln -s $root_dir/data/ucsc/genomes/mm8/mm_mm8.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_mm8/mm8.fa
ln -s $root_dir/data/ucsc/genomes/hg18/hs_hg18.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_hg18/hg18.fa






##1.2 bowtie2 index 12/29/2011
bowtie2-build [options]* <reference_in> <bt2_base>
##1.2.1 genome sequence (or with rDNA sequences combined)
bowtie2-build $root_dir/data/ucsc/genomes/mm8/mm_mm8.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_mm8/mm8 &
bowtie2-build $root_dir/data/ucsc/genomes/hg18/hs_hg18.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_hg18/hg18 &
ln -s $root_dir/data/ucsc/genomes/mm8/mm_mm8.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_mm8/mm8.fa
ln -s $root_dir/data/ucsc/genomes/hg18/hs_hg18.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_hg18/hg18.fa

#galGal4
geno=rheMac2 #sacCer3 galGal4 rheMac3 rheMac2
mkdir $root_dir/data/ucsc/genomes/$geno
mkdir -p $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_$geno

fa_file=$root_dir/data/ucsc/genomes/$geno/$geno.fa
index_f=$root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_$geno/$geno
#ln -s /Home/Data/UCSC/$geno/chromosomes/$geno.fa $fa_file
#ln -s /Home/Data/UCSC/$geno/chromosomes/$geno.fa $index_f.fa
ln -s $fa_file $index_f.fa
bowtie2-build $fa_file  $index_f &


##sacCer3 and mm9/hg19 combine
geno1=mm9 #mm9
cd $root_dir/analyze/dep_seq/1mapping/7map_index
mkdir bowtie2_${geno1}_sacCer3
cd bowtie2_${geno1}_sacCer3
grep ">" $root_dir/data/ucsc/genomes/sacCer3/sacCer3.fa
sed -e 's/^>/>y/' $root_dir/data/ucsc/genomes/sacCer3/sacCer3.fa  >sacCer3.2.fa
grep ">" sacCer3.2.fa
cat $root_dir/data/ucsc/genomes/${geno1}/${geno1}.fa sacCer3.2.fa  >${geno1}_sacCer3.fa
bowtie2-build  ${geno1}_sacCer3.fa  ${geno1}_sacCer3 &


#####PAL-seq tail length standard sequences and mm9 or hg19 combined index
geno1=hg19 #mm9 hg19
cd $root_dir/analyze/dep_seq/1mapping/7map_index
addFasta_f=/HPCTMP_NOBKUP/wl314/data/other/PAL-seq/Tail_length_standard_sequences.fa
mkdir bowtie2_${geno1}_PALseq
cd bowtie2_${geno1}_PALseq
cat $root_dir/data/ucsc/genomes/${geno1}/${geno1}.fa $addFasta_f  >${geno1}_PALseq.fa
bowtie2-build  ${geno1}_PALseq.fa  ${geno1}_PALseq &



##1.2.2 build refseq transcript index 06/06/2012
bowtie2-build $root_dir/analyze/club/23gene_subRegSeq/mm9/refFlat.exons.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_mm9/refseq &


##1.2.3  rDNA sequences alone
cd  $root_dir/data/ncbi/rRNA/mouse/
cat $root_dir/data/ncbi/rRNA/mouse/*.fa >mm_rDNA.fa 
bowtie2-build $root_dir/data/ncbi/rRNA/mouse/mm_rDNA.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_ribo/mm_rDNA &
bowtie2-inspect -s  $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_ribo/mm_rDNA

bowtie2-build $root_dir/data/ucsc/genomes/hg18/U13369.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_ribo/hs_rDNA &

#check index
bowtie2-inspect -s  $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_mm8/mm8

#yeast S. cerevisiae
species=S_cerevisiae
fastafile=$root_dir/data/ncbi/rRNA/$species/rDNA.less_redundant.fa
fastafile=$root_dir/data/ncbi/rRNA/S_cerevisiae/rDNA.representative.fa
fastafile=$root_dir/data/ncbi/rRNA/S_cerevisiae/JQ277730.fa
out_f=$root_dir/analyze/dep_seq/1mapping/7map_index/bowtie2_ribo/${species}_rDNA
bowtie2-build $fastafile $out_f 
bowtie2-inspect -s  $out_f





##2, build mouse LTR index from Genome Biology 2004 paper
bowtie-build $root_dir/data/other/LTR/GB2004_mm.LTR.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_mm/LTR &
ln -s $root_dir/data/other/LTR/GB2004_mm.LTR.fa $root_dir/analyze/dep_seq/1mapping/7map_index/bowtie_mm/LTR.fa






##2, novoalign index
cd $root_dir/analyze/dep_seq/1mapping
novoindex 7map_index/novo_mm8/mm8 $root_dir/data/ucsc/genomes/mm8/mm_mm8.fa
mkdir 7map_index/novo_mm9
novoindex 7map_index/novo_mm9/mm9 $root_dir/data/ucsc/genomes/mm9/mm_mm9.fa




##3, STAR index:
cd $pipeline_root/ReferenceDB/map_index/STAR
geno=hg19 #hg19 mm9
junc_f=$root_dir/analyze/club/0remove_club_redundance/8club_junc_info/1109/$geno.1109.STAR.junc
geno=canFam3 # rheMac2 rn6 canFam3
junc_f=$root_dir/analyze/club/0remove_club_redundance/8club_junc_info/$geno.ensGene.STAR.junc
geno=rheMac10 # rheMac10
junc_f=""

mkdir -p $geno/
STAR --runMode genomeGenerate --genomeDir $geno/ --genomeFastaFiles ../../ucsc/genomes/$geno/$geno.fa --runThreadN 8 \
	--sjdbFileChrStartEnd $junc_f 
	
mkdir -p $geno.noJunc/
STAR --runMode genomeGenerate --genomeDir $geno.noJunc/ --genomeFastaFiles ../../ucsc/genomes/$geno/$geno.fa --runThreadN 8


##create hg19 + plasmid plasmid_UGA20 index
cd $root_dir/analyze/dep_seq/1mapping
geno=hg19 #hg19 mm9
out_folder=7map_index/STAR/${geno}_UGA20/
junc_f=$root_dir/analyze/club/0remove_club_redundance/8club_junc_info/1109/$geno.1109.STAR.junc
addi_fa=/drive2/wli/analyze/dep_seq/1mapping/7map_index/STAR/hg19_UGA20/plasmid_UGA20.fa

STAR --runMode genomeGenerate --genomeDir $out_folder --genomeFastaFiles $root_dir/data/ucsc/genomes/$geno/$geno.fa $addi_fa --runThreadN 8 \
  --sjdbFileChrStartEnd $junc_f 
