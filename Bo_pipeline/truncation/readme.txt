##step1
/drive2/bwang/software/bedtools2/bin/bedtools bamtobed -i scAAV-SCA3-4.bam -tag NM >scAAV-SCA3-4.bed
 awk '$1~/ssAAV-SCA3-3_shifted/ {print }' ssAAV-SCA3-3.bed >ssAAV-SCA3-3.step1.same.genome
 awk '{print $4}' scAAV-SCA3-4.step1.same.genome | sort | uniq -c | awk '$1>1 {print $2}' >scAAV-SCA3-4.step2.candidate.reads.ids
perl get.inverted.reads.pl ssAAV-SCA3-3.step1.same.genome ssAAV-SCA3-3.step2.candidate.reads.ids
perl step4_get_plus_minus_mapping_reads.pl  ssAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments >plus.minus.step4.ssAAV-SCA3-3.step1.same.genome.step3.sameread.all.alignments

perl modify.judge.the.read.inverted.or.not.pl plus.minus.step4.ssAAV-SCA3-4.step1.same.genome.step3.sameread.all.alignments 10

##step2
bedtools genomecov -ibam ID8_scAAV-SCA3-4_plasmid.sorted.bam -bga >ID8_scAAV-SCA3-4_plasmid.bedgraph
sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" change.bedgraph.to.continous.position.pl
perl change.bedgraph.to.continous.position.pl ID8_scAAV-SCA3-4_plasmid.bedgraph >tmp
mv tmp ID8_scAAV-SCA3-4_plasmid.bedgraph.continuous
bedtools genomecov -ibam ID8_scAAV-SCA3-4_plasmid.sorted.bam -bg -strand + >ID8_scAAV-SCA3-4_plasmid.bedgraph.plus
perl change.bedgraph.to.continous.position.pl ID8_scAAV-SCA3-4_plasmid.bedgraph.plus >tmp
mv tmp ID8_scAAV-SCA3-4_plasmid.bedgraph.plus.continuous
bedtools genomecov -ibam ID8_scAAV-SCA3-4_plasmid.sorted.bam -bg -strand - >ID8_scAAV-SCA3-4_plasmid.bedgraph.minus
perl change.bedgraph.to.continous.position.pl ID8_scAAV-SCA3-4_plasmid.bedgraph.minus >tmp
mv tmp ID8_scAAV-SCA3-4_plasmid.bedgraph.minus.continuous
awk '{print $1"\t"$2"\t"$3"\t"$4"\t""+"}' ID8_scAAV-SCA3-4_plasmid.bedgraph.plus.continuous >tmp
mv tmp ID8_scAAV-SCA3-4_plasmid.bedgraph.plus
 awk '{print $1"\t"$2"\t"$3"\t"$4"\t""-"}' ID8_scAAV-SCA3-4_plasmid.bedgraph.minus.continuous >tmp
mv tmp ID8_scAAV-SCA3-4_plasmid.bedgraph.minus
 cat ID8_scAAV-SCA3-4_plasmid.bedgraph.plus ID8_scAAV-SCA3-4_plasmid.bedgraph.minus >ID8_scAAV-SCA3-4_plasmid.bedgraph.strand

##step3

awk '$4~/+/ {print}' uniq.plus.minus.step4.ID8_scAAV-SCA3-4_plasmid.step1.same.genome.step3.sameread.all.alignments.sense.truncation.point >uniq.plus.minus.step4.ID8_scAAV-SCA3-4_plasmid.step1.same.genome.step3.sameread.all.alignments.sense.truncation.point_step2

awk '$4~/+/ {print}' uniq.plus.minus.step4.ID8_scAAV-SCA3-4_plasmid.step1.same.genome.step3.sameread.all.alignments.antisense.truncation.point >uniq.plus.minus.step4.ID8_scAAV-SCA3-4_plasmid.step1.same.genome.step3.sameread.all.alignments.antisense.truncation.point_step2

## check reads position and choose continuous reads to support IR reads
samtools view -h -o ID23_scAAV-SCA3-4E10x2_AAV.highQuality.sam ID23_scAAV-SCA3-4E10x2_AAV.highQuality.bam
paftools.js sam2paf ID23_scAAV-SCA3-4E10x2_AAV.highQuality.sam >ID23.paf

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' ID23_scAAV-SCA3-4E10x2_AAV.highQuality.paf >ID23_scAAV-SCA3-4E10x2_AAV.reads.coordinates

 perl get.each.truncation.reads.from.paf.pl ID23_scAAV-SCA3-4E10x2_AAV.reads.coordinates uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.sense.truncation_point_step2 uniq.plus.minus.step4.ID23_scAAV-SCA3-4E10x2_AAV.step1.same.genome.step3.sameread.all.alignments.antisense.truncation_point_step2

  perl judge.reads.position.continous.or.not.sense.pl IR.sense.reads.position 10
  perl judge.reads.position.continous.or.not.antisense.pl IR.antisense.reads.position 10

##step4

 sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" 05_tmp_get_truncation.point.sense.max.pl
 sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" 06_tmp_get_truncation.point.antisense.max.pl
 perl 05_tmp_get_truncation.point.sense.max.pl reads/IR.sense.reads.position.10.both.confirmed.continuousfor.downstream.analysis ID8_scAAV-SCA3-4_plasmid.bedgraph.continuous ID8_scAAV-SCA3-4_plasmid.truncation.sense.point.bed
 perl 06_tmp_get_truncation.point.antisense.max.pl reads/IR.antisense.reads.position.10.both.confirmed.continuousfor.downstream.analysis ID8_scAAV-SCA3-4_plasmid.bedgraph.continuous ID8_scAAV-SCA3-4_plasmid.truncation.antisense.point.bed

##step5

sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" 07_count_truncation_read_each_point.pl

perl code/revised_07_count_truncation_read_each_point_sense.pl reads/IR.sense.reads.position.10.both.confirmed.continuousfor.downstream.analysis | awk '{print $2}' | sort | uniq -c >number_of_truncation_reads_sense_per_point 
 
perl code/revised_07_count_truncation_read_each_point_antisense.pl reads/IR.antisense.reads.position.10.both.confirmed.continuousfor.downstream.analysis | awk '{print $2}' | sort | uniq -c >number_of_truncation_reads_antisense_per_point

 perl merge.truncation.reads.and.total.reads.for.all.pl number_of_truncation_reads_sense_per_point ID8_scAAV-SCA3-4_plasmid.truncation.sense.point.bed ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.sense

 perl merge.truncation.reads.and.total.reads.for.all.pl number_of_truncation_reads_antisense_per_point ID8_scAAV-SCA3-4_plasmid.truncation.antisense.point.bed ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.antisense

perl /drive2/bwang/analysis/truncation/revised_09_merge.plus.minus.table.pl ID23.truncation.vs.total.reads.sense ID23.truncation.vs.total.reads.antisense

 cat ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plus ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.minus >ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus
cat ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus | sort | uniq >uniq.ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus

sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" combine.truncation.with.normal.pl
 perl combine.truncation.with.normal.pl uniq.ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus ID8_scAAV-SCA3-4_plasmid.bedgraph.continuous

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""truncation"}' uniq.ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus >uniq.ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus.truncationlabel
cat uniq.ID8_scAAV-SCA3-4_plasmid.truncation.vs.total.reads.plusandminus.truncationlabel newnormal.reads.without.truncations.coverage >merged.truncation.normal.reads.coverage
 sort -k2,2n merged.truncation.normal.reads.coverage >sorted.merged.truncation.normal.reads.coverage

sed -i "s/ID23_scAAV-SCA3-4E10x2/ID8_scAAV-SCA3-4/g" get.ref.start.end.counts.by.strand.pl
nohup perl get.ref.start.end.counts.by.strand.pl ID8_scAAV-SCA3-4_plasmid.bed sorted.merged.truncation.normal.reads.coverage &

perl add.truncation.strand.plus.minus.pl plus.sorted.merged.truncation.normal.reads.coverage.withref.startend >final.summary.txt

perl IR.reads.plus.minus.2colum.pl final.summary.txt >ID23-Fraction8.final.summary.txt


