##gather data and report the following: 
#1, sequencing QC; 
#2, demultiplexing results; 
#3, read coverage plots; 
#4, read mapping percentage and length data in different plasmids/chromosomes


args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../../../Pipeline/SharedCodes/Rfunc.inc.R")
library("yaml")
library("openxlsx")

if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }
if(!is.na(args_v["samples"])){ samples= unlist(strsplit(as.character(args_v["samples"])," "))}
if(!is.na(args_v["barcodes"])){ barcodes= unlist(strsplit(as.character(args_v["barcodes"])," "))}
AAV_samples=c("All",samples)
if(!is.na(args_v["AAV_samples"])){ AAV_samples= unlist(strsplit(as.character(args_v["AAV_samples"])," "))}
project_root=ifelse(is.na(args_v["project_root"]), paste0(study_name,"/"), args_v["project_root"]) 
peakAnno_out_prefix=paste0(project_root, "/08_ClusterAnnotation/hostCellDNA_anno")
if(!is.na(args_v["peakAnno_out_prefix"])){ peakAnno_out_prefix= as.character(args_v["peakAnno_out_prefix"]) }

mapping_data_root<- ifelse(is.na(args_v["mapping_data_root"]), project_root, args_v["mapping_data_root"]) 
demultiplexing_stat_f<- ifelse(is.na(args_v["demultiplexing_stat_f"]), paste0(mapping_data_root,"/fastq/renamed/readnum.txt"), args_v["demultiplexing_stat_f"]) 
MinIONQC_out_dir<- ifelse(is.na(args_v["MinIONQC_out_dir"]), paste0(mapping_data_root,"/MinIONQC/", study_name,"/"), args_v["MinIONQC_out_dir"]) 
sample_info_f<- ifelse(is.na(args_v["sample_info_f"]), paste0("/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/",study_name,"/", study_name,"-Samples.xlsx"), args_v["sample_info_f"]) 

min_used_readNum=ifelse(is.na(args_v["min_used_readNum"]), 500, as.numeric(args_v["min_used_readNum"]) )

report_out_dir=ifelse(is.na(args_v["report_out_dir"]), paste(project_root,"/05.Report/",sep=""), args_v["report_out_dir"]) 
out_log_f=paste(report_out_dir,"report.log",sep="")
out_R_inc_f=paste(report_out_dir,"report.inc.R",sep="")
out_data_img_f=paste(report_out_dir,"report.R.img",sep="")

chimeric_types=c("Inversion","Deletion","Insertion","Reversion","fusion")

mkdir_if_not_exist(out_data_img_f)
system( paste0("cp -f NanoporeSeq_report.Rmd",  " ", report_out_dir,"NanoporeSeq_report.Rmd") )

save.image(out_data_img_f)
total_read_num=NA

##0, sample information table:
sample_info_d=read.xlsx(sample_info_f)
sample_info_d=re_format_tb(sample_info_d, back_headers=c("Sample.description","Note"), ch_header_from="renamed.sample.name", ch_header_to="sample" )
report_impurity_references=c("GOI","GOI","pHelper","pRepCap","Human.Chromosomes","Human_Ad5_first5k","","")
report_impurity_regions=c("All","GOI_full","All","All","All","All","NeoR/KanR","AmpR")
report_impurity_headers=change_values(paste0(report_impurity_references,".",report_impurity_regions,"%"),c("GOI.GOI_full","GOI.All", ".All","^\\."),c("ITR-GOI-ITR","GOIplasmid","","")); report_impurity_headers
report_genes=c("NeoR/KanR","AmpR","KanR")

##1, MinIONQC results:
list.files(MinIONQC_out_dir)
MinIONQC_sum_f=paste0(MinIONQC_out_dir,"summary.yaml")

if(file.exists(MinIONQC_sum_f)){
	MinIONQC_sum_l <- yaml.load_file(MinIONQC_sum_f)
	MinIONQC_sum_d= data.frame(var=names(MinIONQC_sum_l$"All reads")[1:8], All.Reads=as.character(unlist(MinIONQC_sum_l$"All reads")[1:8]), 
		"Q>=7"=as.character(unlist(MinIONQC_sum_l$"Q>=7")[1:8]), check.names=F )
	MinIOQC_image_fs=list.files(MinIONQC_out_dir,".png")
	total_read_num=as.numeric(as.character(MinIONQC_sum_d[MinIONQC_sum_d$var=="total.reads","Q>=7"])); total_read_num

	MinIOQC_plot_names=c("length_histogram","q_histogram","length_vs_q","reads_per_hour","yield_over_time","length_by_hour","q_by_hour","yield_by_length","channel_summary","gb_per_channel_overview") # "flowcell_overview" not used
	#copy plot files to report folder:
	copy_files=intersect(paste0(MinIOQC_plot_names,".png"), MinIOQC_image_fs); copy_files
	mkdir_if_not_exist(paste0(report_out_dir,"MinIONQC/"))
	for(f1 in copy_files){
		system( paste0("cp -f ", MinIONQC_out_dir,f1, " ", report_out_dir,"MinIONQC/",f1) )
	}
}

##2, demultiplexing results:
if(file.exists(demultiplexing_stat_f)){
	demultiplexing_stat_d=read.table(demultiplexing_stat_f, header=T, sep="\t", stringsAsFactors=F, quote="", check.names=F)
	demultiplexing_stat_d$Barcode=barcodes[match(demultiplexing_stat_d[,1], samples)]
	total_d=demultiplexing_stat_d[1,]
	total_d[1,c(1,3)]=c("All", "All")
	if(is.na(total_read_num)){total_read_num=sum(demultiplexing_stat_d[,2])}
	total_d[1,2]=total_read_num
	demultiplexing_stat_d=rbind(demultiplexing_stat_d, total_d )
	demultiplexing_stat_d[,"Reads (%)"]=round( demultiplexing_stat_d[,2] / total_read_num *100, 2 )
}else{
	print(paste("Error: demultiplexing_stat_f",demultiplexing_stat_f, "does not exist!"))
	quit("no")
}


##3, Mapping overall results:
mapping_stats_f=paste0(project_root,"/02_Bam2ReadCounts/All.read.counts.txt")
mapping_stats_d=read.table(mapping_stats_f, header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
discard_cols=NULL
col_white_list=grep("Used",names(mapping_stats_d),value=T); col_white_list
for(i in 3:ncol(mapping_stats_d) ){
	if(all(mapping_stats_d[,i]==mapping_stats_d[,i-1]) & ! (names(mapping_stats_d)[i] %in% col_white_list) ){
		discard_cols=c(discard_cols,i)
	}
}
discard_cols
if(length(discard_cols)>0){
	mapping_stats_d=mapping_stats_d[,-discard_cols]
}
mapping_stats_d=re_format_tb(mapping_stats_d, ch_header_from=c(".reads",".clip<=none.mate1","primaryAli.",".exon","\\d\\."), 
	ch_header_to=c("","","","",""), del_pattern="collapsed")
mapping_stats_d=mapping_stats_d[match(mapping_stats_d[,1],samples), ]
mapping_stats_d$RawRead=demultiplexing_stat_d[match(mapping_stats_d[,1], demultiplexing_stat_d[,1]),2]
Used_read_headers=intersect(c("Used","Used.mate1","Used.mate2"),names(mapping_stats_d)); Used_read_headers

mapping_stats_d[,paste0(Used_read_headers,"/Read")]=round(unlist(mapping_stats_d[,Used_read_headers]) / mapping_stats_d$RawRead, 3)
mapping_stats_d=re_format_tb(mapping_stats_d,front_headers=c("Sample","RawRead"))
report_samples=intersect(samples, mapping_stats_d[rowSums(mapping_stats_d[,Used_read_headers,drop=F])>min_used_readNum,1]); report_samples

counts_names=paste("Num",report_samples, sep="_")
Length_names=paste("Len",report_samples, sep="_")
mismatch_names=paste("Mismatch",report_samples, sep="_")
PercCount_names=paste("PercCount",report_samples, sep="_")
PercLen_names=paste("PercLen",report_samples, sep="_")
PercMismatch_names=paste("PercMM",report_samples, sep="_")

report_sampleInfo_ids=match(report_samples,sample_info_d$sample); report_sampleInfo_ids
report_sampleMap_ids=match(report_samples,mapping_stats_d$Sample); report_sampleMap_ids
report_sample_ids=match(report_samples,samples); report_sample_ids
sampleConclusion_d=data.frame(Sample=report_samples
	,SampleInfo=paste0("For sample ",report_samples,", the GOI plasmid is ",sample_info_d$GOI[report_sampleInfo_ids],". The helper plasmid is ",sample_info_d$pHelper[report_sampleInfo_ids],". The RepCap plasmid is ",sample_info_d$pRepCap[report_sampleInfo_ids],". " )
	,ReadInfo=paste0("NGS generated ",formatC(mapping_stats_d$RawRead[report_sampleMap_ids], format="d", big.mark=",")," raw reads. After mapping, **", 
		apply(formatC(as.matrix(mapping_stats_d[report_sampleMap_ids, Used_read_headers, drop=F]),format="d", big.mark=","), 1, paste,collapse=" and "), "**",
		ifelse(length(Used_read_headers)==2, paste0( " (",paste(sub("Used.","",Used_read_headers),collapse=" and "),")" ),  ""), 
		" mapped properly and used for further analysis. ")
	,stringsAsFactors=F
)
sampleConclusion_d$TotalUsedReads=rowSums(as.matrix(mapping_stats_d[report_sampleMap_ids, Used_read_headers, drop=F]))


##4, stats of read mapping to different chromosome or vectors
mapping_detail_dlist=NULL

##load data to mapping_detail_dlist
mapping_detail_f=paste0(project_root,"/03_ana_readMap/readMap.stats.txt") #a summary for all files (usually they have the same GOI vector)
if(file.exists(mapping_detail_f)){
	mapping_detail_d=read.table(mapping_detail_f, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)
	#if_map2GOI=mapping_detail_d$region %in% "GOI_full"
	#GOI_vector_name=mapping_detail_d$reference[if_map2GOI]; GOI_vector_name
	#purity_summary_d=data.frame(Sample=report_samples, 'Estimated purity (%)'=unlist(mapping_detail_d[if_map2GOI, PercLen_names]),
	#  GOI_vector=GOI_vector_name, check.names=F )
	for(i in 1:length(report_samples)){
		mapping_detail_d1=mapping_detail_d[,c("reference","region","ref_size", counts_names[i], Length_names[i], PercCount_names[i], PercLen_names[i])]
		mapping_detail_dlist[[report_samples[i]]]=mapping_detail_d1
	}
}else{ #look for individual stats files, load to mapping_detail_dlist, 
	mapping_detail_fs=paste0(project_root,"/03_ana_readMap/",report_samples, ".readMap.stats.txt")
	for(i in 1:length(report_samples)){
		mapping_detail_f1=mapping_detail_fs[i]
		mapping_detail_d1=read.table(mapping_detail_f1, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)
		mapping_detail_dlist[[report_samples[i]]]=mapping_detail_d1
		if_map2GOI=mapping_detail_d1$region %in% "GOI_full"
		GOI_vector_name=mapping_detail_d1$reference[if_map2GOI]; GOI_vector_name
		mapping_detail_d1$reference[mapping_detail_d1$reference %in% GOI_vector_name]="GOI_vector" #
		mapping_detail_d1=mapping_detail_d1[,setdiff(names(mapping_detail_d1),"ref_size")]
		if(i==1){
			mapping_detail_d=mapping_detail_d1
		}else{
			mapping_detail_d=merge(mapping_detail_d, mapping_detail_d1, all=T)
		}
	}
	OutCombined_mapping_detail_f=paste0(project_root,"/03_ana_readMap/AllCombined.readMap.stats.txt")
	mapping_detail_d=re_format_tb(mapping_detail_d, back_headers=c(counts_names,PercCount_names, Length_names, PercLen_names,  mismatch_names, PercMismatch_names) )
	mapping_detail_d=mapping_detail_d[order(rowMeans(mapping_detail_d[,PercLen_names],na.rm=T),decreasing=T),]
	write.table(mapping_detail_d, file=OutCombined_mapping_detail_f, col.names=T, row.names=F, sep="\t", quote=F)
}

purity_summary_d=data.frame(Sample=report_samples,  stringsAsFactors=FALSE, check.names=F)
purity_summary_d[,report_impurity_headers]=NA
impurity_summary_d=data.frame(Sample=report_samples, 'Impurity1 (%)'=NA, 'Impurity1 source'=NA, 
	'Impurity2 (%)'=NA, 'Impurity2 source'=NA, 'Impurity3 (%)'=NA, 'Impurity3 source'=NA,  check.names=F)
for(i in 1:length(report_samples)){
	mapping_detail_d1=mapping_detail_dlist[[report_samples[i]]]
	if_map2GOI=mapping_detail_d1$region %in% "GOI_full"
	GOI_vector_name=mapping_detail_d1$reference[if_map2GOI]; GOI_vector_name
	mapping_detail_d1$region=sub("^KanR$|^NeoR$|^KanR.ORF$","NeoR/KanR",mapping_detail_d1$region) #WL added 2022/10/02 (change KanR, NeoR or KanR.ORF to "NeoR/KanR")

	report_impurity_references_i=report_impurity_references
	if_convert_ref_name=report_impurity_references %in% names(sample_info_d)
	report_impurity_references_i[if_convert_ref_name]=unlist(sample_info_d[match(report_samples[i], sample_info_d$sample), report_impurity_references[if_convert_ref_name]])
	
	matched_rows=match(paste(report_impurity_references_i,report_impurity_regions,sep=":"), paste(mapping_detail_d1$reference,mapping_detail_d1$region,sep=":"))
	if_matched_rows=!is.na(matched_rows)
	if_matched_region=is.na(matched_rows) & !(report_impurity_regions %in% c("All")) ; sum(if_matched_region)
	purity_summary_d[i,report_impurity_headers[if_matched_rows]]=mapping_detail_d1[matched_rows[if_matched_rows], PercLen_names[i]]
	purity_summary_d[i,report_impurity_headers[if_matched_region]] <- tapply(mapping_detail_d1[, PercLen_names[i]], mapping_detail_d1$region, sum, na.rm=T)[report_impurity_regions[if_matched_region]]

	if_exclude_4cal_impurity=if_map2GOI | mapping_detail_d1$reference %in% "All" | (mapping_detail_d1$reference %in% GOI_vector_name & mapping_detail_d1$region %in% "All")
	mapping_detail_d2=mapping_detail_d1[!if_exclude_4cal_impurity, ]
	mapping_detail_d2=mapping_detail_d2[order(mapping_detail_d2[,PercLen_names[i]], decreasing=T), ]
	impurity_summary_d[i,c(2,4,6)]=round(mapping_detail_d2[1:3, PercLen_names[i]],3)
	impurity_summary_d[i,c(3,5,7)]=paste(mapping_detail_d2$reference[1:3], mapping_detail_d2$region[1:3], sep=":")
}

sampleConclusion_d$PurityInfo=paste0(
	"The purity of the rAAV sample is defined by the percent of read coverage (ReadCoverage%) mapping inside the ITR-GOI-ITR region (desired full-length rAAV genome noted as \\“GOI_full\\” in the report).  "
	,"The estimated purity of this sample is **", round(purity_summary_d[,"ITR-GOI-ITR%"],2),"%**. "
	,"The estimated impurity from the backbone of the GOI plasmid is **", round(purity_summary_d[,"GOIplasmid%"]-purity_summary_d[,"ITR-GOI-ITR%"],2),"%** (",round(purity_summary_d[,"GOIplasmid%"],3),"% - ",round(purity_summary_d[,"ITR-GOI-ITR%"],3),"%). "
	,ifelse(is.na(purity_summary_d[,"pRepCap%"]),"",paste0("The estimated impurity from the RepCap plasmid is **",round(purity_summary_d[,"pRepCap%"],3),"%**. "))
	,ifelse(is.na(purity_summary_d[,"pHelper%"]),"",paste0("The estimated impurity from the helper plasmid is **",round(purity_summary_d[,"pHelper%"],3),"%**. "))
	,"The estimated impurity from the human chromosome is **",round(purity_summary_d[,"Human.Chromosomes%"],3),"%**. "
	,"The estimated impurity from the first 5kb of human adenovirus 5 (Ad5) genome (AC_000008:1-5000bp, contains E1A and E1B genes) is **",round(purity_summary_d[,"Human_Ad5_first5k%"],6),"%**. "
	,ifelse(is.na(purity_summary_d[,"NeoR/KanR%"]),"",paste0("The estimated impurity from the NeoR/KanR gene is **",round(purity_summary_d[,"NeoR/KanR%"],3),"%**. "))
	,ifelse(is.na(purity_summary_d[,"AmpR%"]),"",paste0("The estimated impurity from the AmpR gene is **",round(purity_summary_d[,"AmpR%"],3),"%**. "))
)


##5.1, read alignment coverage pictures
read_cov_dir=paste0(project_root, "/04_readCov_plot/")
read_cov_pic_fl=sapply(setdiff(list.dirs(read_cov_dir,full.names=F),""), function(folder1){list.files(paste0(read_cov_dir,folder1,"/"), ".readCoverage.png")}, simplify=F ); read_cov_pic_fl
read_cov_pic_fs=list.files(read_cov_dir, ".readCoverage.png", recursive=F); read_cov_pic_fs
if(!("All" %in% names(read_cov_pic_fl)) & length(read_cov_pic_fs)>0){
	read_cov_pic_fl$All=read_cov_pic_fs
}

for(sample1 in names(read_cov_pic_fl)){
	mkdir_if_not_exist(paste0(report_out_dir,"readCov_plot/",sample1,"/"))
	for(f1 in read_cov_pic_fl[[sample1]]){
		system( paste0("cp -f ", read_cov_dir,sample1,"/",f1, " ", report_out_dir,"readCov_plot/",sample1,"/",f1) )
	}
}

##5.2, chimeric read plot
read_cov_dir=paste0(project_root, "/04_readCov_plot/")
chimericRead_pic_fl=sapply(setdiff(list.dirs(read_cov_dir,full.names=F),""), function(folder1){list.files(paste0(read_cov_dir,folder1,"/"), ".chimeric.png")}, simplify=F ); chimericRead_pic_fl
# chimericRead_pic_fs=list.files(read_cov_dir, ".chimeric.png", recursive=F); chimericRead_pic_fs
# if(!("All" %in% names(chimericRead_pic_fl)) & length(chimericRead_pic_fs)>0){
# 	chimericRead_pic_fl$All=chimericRead_pic_fs
# }

for(sample1 in names(chimericRead_pic_fl)){
	mkdir_if_not_exist(paste0(report_out_dir,"readCov_plot/",sample1,"/"))
	for(f1 in chimericRead_pic_fl[[sample1]]){
		system( paste0("cp -f ", read_cov_dir,sample1,"/",f1, " ", report_out_dir,"readCov_plot/",sample1,"/",f1) )
	}
}

##5.3 chimeric read stats tables
chimericRead_stats_fs=paste0(project_root,"/02b_ChimericRead/",report_samples, ".chimericRead.detail.txt.stats.txt")
chimericRead_stats_dlist=NULL
for(i in 1:length(report_samples)){
	chimericRead_stats_f1=chimericRead_stats_fs[i]
	chimericRead_stats_d1=read.table(chimericRead_stats_f1, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)
	ifHighQuality_header=grep("ifHighQuality",names(chimericRead_stats_d1),value=T)
	if(length(ifHighQuality_header)==1){
		chimericRead_stats_d1=chimericRead_stats_d1[chimericRead_stats_d1[,ifHighQuality_header] %in% "Y",]
	}
	chimericRead_stats_d1=chimericRead_stats_d1[order(chimericRead_stats_d1$NumberOfMapping,decreasing=T),]
	#calculate percent of chimeric reads
	total_usedReadsCount=sampleConclusion_d$TotalUsedReads[match(report_samples[i], sampleConclusion_d$Sample)]
	chimericRead_stats_d1$PercChimericReads=chimericRead_stats_d1$NumberOfMapping/total_usedReadsCount*100
	chimericRead_stats_dlist[[report_samples[i]]]=chimericRead_stats_d1
	##save files to report folder
	out_f=paste0(report_out_dir, "ChimericReads/",report_samples[i],".chimericRead.stats.txt")
	mkdir_if_not_exist(out_f)
	write.table(chimericRead_stats_d1, file=out_f, col.names=T, row.names=F, sep="\t", quote=F)
}

##5.4 chimeric event tables
chimericEvent_fs=paste0(project_root,"/02b_ChimericRead/",report_samples, ".chimericReadCount.txt")
chimericEvent_dlist=NULL
for(i in 1:length(report_samples)){
	chimericEvent_f1=chimericEvent_fs[i]
	chimericEvent_d1=read.table(chimericEvent_f1, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)
	if_remove=grepl("Boundary",chimericEvent_d1$chimeric_type)
	chimericEvent_d1=chimericEvent_d1[!if_remove,]
	table(chimericEvent_d1$count)
	chimericEvent_dlist[[report_samples[i]]]=chimericEvent_d1
	##save files to report folder
	out_f=paste0(report_out_dir, "ChimericReads/",report_samples[i],".chimericEvent.txt")
	mkdir_if_not_exist(out_f)
	write.table(chimericEvent_d1, file=out_f, col.names=T, row.names=F, sep="\t", quote=F)
}


##6, read coverage length stats
readCovLen_stats_fs=paste0(project_root,"/03_ana_readMap/",report_samples, ".readCovLen.stats.txt")
readCovLen_stats_dlist=NULL
for(i in 1:length(report_samples)){
	readCovLen_stats_f1=readCovLen_stats_fs[i]
	readCovLen_stats_d1=read.table(readCovLen_stats_f1, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)
	names(readCovLen_stats_d1)=sub(paste0(report_samples[i], "."), "", names(readCovLen_stats_d1))
	readCovLen_stats_d1=readCovLen_stats_d1[order(readCovLen_stats_d1$median),]
	readCovLen_stats_dlist[[report_samples[i]]]=readCovLen_stats_d1
}


readCovLen_pic_fs=paste0(report_samples, ".readCovLen.boxplot.png"); readCovLen_pic_fs
mkdir_if_not_exist(paste0(report_out_dir,"readCovLen_boxplot/"))
for(f1 in readCovLen_pic_fs){
	system( paste0("cp -f ", project_root,"/03_ana_readMap/",f1, " ", report_out_dir,"readCovLen_boxplot/", f1) )
}


##7, read host cell DNA cluster annotation information
mkdir_if_not_exist(paste0(report_out_dir,"hostCellDNA_anno/"))

hostCellAnno_l=list()
for(sample in AAV_samples){
	peak_anno_f=paste0(peakAnno_out_prefix,".",sample,".xlsx")
	system( paste0("cp -f ", peak_anno_f, " ", report_out_dir,"hostCellDNA_anno/") )
	peak_anno_d=read.xlsx(peak_anno_f)
	hostCellAnno_l[[sample]]=peak_anno_d
}
peak_summary_f=paste0(peakAnno_out_prefix,".summary.xlsx")
peak_summary_d=read.xlsx(peak_summary_f)
peak_summary_d$Total=rowSums(peak_summary_d[,4:10],na.rm=T)
peak_summary_d_keys=paste(peak_summary_d$sample, peak_summary_d$GeneType, peak_summary_d$Count,sep=":")

peak_summary_d2=data.frame(Sample=report_samples
	,TotalUsedReads=sampleConclusion_d$TotalUsedReads
	,TotalHCReads=peak_summary_d[match(paste(report_samples, "any", "read#",sep=":"),peak_summary_d_keys), "Total"] 
)
peak_summary_d2$Percent_HCReads=peak_summary_d2$TotalHCReads/peak_summary_d2$TotalUsedReads *100
peak_summary_d2$HC_cluster=peak_summary_d[match(paste(report_samples, "any", "cluster#",sep=":"),peak_summary_d_keys), "Total"] 

peak_summary_d2$ReadMap2OncoGeneEx=rowSums(peak_summary_d[match(paste(report_samples, "any", "read#",sep=":"),peak_summary_d_keys), c("Promoter","5'.UTR","Exon")] )
peak_summary_d2$PercReadMap2OncoGeneEx=peak_summary_d2$ReadMap2OncoGeneEx/peak_summary_d2$TotalUsedReads *100
peak_summary_d2$ClusterMap2OncoGeneEx=rowSums(peak_summary_d[match(paste(report_samples, "any", "cluster#",sep=":"),peak_summary_d_keys), c("Promoter","5'.UTR","Exon")] )

sampleConclusion_d$HostCellDNAinfo=paste0(
	"  \n**Host cell DNA and oncogene analysis: ** Reads mapped to human genome were clustered. Reads mapped to the same chromosome overlapping with each other or has a gap<4kb were clustered. The read clusters were annotated with the nearest gene by using ChIPseeker (Yu, Wang, and He 2015), which reports the distance between the read cluster and the transcriptional start site (TSS) of the nearest annotated gene and also the location of the cluster relative to the gene (eg. promoter, 5' UTR,3' UTR, exon,  intron,  downstream, distal Intergenic). To further estimate the \\“oncogenic\\” potential of the host cell DNA, we used The Cancer Gene Census (CGC) annotation of human genes (Sondka et al. 2018) (https://cancer.sanger.ac.uk/census). There are 723 genes annotated in CGC (576 and 147 in tier 1 and 2 respectively) by their roles in cancer, eg. oncogene, tumor suppressor gene (TSG) or fusion, etc. "
	,"In total, ", formatC(peak_summary_d2$TotalHCReads, format="d", big.mark=","), " reads (",round(peak_summary_d2$Percent_HCReads,3),"%) map to human genome, forming ",formatC(peak_summary_d2$HC_cluster, format="d", big.mark=",")," clusters. "
	,"Among them, **", formatC(peak_summary_d2$ReadMap2OncoGeneEx, format="d", big.mark=","), " reads (",round(peak_summary_d2$PercReadMap2OncoGeneEx,3),"%**) map to promoter, 5’UTR or exonic region of an oncogene (annotated by CGC), forming ",formatC(peak_summary_d2$ClusterMap2OncoGeneEx, format="d", big.mark=",")," clusters. "
)

save.image(out_data_img_f)

