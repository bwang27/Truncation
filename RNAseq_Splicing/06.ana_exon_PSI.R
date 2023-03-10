##analyze the inclustion level and change of exons with GAgt-type 5'ss using exon-exon junction reads
##2016/03/09, can study all exons as well (derived from gene blocks)
##2016/03/17, can calculate p value based on fisher's exact test, do not use arbitory read number as cutoff
##2016/03/17, also include exon_a5ss and exon_a3ss intron_a5ss and intron_a3ss gene block types to the table
##2016/08/08 add what2do=study_ssScore_vs_basalPSI
##2018/08/08 changed the way to search junctions skip the exon, find all junctions in the transcript id cover the whole exon region
##2018/09/13 can set deltaPSI to 0 if P-value is greater than a cutoff (zeroDeltaPSI_Pval)
##2020/03/05 force maximum max_skip_juns smallest junctions skip the exon of interest (default max_skip_juns=20), 
##2021/05/24 for Fisher's exact test, consider calculating p-value for each replicate of test and control group. final p-value can be summarized using the geo-mean of all P-values
##2021/06/03 able to calculate P-value adjusted delta-PSI (eg. cal_PadjDeltaPSI=T; PadjDeltaPSI_ps=c("Pttest","Pfisher"); PadjDeltaPSI_bs=c(-0.5,-0.2); PadjDeltaPSI_es=c(1e-6,1e-30); PadjDeltaPSI_cutoff=10 )
##2021/07/21, added function to output a table of top splicing events (topGeneList_tb) for each treatment
##2021/11/09-10, fixed a bug when sort the top splicing events; output more details in top splicing list (eg. rank, deltaPSI etc); also can output a list for non-adjusted deltaPSI method

##available what2do: all, cal_PSI, calReguType
args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../SharedCodes/Rfunc.inc.R") 
source("05.inc.func.R")

basal_psi_cuts=c(0,20,80,100)
NTs=c("A","C","G","T")
diNTs=apply(expand.grid(NTs,NTs),1,paste,collapse=""); diNTs

geno="hg19";
gene_model="refflat";
transcript_id_name=ifelse(is.na(args_v["transcript_id_name"]) , "refseqid" , args_v["transcript_id_name"]);
study_exon_type="all_exon" #GAgt_annotated_exon all_exon annotated_exon


use_p_type="Pttest"; P_cut=0.05
if(!is.na(args_v["use_p_type"])){ use_p_type= unlist ( strsplit(args_v["use_p_type"]," "))   } #can use >1 P values eg. c("Pttest","Pfisher")
if(!is.na(args_v["P_cut"])){ P_cut=as.numeric( unlist ( strsplit(args_v["P_cut"]," ") ) )  }
Default_zeroDeltaPSI_Pvals=c(Pttest=0.2, Pfisher=0.05, Padj=0.05)

zeroDeltaPSI_Pval=rep(0.1, length(use_p_type))
if(all(use_p_type %in% names(Default_zeroDeltaPSI_Pvals))){ 
	zeroDeltaPSI_Pval=Default_zeroDeltaPSI_Pvals[use_p_type] 
}
if(!is.na(args_v["zeroDeltaPSI_Pval"]) ){
	zeroDeltaPSI_Pval=as.numeric( unlist ( strsplit(args_v["zeroDeltaPSI_Pval"]," ") ) )
}


workers <- ifelse(is.na(args_v["workers"]), 8, as.numeric(args_v["workers"])); 
delta_PSI_cut=ifelse(is.na(args_v["delta_PSI_cut"]), 5, as.numeric(args_v["delta_PSI_cut"])); 
PSI_cal_rnumMin=ifelse(is.na(args_v["PSI_cal_rnumMin"]), 10, as.numeric(args_v["PSI_cal_rnumMin"])); 
max_skip_juns <- ifelse(is.na(args_v["max_skip_juns"]), 20, as.numeric(args_v["max_skip_juns"])); 
#what2do<- ifelse(is.na(args_v["what2do"]), "all", args_v["what2do"])
if(is.na(args_v["what2do"])){
	what2do="all"
}else{
	what2do=unlist ( strsplit(args_v["what2do"]," ") )
}


###project setting (example)
study_name="2016_RO247_RO067_hs"
study_exon_type="all_exon"

sample_repl_l<- string2list("PNNDMSO.pA:PNNDMSO.pA1 PNNDMSO.pA2;PNN247300.pA:PNN247300.pA1 PNN247300.pA2;PNN2473uM.pA:PNN2473uM.pA1 PNN2473uM.pA2;PNN06730.pA:PNN06730.pA1 PNN06730.pA2;PNN0671uM.pA:PNN0671uM.pA1 PNN0671uM.pA2;HDFDMSO.pA:HDFDMSO.pA1 HDFDMSO.pA2;HDF247300.pA:HDF247300.pA1 HDF247300.pA2;HDF2473uM.pA:HDF2473uM.pA1 HDF2473uM.pA2;HDF06730.pA:HDF06730.pA1 HDF06730.pA2;HDF0671uM.pA:HDF0671uM.pA1 HDF0671uM.pA2;HDFDMSO.rz:HDFDMSO.rz1 HDFDMSO.rz2 HDFDMSO.rz3;HDF2473uM.rz:HDF2473uM.rz1 HDF2473uM.rz2 HDF2473uM.rz3")
sample_compare_matrix=rbind(
	unlist(strsplit("PNN247300.pA PNN2473uM.pA PNN06730.pA PNN0671uM.pA HDF247300.pA HDF2473uM.pA HDF06730.pA HDF0671uM.pA HDF2473uM.rz"," ")), 
	unlist(strsplit("PNNDMSO.pA PNNDMSO.pA PNNDMSO.pA PNNDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.rz"," ")) )
all_samples=unlist(strsplit("PNNDMSO.pA1 PNNDMSO.pA2 PNN247300.pA1 PNN247300.pA2 PNN2473uM.pA1 PNN2473uM.pA2 PNN06730.pA1 PNN06730.pA2 PNN0671uM.pA1 PNN0671uM.pA2 HDFDMSO.pA1 HDFDMSO.pA2 HDF247300.pA1 HDF247300.pA2 HDF2473uM.pA1 HDF2473uM.pA2 HDF06730.pA1 HDF06730.pA2 HDF0671uM.pA1 HDF0671uM.pA2 HDFDMSO.rz1 HDFDMSO.rz2 HDFDMSO.rz3 HDF2473uM.rz1 HDF2473uM.rz2 HDF2473uM.rz3"," "))
basal_PSI_samples=c("PNNDMSO.pA","HDFDMSO.pA","HDFDMSO.rz")

#parameters from command line 
study_name<- ifelse(is.na(args_v["study_name"]), "", args_v["study_name"])
geno<- ifelse(is.na(args_v["geno"]), "hg19", args_v["geno"])
study_exon_type<- ifelse(is.na(args_v["study_exon_type"]), "all_exon", args_v["study_exon_type"])
sample_compare_matrix=rbind(unlist(strsplit(as.character(args_v["test_samples"])," ")), unlist(strsplit(as.character(args_v["ref_samples"])," ")))
if(!is.na(args_v["sample_repl_l"])){ sample_repl_l<- string2list(args_v["sample_repl_l"]) }
if(!is.na(args_v["all_samples"])){ all_samples<- unlist(strsplit(as.character(args_v["all_samples"])," ")) }
if(!is.na(args_v["basal_PSI_samples"])){ basal_PSI_samples<- unlist(strsplit(as.character(args_v["basal_PSI_samples"])," ")) }
avg_samples=names(sample_repl_l)
gblock_f=paste("05.DEXSeq/geneBlock/combine_d.tbl", sep="")
if(!is.na(args_v["gblock_f"])){ gblock_f<- args_v["gblock_f"] }

read_num_f <- paste("01.comb_junc_map_info","/combine.junc2gene.tbl",sep="")
if(!is.na(args_v["read_num_f"])){ read_num_f<- args_v["read_num_f"] }

exon_defi_f=NA #  paste("../ss_annotation/02GAgt_ss5_exons/",geno,".",gene_model,".GAgt.exons.txt",sep="")
if(study_exon_type=="annotated_exon"){exon_defi_f=paste("../ReferenceDB/Splicing/03.Annotated_exons/",geno,".",gene_model,".Annotated.exons.txt",sep="")}
if(!is.na(args_v["exon_defi_f"])){ exon_defi_f<- args_v["exon_defi_f"] }

out_root=paste("06.exon_PSI/",study_name,"/",sep="")
if(!is.na(args_v["out_root"])){ out_root<- args_v["out_root"] }



##settings for cal_PadjDeltaPSI;
cal_PadjDeltaPSI<- ifelse(is.na(args_v["cal_PadjDeltaPSI"]), F, args_v["cal_PadjDeltaPSI"] %in% c("T","1","True") )

PadjDeltaPSI_ps=c("Pttest","Pfisher"); PadjDeltaPSI_bs=c(-0.5,-0.2); PadjDeltaPSI_es=c(1e-6,1e-30); PadjDeltaPSI_cutoff=10
if(!is.na(args_v["PadjDeltaPSI_ps"])){ PadjDeltaPSI_ps= unlist ( strsplit(args_v["PadjDeltaPSI_ps"]," "))   } #can use >1 P values eg. c("Pttest","Pfisher")
if(!is.na(args_v["PadjDeltaPSI_bs"])){ PadjDeltaPSI_bs=as.numeric( unlist ( strsplit(args_v["PadjDeltaPSI_bs"]," ") ) )  }
if(!is.na(args_v["PadjDeltaPSI_es"])){ PadjDeltaPSI_es=as.numeric( unlist ( strsplit(args_v["PadjDeltaPSI_es"]," ") ) )  }
PadjDeltaPSI_cutoff=ifelse(is.na(args_v["PadjDeltaPSI_cutoff"]), delta_PSI_cut, as.numeric(args_v["PadjDeltaPSI_cutoff"])); 
if(cal_PadjDeltaPSI){
	delta_PSI_cut=PadjDeltaPSI_cutoff
}

max_geneList_num=ifelse(is.na(args_v["max_geneList_num"]), 100, as.numeric(args_v["max_geneList_num"]) )


species=spe_common_names[as.character(tax_ids[geno])]
spe_latin_name=spe_latin_names[as.character(tax_ids[species])]
suppNum_names=paste("num_",all_samples,sep="")
IncUpsNum_names=paste("IncUpsNum_",all_samples,sep="")
IncDnsNum_names=paste("IncDnsNum_",all_samples,sep="")
ExcNum_names=paste("ExcNum_",all_samples,sep="")

PSI_names=paste("PSI_",all_samples,sep="")
avgPSI_names=paste("PSI_",avg_samples,sep="")
sample_pairs=paste(sample_compare_matrix[1,],sample_compare_matrix[2,],sep="_")
delta_PSI_names=paste("deltaPSI_", sample_pairs,sep="")
adjDeltaPSI_names=sub("deltaPSI_","AdjDeltaPSI_",delta_PSI_names); 
fisher_pval_names=paste("Pfisher_", sample_pairs,sep="")
ttest_pval_names=paste("Pttest_", sample_pairs,sep="")
adjpval_names=paste("Padj_", sample_pairs,sep="")
log2R_names=paste("Log2R_", sample_pairs,sep="")
ReguType_names=paste("ReguType_", sample_pairs,sep="")
basal_PSI_names=paste("PSI_",basal_PSI_samples,sep="")

PadjDeltaPSI_pName_matrix=matrix(paste(rep(PadjDeltaPSI_ps,each=length(sample_pairs)),sample_pairs,sep="_"), nrow=length(sample_pairs) )

ReguTypes=c("UP","DN","NC","na")
out_tb_f=paste(out_root,"/",study_exon_type, "/exons.tbl", sep="")
out_img_f=paste(out_root,"/",study_exon_type,"/exons.tbl.img", sep="")
out_param_img_f=paste(out_root,"/",study_exon_type,"/exons.tbl.parameters.img", sep="")


####RUN
define_parallel_fun(nCores=workers)
mkdir_if_not_exist(out_tb_f)
save.image(out_param_img_f)

##1, add read number to exon table
if(any(what2do %in% c("all")) ){
	print(paste("read", read_num_f))
	read_num_d=read.table(read_num_f, header=T, sep="\t", stringsAsFactors=F, quote="")
	read_num_d$intron_id=paste(read_num_d$contig, read_num_d$strand, read_num_d$juncpos5, read_num_d$juncpos3, sep=":")
	transcript2ss5_l=tapply(read_num_d$juncpos5,read_num_d[,transcript_id_name],function(v){v})
	transcript2ss3_l=tapply(read_num_d$juncpos3,read_num_d[,transcript_id_name],function(v){v})
	read_num_d=read_num_d[!duplicated(read_num_d$intron_id),]
	read_num_d$ss5_id=paste(read_num_d$contig, read_num_d$strand, read_num_d$juncpos5,  sep=":")
	read_num_d$ss3_id=paste(read_num_d$contig, read_num_d$strand, read_num_d$juncpos3,  sep=":")
}

if(any(what2do %in% c("all")) & study_exon_type %in% c("GAgt_annotated_exon","annotated_exon") ){
	print(paste("read", exon_defi_f))
	exon_defi_d=read.table(exon_defi_f, header=T, sep="\t", stringsAsFactors=F, quote="") #refseq defined GA-type exons
	exon_defi_d[,IncUpsNum_names]=read_num_d[match(exon_defi_d$ups_intron,read_num_d$intron_id), suppNum_names]
	exon_defi_d[,IncDnsNum_names]=read_num_d[match(exon_defi_d$dns_intron,read_num_d$intron_id), suppNum_names]
	exon_defi_d[, ExcNum_names]=read_num_d[match(exon_defi_d$skip_intron,read_num_d$intron_id), suppNum_names];
	exon_defi_d[,IncUpsNum_names][is.na(exon_defi_d[,IncUpsNum_names])]=0
	exon_defi_d[,IncDnsNum_names][is.na(exon_defi_d[,IncDnsNum_names])]=0
	exon_defi_d[, ExcNum_names][is.na(exon_defi_d[, ExcNum_names])]=0
	rm(read_num_d)
	substr(exon_defi_d$ss5_seq,11,20)=tolower(substr(exon_defi_d$ss5_seq,11,20))
	substr(exon_defi_d$ss3_seq,1,20)=tolower(substr(exon_defi_d$ss3_seq,1,20))
	names(exon_defi_d)=change_values(names(exon_defi_d), c("ss5_seq","ss3_seq"), c("endSS_seq","startSS_seq") )
	if(!("coordinates" %in% names(exon_defi_d))){
		exon_defi_d$coordinates=paste(exon_defi_d$contig, ":", ifelse(exon_defi_d$strand=="+", exon_defi_d$ups_ss3, exon_defi_d$ss5), "-", 
			ifelse(exon_defi_d$strand=="-", exon_defi_d$ups_ss3, exon_defi_d$ss5), sep="" )
	}
	exon_defi_d$region_ano="exon"
}
	

##or 2, all non-refseq defined exons
if(any(what2do %in% c("all")) & study_exon_type=="all_exon" ){
	print(paste("read", gblock_f))
	gblock_d=read.table(gblock_f, header=T, sep="\t", stringsAsFactors=F, quote="")
	exon_defi_d=gblock_d[gblock_d$region_ano %in% c("exon","exon_a3ss","exon_a5ss","intron_a3ss","intron_a5ss","exon_part"),]; nrow(gblock_d)
	rm(gblock_d)
	#if_annotated=exon_defi_d$startSS_supp=="Refseq" & exon_defi_d$endSS_supp=="Refseq"
	#exon_defi_d=exon_defi_d[!if_annotated, ]
	#find ups_intron dns_intron and skip_intron
	#if_remove=grepl("cs|ce|ts|te",exon_defi_d$startPosAno) | grepl("cs|ce|ts|te",exon_defi_d$endPosAno); table(if_remove)
	#exon_defi_d=exon_defi_d[!if_remove, ]
	if_start_ss5=!is.na(exon_defi_d$startPosType) & exon_defi_d$startPosType=="ss5"
	if_end_ss3=!is.na(exon_defi_d$endPosType) & exon_defi_d$endPosType=="ss3"
	strand_signs=ifelse(exon_defi_d$strand=="-",-1,1)
	exon_defi_d$endPos_id=paste(exon_defi_d$contig, exon_defi_d$strand, exon_defi_d$end_pos+ifelse(if_end_ss3,strand_signs,0),  sep=":")
	exon_defi_d$startPos_id=paste(exon_defi_d$contig, exon_defi_d$strand, exon_defi_d$start_pos-ifelse(if_start_ss5,strand_signs,0),  sep=":")
	exon_defi_d$endPos_id[!(exon_defi_d$endPosType %in% c("ss5","ss3"))]=NA
	exon_defi_d$startPos_id[!(exon_defi_d$startPosType %in% c("ss5","ss3"))]=NA
	if_non_ss=is.na(exon_defi_d$endPos_id) & is.na(exon_defi_d$startPos_id); table(if_non_ss); #no splice site
	exon_defi_d=exon_defi_d[!if_non_ss,]

	exon_defi_d$region_id=paste(exon_defi_d$contig, exon_defi_d$strand, exon_defi_d$start_pos,exon_defi_d$end_pos,  sep=":")
	if(!("coordinates" %in% names(exon_defi_d))){
		exon_defi_d$coordinates=paste(exon_defi_d$contig, ":", ifelse(exon_defi_d$strand=="+", exon_defi_d$start_pos, exon_defi_d$end_pos), "-", 
			ifelse(exon_defi_d$strand=="-", exon_defi_d$start_pos, exon_defi_d$end_pos), sep="" )
	}
	exon_defi_d=re_format_tb(exon_defi_d, del_pattern="^P_|ReguType_" )

	##find all upstream 5'ss and all downstream 3'ss
	ss3toss5_l=tapply(read_num_d$juncpos5,read_num_d$ss3_id, function(v){v[!duplicated(v)]})
	ss5toss3_l=tapply(read_num_d$juncpos3,read_num_d$ss5_id, function(v){v[!duplicated(v)]})
	ss3toJuncNum=sapply(ss3toss5_l,length)
	ss5toJuncNum=sapply(ss5toss3_l,length)
	ss3toJuncReadNum_d=data.frame(ss3_id=names(ss3toss5_l))
	ss3toJuncReadNum_d[,suppNum_names]=myApply(suppNum_names, function(suppNum_name){
		print(suppNum_name)
		tapply(read_num_d[,suppNum_name], read_num_d$ss3_id, sum, na.rm=T)[as.character(ss3toJuncReadNum_d$ss3_id)]
	})
	ss5toJuncReadNum_d=data.frame(ss5_id=names(ss5toss3_l))
	ss5toJuncReadNum_d[,suppNum_names]=myApply(suppNum_names, function(suppNum_name){
		print(suppNum_name)
		tapply(read_num_d[,suppNum_name], read_num_d$ss5_id, sum, na.rm=T)[as.character(ss5toJuncReadNum_d$ss5_id)]
	})

	if_start_ss5=!is.na(exon_defi_d$startPosType) & exon_defi_d$startPosType=="ss5"
	if_end_ss3=!is.na(exon_defi_d$endPosType) & exon_defi_d$endPosType=="ss3"
	if_start_ss3=!is.na(exon_defi_d$startPosType) & exon_defi_d$startPosType=="ss3"
	if_end_ss5=!is.na(exon_defi_d$endPosType) & exon_defi_d$endPosType=="ss5"
	exon_defi_d$ups_junc_num=NA; 
	exon_defi_d$dns_junc_num=NA; 
	exon_defi_d$skip_junc_num=NA
	exon_defi_d$ups_junc_num[if_start_ss3]=ss3toJuncNum[as.character(exon_defi_d$startPos_id[if_start_ss3])]
	exon_defi_d$dns_junc_num[if_end_ss5]=ss5toJuncNum[as.character(exon_defi_d$endPos_id[if_end_ss5])]
	exon_defi_d$skip_junc_num[if_start_ss5]=ss5toJuncNum[as.character(exon_defi_d$startPos_id[if_start_ss5])] 
	exon_defi_d$skip_junc_num[if_end_ss3]=ss3toJuncNum[as.character(exon_defi_d$endPos_id[if_end_ss3])] 

	exon_defi_d[,IncUpsNum_names]=NA; 
	exon_defi_d[,IncDnsNum_names]=NA;
	exon_defi_d[, ExcNum_names]=NA
	exon_defi_d[if_start_ss3,IncUpsNum_names]=ss3toJuncReadNum_d[match(exon_defi_d$startPos_id[if_start_ss3], ss3toJuncReadNum_d$ss3_id),suppNum_names]
	exon_defi_d[if_end_ss5,IncDnsNum_names]=ss5toJuncReadNum_d[match(exon_defi_d$endPos_id[if_end_ss5], ss5toJuncReadNum_d$ss5_id),suppNum_names]
	#update ExcNum_names for a5ss and a3ss
	exon_defi_d[if_start_ss5, ExcNum_names]=ss5toJuncReadNum_d[match(exon_defi_d$startPos_id[if_start_ss5], ss5toJuncReadNum_d$ss5_id),suppNum_names]
	exon_defi_d[if_end_ss3, ExcNum_names]=ss3toJuncReadNum_d[match(exon_defi_d$endPos_id[if_end_ss3], ss3toJuncReadNum_d$ss3_id),suppNum_names]
	rm(ss5toJuncReadNum_d,ss3toJuncReadNum_d)

	##find junctions exclude the exon of interest
	if_exon=exon_defi_d$region_ano=="exon"
	exon_d=exon_defi_d[!duplicated(exon_defi_d$region_id) & if_exon, c("contig","strand", transcript_id_name, "region_id","endPos_id","startPos_id","start_pos","end_pos")]
	exonStrand_signs=ifelse(exon_d$strand=="+",1,-1)
	exon_start_coors=exon_d$start_pos
	exon_end_coors=exon_d$end_pos

	exon_ss3toss5_l=ss3toss5_l[as.character(exon_d$startPos_id)]
	exon_ss5toss3_l=ss5toss3_l[as.character(exon_d$endPos_id)]
	# ups_ss5s_idlist=sapply(1:nrow(exon_d), function(i){
	#  if(is.null(exon_ss3toss5_l[[i]])){
	#  		return(NULL)
	#  	}else{
	#  	return( paste( exon_d$contig[i],exon_d$strand[i],exon_ss3toss5_l[[i]],sep=":")) }
	# })
	# dns_ss3s_idlist=sapply(1:nrow(exon_d), function(i){
	#  if(is.null(exon_ss5toss3_l[[i]])){
	#  		return(NULL)
	#  	}else{
	#  	return( paste( exon_d$contig[i],exon_d$strand[i],exon_ss5toss3_l[[i]],sep=":")) }
	# })

	# exon_to_skipJunc_l=myApply(1:nrow(exon_d), function(i){
	# 	as.vector(outer(exon_ss3toss5_l[[i]],  exon_ss5toss3_l[[i]], paste, sep=":"))
	# })
	exon2gene2allss5_l=transcript2ss5_l[exon_d[,transcript_id_name]]
	exon2gene2allss3_l=transcript2ss3_l[exon_d[,transcript_id_name]]
	exon_to_skipJunc_l=myApply(1:nrow(exon_d), function(i){ #2018/05/10
		if(i %% 5000==0){print (i)}
		#ups_ss5s=ups_ss5s_idlist[[i]]
		#dns_ss3s=dns_ss3s_idlist[[i]]
		#if(is.null(ups_ss5s) | is.null(dns_ss3s)){return(NULL)}
		#upsss5toss3s=unlist(ss5toss3_l[as.character(ups_ss5s)])
		#dnsss3toss5s=unlist(ss3toss5_l[as.character(dns_ss3s)])
		#upsss5toss3s=upsss5toss3s[ (upsss5toss3s*exonStrand_signs[i])> (exon_end_coors[i]*exonStrand_signs[i]) ]
		#dnsss3toss5s=dnsss3toss5s[ (dnsss3toss5s*exonStrand_signs[i])< (exon_start_coors[i]*exonStrand_signs[i]) ]
		#unique(c(as.vector(outer(exon_ss3toss5_l[[i]],  upsss5toss3s, paste, sep=":")), as.vector(outer(dnsss3toss5s,  exon_ss5toss3_l[[i]], paste, sep=":"))))

		#read_num_d$contig, read_num_d$strand, read_num_d$juncpos5, read_num_d$juncpos3
		#if_skip_intron=read_num_d$contig %in% exon_d$contig[i] & read_num_d$strand %in% exon_d$strand[i] & 
		#	read_num_d$juncpos5*exonStrand_signs[i] < exon_start_coors[i]*exonStrand_signs[i] &
		#	read_num_d$juncpos3*exonStrand_signs[i] > exon_end_coors[i]*exonStrand_signs[i]
		#if_skip_intron[is.na(if_skip_intron)]=F
		
		geness5s=exon2gene2allss5_l[[i]] #2018/8/8
		geness3s=exon2gene2allss3_l[[i]]
		if(is.null(geness5s) | is.null(geness3s)){return(NULL)}
		if_skip_intron=geness5s*exonStrand_signs[i] < exon_start_coors[i]*exonStrand_signs[i] &
			geness3s*exonStrand_signs[i] > exon_end_coors[i]*exonStrand_signs[i]
		if_skip_intron[is.na(if_skip_intron)]=F
		if(any(if_skip_intron)){
			skip_intron_d=data.frame(ss5=geness5s[if_skip_intron], ss3=geness3s[if_skip_intron])
			skip_intron_d$span=abs(skip_intron_d$ss5-skip_intron_d$ss3)
			skip_intron_d=skip_intron_d[order(skip_intron_d$span), ]
			skip_intron_d$skip_juncs=paste(skip_intron_d$ss5, skip_intron_d$ss3,sep=":")
			skip_intron_d=skip_intron_d[!duplicated(skip_intron_d$skip_juncs),]
			if(nrow(skip_intron_d)>max_skip_juns){
				skip_intron_d=skip_intron_d[1:max_skip_juns,] ##2020/03/05 force maximum 5 smallest junctions skip the exon of intrest
			}
			return(skip_intron_d$skip_juncs)
		}else{
			return(NULL)
		}
	})
	rm(exonStrand_signs,exon_start_coors,exon_end_coors)
	exonSkip_rnum_d=data.frame(region_id=rep(exon_d$region_id, sapply(exon_to_skipJunc_l,length) ))
	exonSkip_rnum_d[,c("contig","strand")]=exon_d[match(exonSkip_rnum_d$region_id, exon_d$region_id), c("contig","strand")]
	exonSkip_rnum_d$intron_id=paste(exonSkip_rnum_d$contig,exonSkip_rnum_d$strand,  unlist(exon_to_skipJunc_l), sep=":") 
	exonSkip_rnum_d[, suppNum_names]=read_num_d[match(exonSkip_rnum_d$intron_id, read_num_d$intron_id),suppNum_names]
	exonSkip_rnum_d[, suppNum_names][is.na(exonSkip_rnum_d[, suppNum_names])]=0
	exonSkipJuncNum_v=tapply(exonSkip_rnum_d$intron_id %in% read_num_d$intron_id, exonSkip_rnum_d$region_id, sum)
	exon_defi_d$skip_junc_num[if_exon]=exonSkipJuncNum_v[ as.character(exon_defi_d$region_id[if_exon]) ]
	table(exon_defi_d$skip_junc_num==0, exon_defi_d$region_ano)
	exon_d[,suppNum_names]=myApply(suppNum_names, function(suppNum_name){
		print(suppNum_name)
		tapply(exonSkip_rnum_d[,suppNum_name], exonSkip_rnum_d$region_id, sum, na.rm=T)[as.character(exon_d$region_id)]
	})

	
	exon_defi_d[if_exon, ExcNum_names]=exon_d[match(exon_defi_d$region_id[if_exon], exon_d$region_id),suppNum_names]
	table(exon_defi_d$skip_junc_num==0 | is.na(exon_defi_d$skip_junc_num), exon_defi_d$region_ano)
	if_no_skipping_event=exon_defi_d$skip_junc_num==0 | is.na(exon_defi_d$skip_junc_num)
	exon_defi_d=exon_defi_d[!if_no_skipping_event,]; nrow(exon_defi_d)
	if_with_inc_events=(!is.na(exon_defi_d$ups_junc_num) & exon_defi_d$ups_junc_num>0) |  (!is.na(exon_defi_d$dns_junc_num) & exon_defi_d$dns_junc_num>0); table(if_with_inc_events)
	exon_defi_d=exon_defi_d[if_with_inc_events, ]; nrow(exon_defi_d)
	rm(exon_d, read_num_d, ss3toss5_l, ss5toss3_l, exon_ss3toss5_l, exon_ss5toss3_l, exon_to_skipJunc_l, exonSkip_rnum_d)
	exon_defi_d[, IncUpsNum_names][is.na(exon_defi_d[, IncUpsNum_names])]=0
	exon_defi_d[, IncDnsNum_names][is.na(exon_defi_d[, IncDnsNum_names])]=0
	exon_defi_d[, ExcNum_names][is.na(exon_defi_d[, ExcNum_names])]=0
	
	print(paste("output raw data to",out_tb_f))
	write.table(exon_defi_d, file=out_tb_f, sep="\t", col.names=T, row.names=F, quote=F)
	save.image(out_img_f)
}




##3,  calculate PSI (percent-included-in values)
if( ! all(what2do %in% c("all")  ) ){
	print(paste("re-load",out_tb_f))
	exon_defi_d=read.table(out_tb_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
}

if(any(what2do %in% c("all","cal_PSI"))  ){
	print(paste("calculate PSI..."))
	IncSumNum_tb=exon_defi_d[,IncUpsNum_names]+exon_defi_d[,IncDnsNum_names]
	IncNum_tb=IncSumNum_tb
	IncNum_tb[exon_defi_d$region_ano=="exon",]=round(IncNum_tb[exon_defi_d$region_ano=="exon",]/2,1)
	allPSI_tb=IncNum_tb/(IncNum_tb+exon_defi_d[, ExcNum_names])*100
	allPSI_tb[(IncSumNum_tb+exon_defi_d[, ExcNum_names])<PSI_cal_rnumMin]=NA
	ExcSumNum_tb=exon_defi_d[, ExcNum_names]
	colnames(allPSI_tb)=all_samples
	colnames(IncSumNum_tb)=all_samples
	colnames(ExcSumNum_tb)=all_samples
	#calculate average PSI, 
	avgPSI_tb=sapply(avg_samples,function(avg_sample1){
		if(avg_sample1 %in% names(sample_repl_l)){
			if(length(sample_repl_l[[avg_sample1]])==1){
				allPSI_tb[, sample_repl_l[[avg_sample1]]  ]
			}else{
				rowMeans(allPSI_tb[, sample_repl_l[[avg_sample1]]], na.rm=T)
			}
			
		}else{
			allPSI_tb[,avg_sample1]
		}
	})
	exon_defi_d[PSI_names]=round(allPSI_tb,3)
	exon_defi_d[avgPSI_names]=round(avgPSI_tb,3)
	#table(is.na(exon_defi_d$PSI_DMSO))
	exon_defi_d[delta_PSI_names]=round(avgPSI_tb[, sample_compare_matrix[1,]] - avgPSI_tb[, sample_compare_matrix[2,]], 1)
	
	#calculate p value using fisher's exact test
	for(i in 1:ncol(sample_compare_matrix)){
		print(paste("calculate",fisher_pval_names[i]))
		if(sample_compare_matrix[1,i] %in% names(sample_repl_l)){
			sample_1s=sample_repl_l[[sample_compare_matrix[1,i] ]]
		}else{
			sample_1s=sample_compare_matrix[1,i]
		}
		if(sample_compare_matrix[2,i] %in% names(sample_repl_l)){
			sample_2s=sample_repl_l[[sample_compare_matrix[2,i] ]]
		}else{
			sample_2s=sample_compare_matrix[2,i]
		}
		fisher_slogPval_m=NULL
		for(sample_1 in sample_1s){
			for(sample_2 in sample_2s){
				tmp_tb=cbind(IncSumNum_tb[sample_1],  ExcSumNum_tb[sample_1], IncSumNum_tb[sample_2],  ExcSumNum_tb[sample_2] )
				fisher_slogPval1=unlist( myApply(1:nrow(tmp_tb), function(row_i){  do_fisher_test(unlist(tmp_tb[row_i,]))  } ) )
				fisher_slogPval1[abs(fisher_slogPval1)==Inf]= sign(fisher_slogPval1[abs(fisher_slogPval1)==Inf])*100
				fisher_slogPval_m=cbind(fisher_slogPval_m, fisher_slogPval1)
			}
		}
		
		fisher_slogPval_m[abs(fisher_slogPval_m)==Inf]= sign(fisher_slogPval_m[abs(fisher_slogPval_m)==Inf])*100
		fisher_slogP_means=rowMeans(fisher_slogPval_m[,,drop=F])
		exon_defi_d[,fisher_pval_names[i]]=10^(-abs(fisher_slogP_means)) ##geoMean of P-value
		#tmp_tb=cbind(rowSums(IncSumNum_tb[sample_1s]),  rowSums(ExcSumNum_tb[sample_1s]), rowSums(IncSumNum_tb[sample_2s]),  rowSums(ExcSumNum_tb[sample_2s]) )
		#exon_defi_d[,fisher_pval_names[i]]=unlist( myApply(1:nrow(tmp_tb), function(row_i){  fisher.test(matrix(tmp_tb[row_i,],2))$p.value } ) )


	}

	#calculate p-value using t-test
	for(i in 1:ncol(sample_compare_matrix)){
		print(paste("calculate",ttest_pval_names[i]))
		if(sample_compare_matrix[1,i] %in% names(sample_repl_l)){
			sample_1s=sample_repl_l[[sample_compare_matrix[1,i] ]]
		}else{
			sample_1s=sample_compare_matrix[1,i]
		}
		if(sample_compare_matrix[2,i] %in% names(sample_repl_l)){
			sample_2s=sample_repl_l[[sample_compare_matrix[2,i] ]]
		}else{
			sample_2s=sample_compare_matrix[2,i]
		}
		exon_defi_d[,ttest_pval_names[i]]=NA
		if(length(sample_1s)>0 & length(sample_2s)>0 & length(sample_1s)+length(sample_2s)>=3 ){
			tmp_tb1=allPSI_tb[,sample_1s,drop=F]
			tmp_tb2=allPSI_tb[,sample_2s,drop=F]
			comp_nums1=apply(tmp_tb1,1,function(v){ length(unique(v[!is.na(v)]) ) })
			comp_nums2=apply(tmp_tb2,1,function(v){ length(unique(v[!is.na(v)]) ) })
			if_do_ttest=comp_nums1>0 & comp_nums2>0 & comp_nums1+comp_nums2>=3
			exon_defi_d[if_do_ttest,ttest_pval_names[i]]=unlist( myApply( (1:nrow(tmp_tb1))[if_do_ttest], function(row_i){  
				v1=tmp_tb1[row_i,]; v2=tmp_tb2[row_i,]
				v1=v1[!is.na(v1)]; v2=v2[!is.na(v2)]
				pval=NA
				if(length(v1)==1){
					pval=t.test( v2, mu=v1 )$p.value
				}else if(length(v2)==1){
					pval=t.test( v1, mu=v2 )$p.value
				}else{
					pval=t.test( v1, v2 )$p.value
				}
				pval
			} ) )
		}
	}

	rm(IncNum_tb, allPSI_tb, avgPSI_tb, IncSumNum_tb, ExcSumNum_tb)
	save.image(out_img_f)

	##sort
	exon_defi_d=exon_defi_d[ order( -apply(exon_defi_d[ttest_pval_names],1,min,na.rm=T), apply(abs(exon_defi_d[delta_PSI_names]),1,max,na.rm=T), decreasing=T) ,]
	#add gene description
	if(! ("gene_desc" %in% names(exon_defi_d) ) ){
		tax_id <- tax_ids[species]
		gene_info_file=paste("../ReferenceDB/gene/01gene_alias2id/",spe_common_names[as.character(tax_id)],".",tax_id,".tbl",sep="")
		if(!file.exists(gene_info_file)){
			gene_info_file=paste("../ReferenceDB/ncbi/gene/geneinfo/All_Mammalia.gene_info",sep="")
		}
		exon_defi_d[c("gene_id","gene_desc")]=gsb2_gene_desc(exon_defi_d$gene_symbol, spe=species, gene_info_file=gene_info_file)[c("gene_id","gene_desc")]
	}

	#output
	print(paste("output raw data to",out_tb_f))
	write.table(exon_defi_d, file=out_tb_f, sep="\t", col.names=T, row.names=F, quote=F)

	##output all detected exons 
	out_detected_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.detected.Rnum",PSI_cal_rnumMin,".tbl", sep="")
	mkdir_if_not_exist(out_detected_f)
	print(paste("output detected block to",out_detected_f))
	rm_headers=c("genef_row","intron_revnum","ss3_seq","junc_id","ups_intron","dns_intron","GeneUnit_id","row","min_Padj","ss_id","cds_len","transcript_len",
		"region_id","endPos_id","startPos_id","unit_id","startPosType","endPosType","Gblock_id",suppNum_names)
	if_detected=rowSums(!is.na(exon_defi_d[avgPSI_names]))>0; table(if_detected)
	out_d=re_format_tb(exon_defi_d[if_detected,], delete_headers=rm_headers, front_headers=c("gene_symbol","gene_id"), 
		back_headers=c(PSI_names,avgPSI_names,"coordinates",delta_PSI_names,ttest_pval_names,fisher_pval_names,ReguType_names,"gene_desc"), 
		ch_header_from=c("introni"),ch_header_to=c("region_id") )
	write.table(out_d, file=out_detected_f, sep="\t", col.names=T, row.names=F, quote=F)
	rm(out_d)

}



if(any(what2do %in% c("all","calReguType"))  ){
	#calculate ReguType
	if(cal_PadjDeltaPSI){
		Pval_setting=paste(c("AdjP_",paste(PadjDeltaPSI_ps,".b",PadjDeltaPSI_bs,".e",PadjDeltaPSI_es,sep="")),  collapse="")
		Pval_setting2=paste(c("P-value adjusted delta PSI method: ",paste(PadjDeltaPSI_ps,"; b=",PadjDeltaPSI_bs,"; e=",PadjDeltaPSI_es,sep="")),  collapse="")
	}else{
		Pval_setting=paste(paste(use_p_type,P_cut,sep=""), collapse="_")
		Pval_setting2=paste(paste(use_p_type,P_cut,sep="<"), collapse=" ")
	}
	out_format_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".txt", sep="")
	out_allReguType_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".allReguType.txt", sep="")
	out_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".stats.tbl", sep="")
	print(paste("output regulation data to",out_format_f))
	mkdir_if_not_exist(out_format_f)
	rm_headers=c("genef_row","intron_revnum","ss3_seq","junc_id","ups_intron","dns_intron","GeneUnit_id","row","min_Padj","ss_id","cds_len","transcript_len",
		"region_id","endPos_id","startPos_id","unit_id","startPosType","endPosType","Gblock_id",suppNum_names)

	exon_regu_d=cal_PSI_regu_type(exon_defi_d, use_p_type, delta_PSI_names, P_cut, delta_PSI_cut, ReguType_names, sample_pairs, zeroDeltaPSI_Pval=zeroDeltaPSI_Pval)
	if(cal_PadjDeltaPSI){
		#PadjDeltaPSI_ps=c("Pttest","Pfisher"); PadjDeltaPSI_bs=c(-0.5,-0.2); PadjDeltaPSI_es=c(1e-6,1e-30); PadjDeltaPSI_cutoff=10
		exon_regu_d[,adjDeltaPSI_names]=unlist(myApply(1:length(adjDeltaPSI_names), function(i){
			apply(exon_regu_d[,c(delta_PSI_names[i],PadjDeltaPSI_pName_matrix[i,])],1,function(v){cal_Padj_change(v[1], v[-1], PadjDeltaPSI_bs, PadjDeltaPSI_es)})
		}))
		exon_regu_d=cal_PSI_regu_type( exon_regu_d, use_p_type, adjDeltaPSI_names, rep(1,length(use_p_type)), PadjDeltaPSI_cutoff, ReguType_names, sample_pairs, zeroDeltaPSI_Pval=rep(1,length(use_p_type)) )
	}

	#output exon_regu_d and stats
	exon_regu_d=exon_regu_d[!duplicated(exon_regu_d$coordinates), ]; nrow(exon_regu_d)
	if_regulated=rowSums(exon_regu_d[ReguType_names]=="UP" | exon_regu_d[ReguType_names]=="DN")>0; table(if_regulated)
	exon_regu_d=re_format_tb(exon_regu_d, delete_headers=rm_headers, front_headers=c("gene_symbol","coordinates","gene_id"), 
		back_headers=c(log2R_names, adjpval_names, avgPSI_names,delta_PSI_names,adjDeltaPSI_names,ttest_pval_names,fisher_pval_names,ReguType_names,"gene_desc", IncUpsNum_names,IncDnsNum_names,ExcNum_names,PSI_names), 
		ch_header_from=c("introni"),ch_header_to=c("region_id") )
	names(exon_regu_d); nrow(exon_regu_d)
	write.table(exon_regu_d[if_regulated,], file=out_format_f, sep="\t", col.names=T, row.names=F, quote=F)
	
	#write Excel file:
	exon_regu_d2<-exon_regu_d[if_regulated,]
	if("openxlsx" %in% installed.packages()){
		library("openxlsx")
		out_xlsx_f= paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".xlsx", sep="")
		print(paste("write",out_xlsx_f))
		wb <- createWorkbook()
		sheet_name=ifelse(cal_PadjDeltaPSI, paste0("PadjDeltaPSI",delta_PSI_cut), paste0("dPSI",delta_PSI_cut,"_",Pval_setting))
		addWorksheet(wb, sheetName = sheet_name ) #
		freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
		writeDataTable(wb, sheet = 1, x = exon_regu_d2, colNames = T, rowNames = F,  firstColumn=T)
		LeftBorderStyle <- createStyle( border = "left")
		leftBorderHeaders=intersect(c( IncUpsNum_names[1],IncDnsNum_names[1],ExcNum_names[1], avgPSI_names[1], delta_PSI_names[1],adjDeltaPSI_names[1], fisher_pval_names[1], ttest_pval_names[1], adjpval_names[1], log2R_names[1],ReguType_names[1],"gene_desc") ,names(exon_regu_d2))
		leftBorderCols=match(leftBorderHeaders , names(exon_regu_d2) ) 
		for(leftBorderCol in leftBorderCols){
			addStyle(wb, 1, style = LeftBorderStyle, rows = 1:(nrow(exon_regu_d2)+1), cols=leftBorderCol, gridExpand = T)
		}
		headerStyle <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", border = "Bottom", textRotation=90)
		addStyle(wb, 1, style = headerStyle, rows = 1, cols=1:ncol(exon_regu_d2), gridExpand = F, stack=T)
		setColWidths(wb, sheet = 1, cols =1:ncol(exon_regu_d2), widths = 5)
		shortCol_headers=intersect(c("strand", "exon_i_inFrame","exon_i_stopCodon","ups_junc_num","dns_junc_num","skip_junc_num", ReguType_names), names(exon_regu_d2))
		setColWidths(wb, sheet = 1, cols = match(shortCol_headers, names(exon_regu_d2) ), widths = 3)
		wideCol_headers=intersect(c("gene_desc","startSS_seq","endSS_seq"), names(exon_regu_d2))
		setColWidths(wb, sheet = 1, cols = match(wideCol_headers, names(exon_regu_d2) ), widths = 30)
		setCol_headers=intersect(c("gene_symbol"), names(exon_regu_d2))
		setColWidths(wb, sheet = 1, cols = match(setCol_headers, names(exon_regu_d2) ), widths = 8)
		upStyle <- createStyle( bgFill = "red")
		dnStyle <- createStyle( bgFill = "royalblue") 
		conditionalFormatting(wb, 1, cols=match(ReguType_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "contains", rule="UP", style = upStyle)
		conditionalFormatting(wb, 1, cols=match(ReguType_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "contains", rule="DN", style = dnStyle)
		conditionalFormatting(wb, 1, cols=match(delta_PSI_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "colourScale", rule=c(-50,0,50), style = c("royalblue","white","red"))
		if(all(adjDeltaPSI_names %in% names(exon_regu_d2))){
			conditionalFormatting(wb, 1, cols=match(adjDeltaPSI_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "colourScale", rule=c(-30,0,30), style = c("royalblue","white","red"))
		}
		
		sel_L2R_names=intersect(log2R_names, names(exon_regu_d2) )
		if(length(sel_L2R_names)>0 ){
			conditionalFormatting(wb, 1, cols=match(sel_L2R_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "colourScale", rule=c(-3,0,3), style = c("royalblue","white","red"))
		}

		for(P_val_names in list(fisher_pval_names, ttest_pval_names, adjpval_names)){
			P_val_names=intersect(P_val_names, names(exon_regu_d2) )
			if(length(P_val_names)>0){
				conditionalFormatting(wb, 1, cols=match(P_val_names,names(exon_regu_d2)), rows=1+1:(nrow(exon_regu_d2)), type = "colourScale", rule=c(0,0.1), style = c("orange","white"))
			}
		}
		
		saveWorkbook(wb, out_xlsx_f, overwrite = TRUE)
		#write.xlsx(exon_regu_d2, file = out_xlsx_f, asTable = TRUE)
	}

	#output a table with rows as each comparison (sample pair), columns as numbers of UP or DN and sorted gene list AS events of UP or DN
	event_id_headers=c("gene_symbol","region_ano","coordinates","strand")
	tb=as.matrix(apply(exon_regu_d2[ReguType_names],2,function(v){table( factor(v,levels=ReguTypes) )})[ReguTypes,], nrow=length(ReguTypes))
	topGeneList_tb=data.frame(project=rep(study_name,length(sample_pairs)), sample_pair=sample_pairs)
	topGeneList_tb[paste0("Num_",ReguTypes[1:2])]=t(tb[1:2,])
	tmp_tb=t(matrix(unlist(myApply(1:length(sample_pairs), function(i){
		Pval_headers=paste(use_p_type,"_", sample_pairs[i],sep="")
		if(cal_PadjDeltaPSI){
			data_headers=c(adjDeltaPSI_names[i], delta_PSI_names[i]) #the first value of data_headers was used for sorting
		}else{
			data_headers=c(delta_PSI_names[i])
		}
		if_up=  exon_regu_d2[,ReguType_names[i]] %in% ReguTypes[1]
		if_dn=  exon_regu_d2[,ReguType_names[i]] %in% ReguTypes[2]
		return_v=c(0,"","","","")
		if(any(if_up | if_dn)){
			return_v[1]=max(abs(exon_regu_d2[if_up | if_dn, data_headers[1]]),  na.rm=T)
		}
		report_list=list(
			up=list(if_inlist=if_up, sort_decreasing=T, return_field=c(2,4) ), 
			down=list(if_inlist=if_dn, sort_decreasing=F, return_field=c(3,5) )
		)
		for(report_cat in names(report_list)){
			report_list1=report_list[[report_cat]]
			if_inlist=report_list1$if_inlist
			if(any(if_inlist)){
				Topevent_list=exon_regu_d2[if_inlist, c(event_id_headers, data_headers, Pval_headers)]
				Topevent_list=Topevent_list[order(Topevent_list[,data_headers[1]],decreasing=report_list1$sort_decreasing),]
				if(nrow(Topevent_list)>max_geneList_num){Topevent_list=Topevent_list[1:max_geneList_num,]}
				Topevent_list[,data_headers]=round(Topevent_list[,data_headers],1)
				Topevent_list[,Pval_headers]=format(Topevent_list[,Pval_headers],scientific=T, digits=2)
				return_v[report_list1$return_field[1]]=paste(apply(Topevent_list[,event_id_headers],1,paste,collapse=":"), collapse=" ")
				return_v[report_list1$return_field[2]]=paste(paste(1:nrow(Topevent_list), apply(Topevent_list[,event_id_headers],1,paste,collapse=":"), apply(Topevent_list[,c(data_headers,Pval_headers)],1,paste,collapse=":"), sep="|"), collapse=" ")
			}
		}
		return(return_v)
	})), nrow=5))
	tmp_tb[,4]=gsub(": +",":", gsub("\\| +","|", tmp_tb[,4]))
	tmp_tb[,5]=gsub(": +",":", gsub("\\| +","|", tmp_tb[,5]))
	topGeneList_tb[,c("max_absDeltaPSI",paste0("List_",ReguTypes[1:2]), paste0("ListDetail_",ReguTypes[1:2]) )] = tmp_tb
	topGeneList_tb$max_absDeltaPSI=as.numeric(unlist(topGeneList_tb$max_absDeltaPSI))
	out_topGeneList_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".topGeneList.csv", sep="")
	write.csv(topGeneList_tb, file=out_topGeneList_f,  row.names=F)

	
	#write a table with all exons and regulation type (used for cis elements analysis)
	mkdir_if_not_exist(out_allReguType_f)
	exon_regu_d$if_gt_ag_ss_exon=as.numeric( exon_regu_d$region_ano=="exon" & substring(exon_regu_d$startSS_seq,49,50)=="ag" & substring(exon_regu_d$endSS_seq,11,12)=="gt" )
	exon_regu_d$if_GAss5_exon=as.numeric( exon_regu_d$region_ano=="exon" &  substring(exon_regu_d$endSS_seq,9,10)=="GA" )
	print(table(exon_regu_d$if_gt_ag_ss_exon, exon_regu_d$region_ano))
	out_format_d=re_format_tb(exon_regu_d, 
		delete_headers=c("GeneUnit_id","row","gene_desc","min_Padj","ss_id","cds_len","transcript_len","oordinates",PSI_names,delta_PSI_names,adjDeltaPSI_names,adjpval_names,ttest_pval_names,fisher_pval_names,"gene_desc"), 
		del_pattern="num_|^P_|Log2R_|Padj_|\\.P_|Num_|PSI_|\\..*_pos|\\.reg_len|\\.region_ano|_junc_num|_ConsSco|SS_seq",
		back_headers=c(ReguType_names) )
	out_format_d$Gblock_id2=paste(out_format_d$contig, out_format_d$strand, format(out_format_d$start_pos,scientific=F,trim=T), format(out_format_d$end_pos,scientific=F,trim=T),  sep=":")
	write.table(out_format_d, file=out_allReguType_f, col.names=T, row.names=F, sep="\t", quote=F)
	
	###stats
	#regulation number for individual samples
	library("fmsb")
	tb=apply(exon_regu_d[ReguType_names],2,function(v){table(factor(v, levels=ReguTypes) )[ReguTypes] })
	write(paste("criteria: |delta_PSI|>",delta_PSI_cut,"; ",Pval_setting2,  sep=""), file=out_stats_f)
	write.table(cbind(ReguTypes,tb), file=out_stats_f, sep="\t", col.names=T, row.names=F, quote=F, append=T)

	region_anos=unique(exon_regu_d$region_ano) ; region_anos=region_anos[!is.na(region_anos)]
	for(region_ano in region_anos){
		print(region_ano)
		ifsel=!is.na(exon_regu_d$region_ano) & exon_regu_d$region_ano==region_ano
		tb=as.matrix( apply(as.matrix(exon_regu_d[ifsel, ReguType_names]), 2, function(v){ table(factor(v,levels=ReguTypes))} ) [ReguTypes,] , nrow=length(ReguTypes) )
		colnames(tb)=ReguType_names
		write(paste("\nRegion annotation=",region_ano),  file=out_stats_f, append=T)
		write.table(cbind(ReguTypes,tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
	}
	
	if(length(ReguType_names)>1){
		type_headers1=ReguType_names  ;  type_headers2=ReguType_names
		for(region_ano in region_anos){
			print(region_ano)
			ifsel=!is.na(exon_regu_d$region_ano) & exon_regu_d$region_ano==region_ano
			compare_types_among_vars(out_stats_f, exon_regu_d[ifsel,], type_headers1=type_headers1, Types=ReguTypes, type_headers2=type_headers2, 
				prefix=paste("\ncross table: region_ano=",region_ano) )
		}
	}

}

	
if(any(what2do %in% c("plot_delta_vs_basalPSI"))  ){
	##draw scatter plot of delta-PSI vs. basal-PSI
	basal_psi_cuts=c(0,20,90,100)
	out_image_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,".deltaVSbasalPSI.scatter.pdf", sep="")
	out_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,".deltaVSbasalPSI.stats.tbl", sep="")
	write(paste("\nstudy basal PSI and PSI change:"), file=out_stats_f, append=F)
	pdf(out_image_f, width=15,height=10); layout(matrix(1:6,2))
	for(i in 1:length(ReguType_names)){
		basal_psi_name=paste("PSI_",sample_compare_matrix[2,i],sep="")
		basal_psi_grp=cut(exon_regu_d[,basal_psi_name], breaks=basal_psi_cuts, include.lowest=T, ordered_result=T)
		tb=table(basal_psi_grp, factor(exon_regu_d[,ReguType_names[i]],levels=ReguTypes))[, ReguTypes]
		write(paste("basal PSI (rows)",basal_psi_name,"vs. PSI change", ReguType_names[i]), file=out_stats_f, append=T)
		write.table(cbind(rownames(tb),tb), file=out_stats_f, sep="\t", col.names=T, row.names=F, quote=F, append=T)
		x=exon_regu_d[,basal_psi_name]
		y=exon_regu_d[,delta_PSI_names[i]]
		plot(x, y, type="n", xlim=c(0,100), ylim=c(-100,100), 
			xlab=paste(basal_psi_name,"(%)"), ylab=paste(delta_PSI_names[i],"(%)"), main=delta_PSI_names[i])
		abline(a=100,b=-1, col="gray")
		abline(a=0,b=-1, col="gray")
		points(x[exon_regu_d[,ReguType_names[i]]=="NC"], y[exon_regu_d[,ReguType_names[i]]=="NC"], col="gray")
		text(x[exon_regu_d[,ReguType_names[i]]=="UP"], y[exon_regu_d[,ReguType_names[i]]=="UP"], col="red", labels=exon_regu_d$gene_symbol[exon_regu_d[,ReguType_names[i]]=="UP"],cex=1)
		text(x[exon_regu_d[,ReguType_names[i]]=="DN"], y[exon_regu_d[,ReguType_names[i]]=="DN"], col="blue", labels=exon_regu_d$gene_symbol[exon_regu_d[,ReguType_names[i]]=="DN"],cex=1)
	}
	dev.off()
}

	
if(any(what2do %in% c("ana_m4m3_diNT_freq","ana_deltaPSI_vs_m4m3_diNT","ana_deltaPSItype_vs_m4m3_diNT") ) ){
	##study the -4 and -3 position NT frequency for GA-type exons:
	exon_defi_d$avg_basal_PSI=rowMeans(exon_defi_d[basal_PSI_names],na.rm=T)
	basal_psi_cuts=c(0,10,95,100)
	exon_defi_d$avg_basal_PSI_bin=cut(exon_defi_d$avg_basal_PSI, breaks=basal_psi_cuts, include.lowest=T); table(exon_defi_d$avg_basal_PSI_bin)
	basal_PSI_bins=levels(exon_defi_d$avg_basal_PSI_bin)

	exon_defi_d$ss5m4=substr(exon_defi_d$endSS_seq,7,7)
	exon_defi_d$ss5m3=substr(exon_defi_d$endSS_seq,8,8)
	exon_defi_d$ss5m4m3=substr(exon_defi_d$endSS_seq,7,8)
	exon_defi_d$ss5m2m1=substr(exon_defi_d$endSS_seq,9,10)
}

##study the diNT frequency in GAgt exons with different basal PSI level
if(any(what2do %in% c("ana_m4m3_diNT_freq") ) ){
	ifsel=exon_defi_d$region_ano=="exon" & substr(exon_defi_d$endSS_seq,11,12)=="gt" & exon_defi_d$ss5m2m1=="GA" ; table(ifsel)
	exon_defi_d2=exon_defi_d[ifsel,]
	sum_tb=NULL
	for(basal_PSI_bin in basal_PSI_bins){
		ifsel=exon_defi_d2$avg_basal_PSI_bin==basal_PSI_bin
		exon_d2=exon_defi_d2[ifsel,]
		tb1=table(factor(exon_d2[ , "ss5m4"],levels=NTs), factor(exon_d2[ , "ss5m3"],levels=NTs) ) [NTs,NTs] #each value is a di-NT frequency
		m4Nums1=rowSums(tb1); m3Nums1=colSums(tb1)
		
		sum_tb=rbind(sum_tb, c(m4Nums1, m4Nums1/sum(m4Nums1), m3Nums1, m3Nums1/sum(m3Nums1), 
			tb1, sum(tb1), tb1/sum(tb1) ) )
	}
	colnames(sum_tb)=c(paste("m4",NTs,sep="."),paste("m4",NTs,"perc",sep="."), paste("m3",NTs,sep="."),paste("m3",NTs,"perc",sep="."),
		 diNTs,"total",  paste(diNTs,"perc",sep=".") )
	rownames(sum_tb)=basal_PSI_bins

	out_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.basal_PSI_vs_m4m3NT",".stats.tbl", sep="")
	write.table(cbind(rownames(sum_tb),sum_tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F)
}

if(any(what2do %in% c("ana_deltaPSItype_vs_m4m3_diNT") ) ){
	Pval_setting=paste(paste(use_p_type,P_cut,sep=""), collapse="_");Pval_setting
	ifsel=exon_defi_d$region_ano=="exon" & substr(exon_defi_d$endSS_seq,11,12)=="gt" & substr(exon_defi_d$startSS_seq, 49, 50)=="ag" & exon_defi_d$reg_len<200; table(ifsel)
	exon_defi_d2=exon_defi_d[ifsel,]
	exon_defi_d2$region_annotation="Either"
	exon_defi_d2$region_annotation[exon_defi_d2$endSS_supp %in% c("Refseq","KnG_Ens") & exon_defi_d2$startSS_supp %in% c("Refseq","KnG_Ens")]="annotated"
	exon_defi_d2$region_annotation[!(exon_defi_d2$endSS_supp %in% c("Refseq","KnG_Ens")) & !(exon_defi_d2$startSS_supp %in% c("Refseq","KnG_Ens"))]="None"
	table(exon_defi_d2$region_annotation)
	exon_defi_d2=cal_PSI_regu_type(exon_defi_d2, use_p_type, delta_PSI_names, P_cut, delta_PSI_cut, ReguType_names, sample_pairs)

	selReguType_names=ReguType_names
	out_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.m4m3NT_vs_PSIreguType.delta_PSI",delta_PSI_cut,"_",Pval_setting,".stats.tbl", sep="")
	#study "5\'ss motif vs. AS regulation"
	write(paste("#5\'ss motif vs. AS regulation (delta_PSI",delta_PSI_cut," ",Pval_setting,")",sep=""), file=out_stats_f, append=F)
	sum_tb=NULL
	for(ss5_ano in c("None","annotated")){
		ifsel1=!is.na(exon_defi_d2$region_annotation) & exon_defi_d2$region_annotation==ss5_ano
		for(ss5m2m1 in c("GA","AG") ){
			ifsel2=ifsel1 & exon_defi_d2$ss5m2m1==ss5m2m1
			for(ReguType_name in selReguType_names){
				for(ReguType in c("UP","NC")){
					ifsel3=ifsel2 & exon_defi_d2[,ReguType_name] %in% ReguType & !is.na(exon_defi_d2[,ReguType_name])
					exon_defi_d3=exon_defi_d2[ifsel3,]
					nums=table(factor(exon_defi_d3$ss5m4m3,levels=diNTs) )[diNTs]
					total_num=sum(nums)
					percs=nums/total_num
					sum_tb=rbind(sum_tb, c(ss5_ano, ss5m2m1, ReguType_name, ReguType, total_num,nums,percs) )
				}
			}
		}
	}
	colnames(sum_tb)=c("ss5_ano","ss5m2m1","ReguType_name","ReguType", "total_num", 
		diNTs, paste("Perc_",diNTs,sep="") )
	write.table(sum_tb, file=out_stats_f, append=T, col.names=T, row.names=F, sep="\t", quote=F)
}

if(any(what2do %in% c("ana_deltaPSI_vs_m4m3_diNT") ) ){
	##study change of splicing in different diNT and basel PSI groups
	ifsel=exon_defi_d$ss5m2m1=="GA"
	exon_defi_d2=exon_defi_d[ifsel,]
	table(exon_defi_d2$avg_basal_PSI_bin)
	out_image_f=paste(out_root,"/",study_exon_type,"/exons.deltaPSI_vs_m4m3NT.Rnum",PSI_cal_rnumMin,".boxplot.pdf", sep="")
	print(paste("ana_deltaPSI_vs_m4m3_diNT",out_image_f))
	NTcolors=rainbow(length(diNTs))

	selbasal_PSI_bins=basal_PSI_bins[c(1,3)]
	ylim_l=list(c(-5,60), c(-10,10))
	names(ylim_l)=selbasal_PSI_bins
	min_delta_PSI=2
	pdf(out_image_f, width=20, height=8)
	for(basal_PSI_bin in selbasal_PSI_bins){
		layout(matrix(1:4,2))
		if_regulatable=rowSums(abs(exon_defi_d2[delta_PSI_names])>min_delta_PSI)>0
		exon_defi_d3=exon_defi_d2[exon_defi_d2$avg_basal_PSI_bin==basal_PSI_bin & if_regulatable,]
		for(i in 1:length(delta_PSI_names)){
			# draw_grped_boxplot(exon_defi_d3,"ss5m4m3", diNTs, delta_PSI_names[i], title=paste(basal_PSI_bin,delta_PSI_names[i],"\n"), 
			# 	 xlabel="-4 -3 position of ss5", ylabel=delta_PSI_names[i], hline_at=0 , colorBy_vars=F, spaceBwGrp=0, boxwex=1.5,
			# 	 ifoutline=F)
			exon_defi_d4=exon_defi_d3[!is.na(exon_defi_d3[,delta_PSI_names[i] ]) , ]
			NTnums=table( factor(exon_defi_d4$ss5m4m3,levels=diNTs) )[diNTs]
			boxplot(exon_defi_d4[,delta_PSI_names[i]] ~ factor(exon_defi_d4$ss5m4m3,levels=diNTs),  outpch = NA, xlab="-4 -3 position of ss5", ylab=delta_PSI_names[i], 
				main=paste("baselPSI Bin",basal_PSI_bin,delta_PSI_names[i],"\n", paste(paste(diNTs,NTnums,sep="="), collapse=" ") ), ylim=ylim_l[[basal_PSI_bin]] )
			stripchart(exon_defi_d4[,delta_PSI_names[i]] ~ factor(exon_defi_d4$ss5m4m3,levels=diNTs), vertical = TRUE, method = "jitter", jitter=0.3, pch = 1, cex=1, 
				col=adjustcolor(NTcolors, alpha.f = 0.5),  add = TRUE) 
			abline(h=0, col="gray")

		}
	}
	dev.off()
}

	
if(any( what2do %in% c("study_ssScore_vs_basalPSI") ) ){ 
	##draw scatter plot of delta-PSI vs. basal-PSI
	basal_psi_cuts=c(0,20,95,100); basal_psi_cut_bins=c("below20","mid","above95"); 
	study_ssSco_names=c("ups_ss5_sco","ups_ss3_sco","ss5_sco","ss3_sco") # "ups_ss3_sco","ss5_sco" are the two ss for the alt-exon
	out_image_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.ssScoVSbasalPSI.scatter.pdf", sep="")
	out_data_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.ssScoVSbasalPSI.data.tbl", sep="")
	basal_PSI_names
	PSI_cal_rnumMin2=10
	#re-calculate PSI using new PSI_cal_rnumMin2
	IncSumNum_tb=exon_defi_d[,IncUpsNum_names]+exon_defi_d[,IncDnsNum_names]
	IncNum_tb=IncSumNum_tb
	IncNum_tb[exon_defi_d$region_ano=="exon",]=round(IncNum_tb[exon_defi_d$region_ano=="exon",]/2,1)
	allPSI_tb=round(IncNum_tb/(IncNum_tb+exon_defi_d[, ExcNum_names])*100,1)
	allPSI_tb[(IncSumNum_tb+exon_defi_d[, ExcNum_names])<PSI_cal_rnumMin2]=NA
	ExcSumNum_tb=exon_defi_d[, ExcNum_names]
	colnames(allPSI_tb)=all_samples
	colnames(IncSumNum_tb)=all_samples
	colnames(ExcSumNum_tb)=all_samples
	#calculate average PSI, 
	avgPSI_tb=sapply(avg_samples,function(avg_sample1){
		if(avg_sample1 %in% names(sample_repl_l)){
			rowMeans(allPSI_tb[, sample_repl_l[[avg_sample1]]], na.rm=T)
		}else{
			allPSI_tb[,avg_sample1]
		}
	})
	exon_defi_d[PSI_names]=round(allPSI_tb,0)
	exon_defi_d[avgPSI_names]=round(avgPSI_tb,1)
	exon_defi_d2=exon_defi_d[exon_defi_d$region_ano=="exon",]
	exon_defi_d2$avg_basal_PSI=rowMeans(exon_defi_d2[,basal_PSI_names],na.rm=T)
	exon_defi_d2$avg_basal_PSI_bin=cut(exon_defi_d2$avg_basal_PSI, breaks=basal_psi_cuts, include.lowest = T, labels=basal_psi_cut_bins)
	print(table(exon_defi_d2$avg_basal_PSI_bin))
	write.table(exon_defi_d2,file=out_data_f, col.names=T, row.names=F, sep="\t", quote=F)

	pdf(out_image_f, width=15,height=10); layout(matrix(1:6,2,byrow=T))
	draw_ecdf_density(exon_defi_d2, study_ssSco_names, "avg_basal_PSI_bin", rev(basal_psi_cut_bins), mycols=c("red","green","blue"), test_sample_pairs = "first_oths", 
		test_method = "wilcox.test", verticals=T, pch=NA, lwd=2, method = "ecdf")
	layout(1)
	draw_grped_boxplot(exon_defi_d2, "avg_basal_PSI_bin", rev(basal_psi_cut_bins), study_ssSco_names, mycolor=c("red","green","blue"), test_sample_pairs = "first_oths", 
		test_method = "wilcox.test", colorBy_vars=F, join_by_grps=T, ymin=0)
	dev.off()
}
