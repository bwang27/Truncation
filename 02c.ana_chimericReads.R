##further analyze the output of 02b.ana_chimeric_reads.pl, find the best threshold to define a chimeric-reads
##2022/11/01 added function to annotate potential fusion events generating rcAAV (in MergedReadCatCount_d), need to define geneFeature_bed_f



args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../../../Pipeline/SharedCodes/Rfunc.inc.R")

study_name="2022-04-19_AAV_SEQ"
samples=unlist(strsplit("ID23_scAAV-SCA3-4E10x2_AAV SC-AAV9-SCA3-E10x2-CAGDel7-ID43-01-02"," "))

if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }
if(!is.na(args_v["samples"])){ samples= unlist(strsplit(as.character(args_v["samples"])," "))}

chimeric_in_dir=paste0("/home/Research_Archive/ProcessedData/Nanopore_seq/",study_name, "/02b_ChimericRead/")
if(!is.na(args_v["chimeric_in_dir"])){ chimeric_in_dir= as.character(args_v["chimeric_in_dir"]) }
out_dir=chimeric_in_dir
if(!is.na(args_v["out_dir"])){ out_dir= as.character(args_v["out_dir"]) }

chimeric_simple_types=c("Boundary","Inversion","Deletion","Insertion","Reversion","fusion")
chimeric_types_colors=c("black","red","blue","pink","orange","gray")
names(chimeric_types_colors)=chimeric_simple_types

min_readCount=ifelse(!is.na(args_v["min_readCount"]), as.numeric(args_v["min_readCount"]), 1) #minimal read count for output
min_readCount2max_ratio=ifelse(!is.na(args_v["min_readCount2max_ratio"]), as.numeric(args_v["min_readCount2max_ratio"]), 0) #minimal percent of read counts relative the the maximal read count for the same sample
readCounts_names=paste0("Num_",samples)

#geneFeature_bed_f=paste0("/home/Research_Archive/ProcessedData/Nanopore_seq/target_index/", "pDG-KanR.NA.pAAV-hDDC-KanR-shift/pDG-KanR.NA.pAAV-hDDC-KanR-shift.bed")
index_folder="/home/Research_Archive/ProcessedData/Nanopore_seq/target_index"
if(!is.na(args_v["index_folder"])){ index_folder= as.character(args_v["index_folder"]) }
if(!is.na(args_v["all_indexes"])){ all_indexes= unlist(strsplit(as.character(args_v["all_indexes"])," "))}
indexes2sample_l=tapply(samples,all_indexes,function(v){v})

ITR_feature_patt="ITR"
Rep_feature_patt="Rep"
Cap_feature_patt="Cap"
rcAAV_fusion_max_dist=250 #maximal distance (bp) allowed to call a fusion between an ITR to Rep/Cap a rcAAV fusion event
missing_Seq_max_size=100; #maximal size of nucleotides allowed to be deleted because of a fusion event (eg, part of the ITR, Rep or Cap sequences)

nCores=ifelse(!is.na(args_v["nCores"]), as.numeric(args_v["nCores"]), 6) #minimal read count for output
define_parallel_fun(nCores=nCores)


out_R_img_f=paste0(out_dir,  "02c.ana_chimericReads.R.img")
out_MergedReadCatCount_f=paste0(out_dir,  "All.ChimericReadCatCount.txt")
save.image(out_R_img_f)


##1, read and load all .chimericRead.detail.txt files
mergedReadDetail_d=NULL
for (sample in samples){
	inf=paste0(chimeric_in_dir, sample, ".chimericRead.detail.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F)
	d1=data.frame(sample=rep(sample,nrow(d)))
	d1=cbind(d1,d)
	mergedReadDetail_d=rbind(mergedReadDetail_d, d1)
}
table(mergedReadDetail_d$sample)
table(mergedReadDetail_d$sample,mergedReadDetail_d$chimeric_type)

##2, plot the distribution of min_mapq min_refCovLen min_ali_sco max_percMM twoBlockDist for different chimeric_type
mergedReadDetail_d$chimeric_type_simple=sub("\\+|\\-","",mergedReadDetail_d$chimeric_type); table(mergedReadDetail_d$chimeric_type_simple)
table(mergedReadDetail_d$chimeric_Ref_type,mergedReadDetail_d$chimeric_type_simple)
unique(mergedReadDetail_d$chimeric_type_simple)

mergedReadDetail_d$twoBlockDist[mergedReadDetail_d$chimeric_type_simple =="fusion"]=NA
##make density plot
plot_vars=unlist(strsplit("min_mapq min_refCovLen min_ali_sco twoBlockDist max_percMM"," ")); plot_vars
#xmaxs=c(2000); names(xmaxs)="twoBlockDist"
xmaxs=NULL

#par(las=2, mar=c(9,12,2,0.5)+0.1, xaxt="n")
plot_sampleName=ifelse(length(samples)>1,"AllSamples",samples)
out_chimericReadQualityPlot_f=paste0(out_dir, "density/", plot_sampleName,".chimericReads.parameters.densityPlot.pdf")
mkdir_if_not_exist(out_chimericReadQualityPlot_f)
pdf(out_chimericReadQualityPlot_f, width=20, height=10)
layout(matrix(1:6,2, byrow=T))
	for(plot_var in plot_vars){
		xmax=NULL; if(plot_var %in% names(xmaxs)){xmax=xmaxs[plot_var]}; 
		plot_chimeric_simple_types=intersect(chimeric_simple_types, mergedReadDetail_d$chimeric_type_simple[!is.na(mergedReadDetail_d[,plot_var])]) 
		draw_ecdf_density(mergedReadDetail_d, plot_var, "chimeric_type_simple", plot_chimeric_simple_types, desc=plot_sampleName, method="density", test_sample_pairs="first_oths", lwd=3, scale_Start_quantile=0.01,
			mycols=chimeric_types_colors[plot_chimeric_simple_types], xmax=xmax, iflegend=plot_var==plot_vars[length(plot_vars)])
	}
dev.off()

rm(mergedReadDetail_d,d1)

##3, combine chimeric reads results (eg. fusion)

ReadCatCount_dlist=list()
readMin_v=NULL
for (sample in samples){
	inf=paste0(chimeric_in_dir, sample, ".chimericReadCount.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F)
	readMin=max(c(min_readCount, max(d$count, na.rm=T)*min_readCount2max_ratio ))
	readMin_v=c(readMin_v, readMin)
	names(d)=sub("^count$",paste0("Num_",sample),names(d))
	ReadCatCount_dlist[[sample]]=d
}

#merge count data
readMin_v
readMin4All=min(readMin_v); readMin4All
MergedReadCatCount_d=NULL
for (sample in samples){
	d=ReadCatCount_dlist[[sample]]
	if_sel=d[,paste0("Num_",sample)]>=readMin4All
	if(is.null(MergedReadCatCount_d)){
		MergedReadCatCount_d=d[if_sel,]
	}else{
		MergedReadCatCount_d=merge(MergedReadCatCount_d, d[if_sel,], all=T)
	}
	print(dim(MergedReadCatCount_d))
}
MergedReadCatCount_d=re_format_tb(MergedReadCatCount_d, back_headers=c(readCounts_names, "template_seq"))
MergedReadCatCount_d[readCounts_names][is.na(MergedReadCatCount_d[readCounts_names])]=0
MergedReadCatCount_d=MergedReadCatCount_d[order(apply(MergedReadCatCount_d[,readCounts_names,drop=F],1,max), decreasing=T), ]

judge_rcAAV_fusion<-function(ref,pos,str,if_aliBlock1){
	if(!if_aliBlock1){
		str=ifelse(str %in% "+","-","+")
	}
	for(i in 1:length(rcAAV_fusion_l)){
		#print(names(rcAAV_fusion_l)[i])
		if(ref == rcAAV_fusion_l[[i]][[1]]){
			if(str == rcAAV_fusion_l[[i]][[4]]){
				if(pos>=rcAAV_fusion_l[[i]][[2]] & pos<=rcAAV_fusion_l[[i]][[3]]){
					dist=(pos-rcAAV_fusion_l[[i]][[5]])*ifelse(str %in% "+",1,-1)
					return(c(names(rcAAV_fusion_l)[i],dist))
				}
			}
		}
	}
	return( c("",NA) )
}

for(i in 1:length(indexes2sample_l)){
	index_name=names(indexes2sample_l)[i]
	remove_readCount_names=setdiff(readCounts_names, paste0("Num_", indexes2sample_l[[i]]))
	geneFeature_bed_f=paste0(index_folder,"/",index_name,"/",index_name,".bed")
	geneFeature_bed_d=read.table(geneFeature_bed_f, header=F, sep="\t", quote="", stringsAsFactors=F)
	names(geneFeature_bed_d)=c("chromosome","start","end","name","score","strand")
	geneFeature_bed_d=geneFeature_bed_d[order(geneFeature_bed_d$chromosome,geneFeature_bed_d$start,-geneFeature_bed_d$end),]
	ITR_d=geneFeature_bed_d[grepl(ITR_feature_patt,geneFeature_bed_d$name) & !grepl("Ad5_ITR",geneFeature_bed_d$name),]
	Rep_d=geneFeature_bed_d[grepl(Rep_feature_patt,geneFeature_bed_d$name ),]
	Cap_d=geneFeature_bed_d[grepl(Cap_feature_patt,geneFeature_bed_d$name ),]
	Rep_d=Rep_d[order(Rep_d$end-Rep_d$start, decreasing=T),]
	Cap_d=Cap_d[order(Cap_d$end-Cap_d$start, decreasing=T),]
	if(nrow(ITR_d)>2){
		ITR_d=ITR_d[c(1, nrow(ITR_d)),]
	}

	if(nrow(ITR_d)==2 & nrow(Rep_d)>0 & nrow(Cap_d)>0){
		rcAAV_fusion_l=list(
			leftITR=list(ITR_d$chromosome[1], ITR_d$end[1]-missing_Seq_max_size,  ITR_d$end[1]+rcAAV_fusion_max_dist, "+", ITR_d$end[1])
			,rightITR=list(ITR_d$chromosome[2], ITR_d$start[2]-rcAAV_fusion_max_dist,  ITR_d$start[2]+missing_Seq_max_size, "-", ITR_d$start[2])
			,UpsRep=list(Rep_d$chromosome[1], ifelse(Rep_d$strand[1] %in% c("+","."), Rep_d$start[1]-rcAAV_fusion_max_dist, Rep_d$end[1]-missing_Seq_max_size),
				ifelse(Rep_d$strand[1] %in% c("+","."), Rep_d$start[1]+missing_Seq_max_size, Rep_d$end[1]+rcAAV_fusion_max_dist),
				ifelse(Rep_d$strand[1] %in% c("+","."), "-", "+"),
				ifelse(Rep_d$strand[1] %in% c("+","."), Rep_d$start[1], Rep_d$end[1]) )
			,DnsCap=list(Cap_d$chromosome[1], ifelse(Cap_d$strand[1] %in% c("+","."), Cap_d$end[1]-missing_Seq_max_size, Cap_d$start[1]-rcAAV_fusion_max_dist ),
				ifelse(Cap_d$strand[1] %in% c("+","."), Cap_d$end[1]+rcAAV_fusion_max_dist, Cap_d$start[1]+missing_Seq_max_size),
				ifelse(Cap_d$strand[1] %in% c("+","."), "+", "-"),
				ifelse(Cap_d$strand[1] %in% c("+","."), Cap_d$end[1], Cap_d$start[1]) )
			,leftITRrev=list(ITR_d$chromosome[1], ITR_d$start[1],  ITR_d$start[1]+missing_Seq_max_size, "-", ITR_d$start[1])
			,rightITRrev=list(ITR_d$chromosome[2], ITR_d$end[2]-missing_Seq_max_size,  ITR_d$end[2], "+", ITR_d$end[2])
		)
		rcAAV_plasmids=c(ITR_d$chromosome,Rep_d$chromosome)
		if_potential_rcAAV=MergedReadCatCount_d$chimeric_type=="fusion" & MergedReadCatCount_d$chimeric_Ref_type=="inter-GOI-others" & MergedReadCatCount_d$refName1 %in% rcAAV_plasmids & MergedReadCatCount_d$refName2 %in% rcAAV_plasmids ; table(if_potential_rcAAV)
		rc_AAV_d=MergedReadCatCount_d[if_potential_rcAAV, ]; nrow(rc_AAV_d)
		if(length(remove_readCount_names)>0){
			rc_AAV_d=rc_AAV_d[,setdiff(names(rc_AAV_d),remove_readCount_names)]
		}
		rc_AAV_d[c("ref1.rcAAV.ano","ref1.dist")]=matrix(unlist(myApply(1:nrow(rc_AAV_d),function(rowi){
			judge_rcAAV_fusion(ref=rc_AAV_d$refName1[rowi], pos=rc_AAV_d$refEnd1[rowi], str=rc_AAV_d$strand1[rowi], if_aliBlock1=TRUE)
		})), ncol=2, byrow=T)
		rc_AAV_d[c("ref2.rcAAV.ano","ref2.dist")]=matrix(unlist(myApply(1:nrow(rc_AAV_d),function(rowi){
			judge_rcAAV_fusion(ref=rc_AAV_d$refName2[rowi], pos=rc_AAV_d$refStart2[rowi], str=rc_AAV_d$strand2[rowi], if_aliBlock1=FALSE)
		})), ncol=2, byrow=T)
		table(rc_AAV_d$ref1.rcAAV.ano=="", rc_AAV_d$ref2.rcAAV.ano=="")
		out_rcAAV_fusion_f=paste0(out_dir, index_name,  ".rcAAV_fusion.txt")
		write.table(rc_AAV_d, file=out_rcAAV_fusion_f, col.names=T, row.names=F, sep="\t", quote=F, na="")
	}

}



write.table(MergedReadCatCount_d, file=out_MergedReadCatCount_f, col.names=T, row.names=F, sep="\t", quote=F)
rm(ReadCatCount_dlist,d) 

save.image(out_R_img_f)
