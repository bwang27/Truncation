##further analyze the output of 02b.ana_chimeric_reads.pl, find the best threshold to define a chimeric-reads




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
out_dir=chimeric_in_dir
if(!is.na(args_v["chimeric_in_dir"])){ chimeric_in_dir= as.character(args_v["chimeric_in_dir"]) }
if(!is.na(args_v["out_dir"])){ out_dir= as.character(args_v["out_dir"]) }

chimeric_simple_types=c("Boundary","Inverted","Deletion","Insertion","Reversion","unkn")
chimeric_types_colors=c("black","red","blue","pink","orange","gray")
names(chimeric_types_colors)=chimeric_simple_types

##1, read and load all .chimericRead.category.txt files
merged_d=NULL
for (sample in samples){
	inf=paste0(chimeric_in_dir, sample, ".chimericRead.category.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F)
	d1=data.frame(sample=rep(sample,nrow(d)))
	d1=cbind(d1,d)
	merged_d=rbind(merged_d, d1)
}
table(merged_d$sample)
table(merged_d$sample,merged_d$chimeric_type)

##2, plot the distribution of min_mapq min_refCovLen min_ali_sco max_percMM twoBlockDist for different chimeric_type
merged_d$chimeric_type_simple=sub("\\+|\\-","",merged_d$chimeric_type); table(merged_d$chimeric_type_simple)
table(merged_d$chimeric_Ref_type,merged_d$chimeric_type_simple)
unique(merged_d$chimeric_type_simple)

merged_d$twoBlockDist[merged_d$chimeric_type_simple =="unkn"]=NA
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
		plot_chimeric_simple_types=intersect(chimeric_simple_types, merged_d$chimeric_type_simple[!is.na(merged_d[,plot_var])]) 
		draw_ecdf_density(merged_d, plot_var, "chimeric_type_simple", plot_chimeric_simple_types, desc=plot_sampleName, method="density", test_sample_pairs="first_oths", lwd=3, scale_Start_quantile=0.01,
			mycols=chimeric_types_colors[plot_chimeric_simple_types], xmax=xmax, iflegend=plot_var==plot_vars[length(plot_vars)])
	}
dev.off()


