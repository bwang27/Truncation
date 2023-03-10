##analyze read mapping results (output of 02.count_AAV_read_fromBam.pl)



args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../../../Pipeline/SharedCodes/Rfunc.inc.R")


study_name="2019-11-11_AAV9-CMV-CDKL5-03"
samples=unlist(strsplit("Lam500ng_wDNase Vg50ulLam25ng Vg50ul H2O Lam500ng_woDNase"," "))

if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }
if(!is.na(args_v["samples"])){ samples= unlist(strsplit(as.character(args_v["samples"])," "))}

read_in_dir=paste0(study_name, "/02_Bam2ReadCounts/")
out_dir=paste0(study_name, "/03_ana_readMap/")
if(!is.na(args_v["read_in_dir"])){ read_in_dir= as.character(args_v["read_in_dir"]) }
if(!is.na(args_v["out_dir"])){ out_dir= as.character(args_v["out_dir"]) }

ref_sizes_f="target_index/hg19_CDKL5/hg19_CDKL5.fa.sizes"
if(!is.na(args_v["ref_sizes_f"])){ ref_sizes_f= as.character(args_v["ref_sizes_f"]) }

GOI_info_f="target_index/hg19_CDKL5/GOI_info.txt"
if(!is.na(args_v["GOI_info_f"])){ GOI_info_f= as.character(args_v["GOI_info_f"]) }


out_readMap_f=paste0(out_dir, "readMap.stats.txt")
out_data_img_f=paste0(out_dir, "readMap.stats.R.image")

mkdir_if_not_exist(out_readMap_f)

counts_names=paste("Num",samples, sep="_")
Length_names=paste("Len",samples, sep="_")
mismatch_names=paste("Mismatch",samples, sep="_")
PercCount_names=paste("PercCount",samples, sep="_")
PercLen_names=paste("PercLen",samples, sep="_")
PercMismatch_names=paste("PercMM",samples, sep="_")
AvgDensity_names=paste("AvgDensity",samples, sep="_")
AvgDensityVsGOI_Ratio_names=paste("AvgDensityVsGOI",samples, sep="_")

##2, combine all read counts files and combine, calculate percentage etc
for (sample in samples){
	inf=paste0(read_in_dir, sample, ".read.counts.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F)
	if(sample==samples[1]){
		comb_d=d
	}else{
		comb_d=merge(comb_d,d,all=T)
	}
}

##add reference/chromosome ref_size information
ref_sizes_d=read.table(ref_sizes_f, header=F, sep="\t", quote="", stringsAsFactors=F)
names(ref_sizes_d)=c("reference","ref_size")
ref_sizes_d$region="All"
ref_sizes_d=ref_sizes_d[c("reference","region","ref_size")]

GOI_info_d=read.table(GOI_info_f, header=T, sep="\t", quote="", stringsAsFactors=F)
GOI_info_d2=data.frame(reference=GOI_info_d$chromosome)
GOI_info_d2$region=GOI_info_d$name
GOI_info_d2$ref_size=GOI_info_d$to-GOI_info_d$from+1
ref_sizes_d=rbind(ref_sizes_d, GOI_info_d2, c("All","All", sum(ref_sizes_d$ref_size)) )
ref_sizes_d$ref_size=as.numeric(ref_sizes_d$ref_size)

comb_d=merge(comb_d,ref_sizes_d,all=T)
#comb_d$ref_size=ref_sizes_d$ref_size[match(paste(comb_d$reference, comb_d$region), paste(ref_sizes_d$reference, ref_sizes_d$region) )]
write.table(comb_d, file=out_readMap_f, sep="\t", row.names=F, col.names=T, quote=F)
save.image(out_data_img_f)

##calculate percentage etc.
val_headers=c(counts_names, Length_names, mismatch_names)
comb_d[, val_headers][ is.na(comb_d[, val_headers]) ]=0

if_human_chrom=comb_d$region %in% "All" & grepl("chr", comb_d$reference); comb_d$reference[if_human_chrom]
human_sum_d=comb_d[1,]
human_sum_d[,c("reference","region")]=c("Human.Chromosomes","All")
human_sum_d[,val_headers]=colSums(comb_d[if_human_chrom, val_headers])
comb_d=rbind(comb_d, human_sum_d)

if_sum_row=comb_d$reference=="All" & comb_d$region=="All" ; sum(if_sum_row)
sum_d=comb_d[if_sum_row, val_headers]

comb_d[PercCount_names]=t(t(comb_d[,counts_names])/unlist(sum_d[1,counts_names]))*100
comb_d[PercLen_names]=t(t(comb_d[,Length_names])/unlist(sum_d[1,Length_names]))*100
comb_d[PercMismatch_names]=comb_d[,mismatch_names] / comb_d[,Length_names]*100
comb_d[AvgDensity_names]=comb_d[,Length_names]/comb_d$ref_size
if_GOI=comb_d$region %in% "GOI_full"
if(sum(if_GOI)==1){
	comb_d[AvgDensityVsGOI_Ratio_names]=t(t(comb_d[,AvgDensity_names])/unlist(comb_d[if_GOI,AvgDensity_names]))
}

all_region_names=paste(comb_d$reference, comb_d$region,sep="."); all_region_names
all_human_regions=all_region_names[if_human_chrom]

##study the read coverage length distribution
all_readCov_l=list()
for (sample in samples){
	inf=paste0(read_in_dir, sample, ".read.Cov.distr.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F )
	region_names=names(d)[-1]
	d[,-1][is.na(d[,-1])]=0
	d$Human.Chromosomes.All=rowSums(d[intersect(names(d), all_human_regions)])
	d=d[,setdiff(names(d), all_human_regions)]
	groupped_regions=names(d)[-1]
	for( region1 in groupped_regions ){
		all_readCov_l[[sample]][[region1]]=rep(d[,1], d[,region1])
	}
}
#generate read coverage length stats table
readCov_stats_d=data.frame(reference.name=groupped_regions)
for(anaFun in c("mean","sd", "median", "max")){
	for(sample in samples){
		readCov_stats_d[, paste(sample,anaFun,sep=".")]=sapply(all_readCov_l[[sample]], anaFun)[groupped_regions]
	}
}

out_readCovStats_f=paste0(out_dir, "readCovLen.stats.txt")
write.table(readCov_stats_d, file=out_readCovStats_f, sep="\t", col.names=T, row.names=F, quote=F)

##make plot
out_readCovPlot_f=paste0(out_dir, "readCovLen.boxplot.pdf")
out_readCovPlot_pngf=paste0(out_dir, "readCovLen.boxplot.png")
data_range=c(min(unlist(all_readCov_l)), max(unlist(all_readCov_l)))

if(length(samples)>1){
	pdf(out_readCovPlot_f, width=18, height=8)
	layout(matrix(1:6,2, byrow=T))
}else{
	pdf(out_readCovPlot_f, width=10, height=5)
}
par(las=2, mar=c(9,12,2,0.5)+0.1, xaxt="n")
for(sample in samples){
	sorted_regions=groupped_regions[ order( readCov_stats_d[,paste(sample,"median",sep=".")], decreasing=F) ]
	GOI_region_name=grep(".GOI_full$", sorted_regions, value=T)
	sorted_regions=setdiff(sorted_regions, c("All.All",sub(".GOI_full$",".All",GOI_region_name) ) )
	median_lens=round(readCov_stats_d[match(sorted_regions, readCov_stats_d[,1]),paste(sample,"median",sep=".")])
	#stripchart(all_readCov_l[[sample]][sorted_regions], vertical = TRUE, method = "jitter", jitter=0.3, pch = ".", cex=1, col = adjustcolor("gray", alpha.f = 0.3),  add = TRUE) 
	#boxplot(all_readCov_l[[sample]][sorted_regions], horizontal=F, main=sample, log="y", ylim=data_range, ylab="Read Coverage Length")
	stripchart(all_readCov_l[[sample]][sorted_regions], vertical = TRUE, method = "jitter", jitter=0.3, pch = ".", cex=1, 
		col = adjustcolor("blue", alpha.f = 0.3), main=paste(c(sample,"\nMedian=",median_lens),collapse=" "), log="y", ylim=data_range, ylab="Read Coverage Length") 
	boxplot(all_readCov_l[[sample]][sorted_regions],   add = TRUE, outline = F)
	abline(h=median(all_readCov_l[[sample]][[GOI_region_name]]), lty=2, col="red")
	axis(1, at=seq(1, length(sorted_regions), by=1), labels = FALSE)
	text(seq(1, length(sorted_regions), by=1), 10^par("usr")[3] -1, labels = abbreviate(sorted_regions, minlength = 50, method="both.sides"), 
		srt = 25, adj=c(1,1),  xpd = TRUE)
}
dev.off()

convert_cmd=paste("convert  -density 150 ", out_readCovPlot_f, out_readCovPlot_pngf)
system(convert_cmd)

comb_d=re_format_tb( comb_d, front_headers=c("reference","region","ref_size",val_headers,PercCount_names,PercLen_names, PercMismatch_names ) )
write.table(comb_d, file=out_readMap_f, sep="\t", row.names=F, col.names=T, quote=F)
save.image(out_data_img_f)
