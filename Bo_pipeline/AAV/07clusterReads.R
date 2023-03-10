##cluster reads mapped to human genome and calculate enrichment and p-value


args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../../../Pipeline/SharedCodes/Rfunc.inc.R")

study_name="2020-06-Research7_CMC5_AAVs_Seq"
samples="SS-AAV9GTX-CMV-ReelinR3p6V5-03"
if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }
if(!is.na(args_v["samples"])){ samples= unlist(strsplit(as.character(args_v["samples"])," "))}

read_in_dir=paste0(study_name, "/02_Bam2ReadCounts/")
out_dir=paste0(study_name, "/07_readCluster/", ifelse(length(samples)==1, paste0(samples,"."), ""))

if(!is.na(args_v["read_in_dir"])){ read_in_dir= as.character(args_v["read_in_dir"]) }
if(!is.na(args_v["out_dir"])){ out_dir= as.character(args_v["out_dir"]) }

if(!is.na(args_v["ref_sizes_f"])){ ref_sizes_f= as.character(args_v["ref_sizes_f"]) }

out_readCluster_f=paste0(out_dir, "readMap.cluster.txt")
out_readCluster_bedf=paste0(out_dir, "readMap.cluster.bed")

out_data_img_f=paste0(out_dir, "readMap.cluster.R.image")

hostGenome_patt="^chr|Ecoli|Human_Ad5"; #pattern to define host genome (used for read clustering)
max_gap_len=ifelse(is.na(args_v["max_gap_len"]), 4000, as.numeric(args_v["max_gap_len"]) ); #maximium gap size (bp) to cluster reads
default_min_cluster_size=ifelse(is.na(args_v["default_min_cluster_size"]), 4700, as.numeric(args_v["default_min_cluster_size"]) ); #default minimal cluster size (aav genome size)


mkdir_if_not_exist(out_readCluster_f)
counts_names=paste("Num",samples, sep="_")
pval_names=paste("Pval",samples, sep="_")
enrichment_names=paste("Enrich",samples, sep="_")

ref_sizes_d=read.table(ref_sizes_f, header=F, sep="\t", quote="", stringsAsFactors=F)
names(ref_sizes_d)=c("reference","ref_size")

##1, combine all read counts
for (sample in samples){
	inf=paste0(read_in_dir, sample, ".read.map.coor.txt"); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, check.names=F)
	d=d[grepl(hostGenome_patt, d$rname), ]
	if("readsnum" %in% names(d)){
		names(d)[names(d)=="readsnum"]=paste0("Num_",sample)
	}
	if(sample==samples[1]){
		comb_d=d
	}else{
		comb_d=merge(comb_d,d,all=T)
	}
}
comb_d[,counts_names][is.na(comb_d[,counts_names])]=0
comb_d[,counts_names][comb_d[,counts_names]>1]=1 #use unique reads
comb_d=comb_d[order(comb_d$rname, comb_d$posfr, comb_d$posto),]
if_new_cluster=rep(T,nrow(comb_d))
if_cluster=comb_d$posfr[-1] - comb_d$posto[-nrow(comb_d)]<max_gap_len & (comb_d$rname[-1] == comb_d$rname[-nrow(comb_d)] )
if_new_cluster[-1][if_cluster]=F
table(if_new_cluster)
comb_d$clusterID=cumsum(if_new_cluster)
range(comb_d$clusterID)

cluster_d=data.frame(clusterID=sort(unique(comb_d$clusterID))) 
cluInfo_l=tapply(1:nrow(comb_d), comb_d$clusterID, function(rowIDs){
	c(comb_d$rname[rowIDs[1]], range(comb_d[rowIDs, c("posfr","posto")]) )
})[as.character(cluster_d$clusterID)]
cluInfo_m=matrix(unlist(cluInfo_l),nrow=3)
cluster_d$rname=cluInfo_m[1,]
cluster_d$clu_from=as.numeric(cluInfo_m[2,])
cluster_d$clu_to=as.numeric(cluInfo_m[3,])
for(counts_name in counts_names){
	cluster_d[, counts_name]=unlist( tapply(comb_d[,counts_name], comb_d$clusterID, sum, na.rm=T)[as.character(cluster_d$clusterID)] )
}

##calculate P-value for each cluster
table(cluster_d$rname)
table(comb_d$rname)
#calculate total read counts for each chromosome, add data to ref_sizes_d
ref_sizes_d=ref_sizes_d[ref_sizes_d$reference %in% comb_d$rname, ]
for(counts_name in counts_names){
	ref_sizes_d[,counts_name]=unlist(tapply(comb_d[,counts_name], comb_d$rname, sum)[as.character(ref_sizes_d$reference)])
}
#calculate p-values based on poisson distribution
cluster_d$cluster_size=cluster_d$clu_to-cluster_d$clu_from+1
table(cluster_d$cluster_size>0)


for(i in 1:length(counts_names)){
	counts_name=counts_names[i]
	pval_name=pval_names[i]
	used_cluster_size=ifelse(cluster_d$cluster_size<default_min_cluster_size, default_min_cluster_size, cluster_d$cluster_size)
	lambdas=(ref_sizes_d[,counts_name] /ref_sizes_d$ref_size) [match(cluster_d$rname, ref_sizes_d$reference)] * used_cluster_size
	if_not_0=cluster_d[,counts_name]>0
	cluster_d[,pval_name]=1
	cluster_d[, enrichment_names[i]]=0
	cluster_d[if_not_0,pval_name]=ppois(cluster_d[if_not_0,counts_name], lambda=lambdas[if_not_0], lower=FALSE)
	cluster_d[if_not_0, enrichment_names[i]]=round(cluster_d[,counts_name]/lambdas,1)[if_not_0]
}

cluster_d$name=paste(sprintf(paste0("clu%0",nchar(nrow(cluster_d)),"d"), 1:nrow(cluster_d) ), 
	"|n=",apply(cluster_d[counts_names],1,max,na.rm=T),
	"|p=", format(apply(cluster_d[pval_names],1,min,na.rm=T), scientific=T, digits=2), sep="" )

#output cluster_d
cluster_d=re_format_tb(cluster_d, back_header=c(counts_names,pval_names,enrichment_names) )
write.table(cluster_d, file=out_readCluster_f, col.names=T, row.names=F, sep="\t", quote=F)

##create and output bed format file
scores=apply(cluster_d[enrichment_names],1,max)/10
scores[scores>1000]=1000
bed_d=data.frame(chr=cluster_d$rname, from=cluster_d$clu_from-1, to=cluster_d$clu_to, name=cluster_d$name,  sco=scores , strand=".")
write.table(bed_d, file=out_readCluster_bedf, col.names=F, row.names=F, sep="\t", quote=F)

save.image(out_data_img_f)
