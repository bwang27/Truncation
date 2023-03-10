#for RNA-seq differential expression (DE) analysis, calculate rpkm change and p-value (fisher's exact test, chi square test, DEseq2)
#3/3/2014, can allow parallal computing when nCores was set >1
#3/31/2015, add function to use gene expression normalization method in edgeR (GexNorm_method)
#9/1/2015, allow analysis using DESeq2
#3/16/2016, draw clustering of raw samples based on RPKM or RPM, identify sample outliers. (what2do=clust_rawSample_exp)
#3/17/2016, can calculate gene expression with a phenotype (eg., tumor size), what2do=calCorWithPhenotype
#1/29/2019, can set fold change to 1, if P value>setFC1_p_cut
#7/15/2019, modified to allow using raw P-value (instead of default adjusted P-value) in DESeq2, set useP_type="pvalue"
#6/7/2021, add function to allowing selection of significantly regulated genes using an higher_avg_RPKM_cutoff (between ) cutoff (eg. for the comparion two groups of samples, the higher one of the average rpkm within a group should be higher than a cutoff)
#7/21/2021, added function to output a table (topGeneList_tb) of top DEGs for each treatment
#11/10/2021, for topGeneList_tb, also output the L2FC values etc as details

#available what2do: all (default), formatGeneTb, clust_rawSample_exp, calCorWithPhenotype, make_plots, rawSample_ScatterPlot

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}
source("../SharedCodes/Rfunc.inc.R")


#parameters from command line 
what2do<- "all"
if(!is.na(args_v["what2do"])){what2do=unlist(strsplit(as.character(args_v["what2do"])," "))}
out_prefix<- ifelse(is.na(args_v["out_prefix"]), "", args_v["out_prefix"])
p_cut <- ifelse(is.na(args_v["p_cut"]), 0.01, as.numeric(args_v["p_cut"])); 
setFC1_p_cut <- ifelse(is.na(args_v["setFC1_p_cut"]), 0.5, as.numeric(args_v["setFC1_p_cut"]));  #set fold change to 1, if P value>setFC1_p_cut
if_correct_pval<- ifelse(is.na(args_v["if_correct_pval"]), FALSE, args_v["if_correct_pval"]=="1" ); if_correct_pval
pval_fun<- ifelse(is.na(args_v["pval_fun"]), "DESeq2", args_v["pval_fun"]); #fisher chisq none DESeq2
foldchange_cut <- ifelse(is.na(args_v["foldchange_cut"]), 1.5, as.numeric(args_v["foldchange_cut"])); foldchange_cut
IfCal_RPKM <- ifelse(is.na(args_v["IfCal_RPKM"]), TRUE, args_v["IfCal_RPKM"]=="1" ); IfCal_RPKM
IfCal_RPM <- ifelse(is.na(args_v["IfCal_RPM"]), TRUE, args_v["IfCal_RPM"]=="1" ); IfCal_RPM  ##calculate RPM, eg. for small RNA-seq; 3'READS set T, otherwise set F
gene_model<- ifelse(is.na(args_v["gene_model"]), "refseqcds", args_v["gene_model"]); gene_model ##club refseq enscds repeat LTR Rcluster pAR2club
tlRnum_using<-ifelse(is.na(args_v["tlRnum_using"]), "colSum", args_v["tlRnum_using"]); tlRnum_using #"colSum", using read number column sum (Read mapped to all gene models) as total read number; "uniqueMapped", using uniquely mapped reads number
comb_sample_l=list()
if(!is.na(args_v["comb_sample_l"])){ comb_sample_l<- string2list(args_v["comb_sample_l"]) }; comb_sample_l
nCores <- ifelse(is.na(args_v["nCores"]), 1, as.numeric(args_v["nCores"])); 
higher_avg_RPKM_cutoff <- ifelse(is.na(args_v["higher_avg_RPKM_cutoff"]), 1, as.numeric(args_v["higher_avg_RPKM_cutoff"])); higher_avg_RPKM_cutoff

geno<- ifelse(is.na(args_v["geno"]), "hg19", args_v["geno"])
study_name<- ifelse(is.na(args_v["study_name"]), NA, args_v["study_name"]) ##project name
output_root<- ifelse(is.na(args_v["output_root"]), paste("05.Gene_DE/",study_name,"/",sep=""), args_v["output_root"]) 
out_plot_root=paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/GexPlot/",sep="");

if( !is.na(args_v["gene_rnum_f"]) ){ 	gene_rnum_f=args_v["gene_rnum_f"] } #gene read counts file
gene_rnum_f;
gene_len_var <- ifelse(is.na(args_v["gene_len_var"]), NA, args_v["gene_len_var"]) #gene length variable (column name)
test_samples<- unlist(strsplit(as.character(args_v["test_samples"])," ")) #must set a value 
ref_samples<- unlist(strsplit(as.character(args_v["ref_samples"])," ")) #must set a value 
if(!is.na(args_v["other_samples"])){ 
	other_samples<- unlist(strsplit(as.character(args_v["other_samples"])," ")) 
}
if( any( c( is.na(test_samples) , is.na(ref_samples) ) ) ){		print ("test_samples or ref_samples is NA!"); quit(); 	}
sel_gex_ids=1:length(test_samples)
if(!is.na(args_v["sel_gex_ids"])){ 
	sel_gex_ids<- as.numeric( unlist(strsplit(args_v["sel_gex_ids"]," ") ) ) 
}

GexNorm_method<-ifelse(is.na(args_v["GexNorm_method"]), "totalMillionRead", args_v["GexNorm_method"]);GexNorm_method; #totalMillionRead(RPM or RPKM), TMM (weighted trimmed mean of M-values, used on edgeR)

out_file <- paste(output_root, out_prefix, gene_model,
	ifelse( GexNorm_method=="totalMillionRead","",paste(".",GexNorm_method,sep="") ),
	".",pval_fun,".p",sep="");
if( !is.na(args_v["out_file"]) ){ out_file=args_v["out_file"] }
g_ano_f<-ifelse(is.na(args_v["g_ano_f"]), NA, args_v["g_ano_f"])
anof_idname<-ifelse(is.na(args_v["anof_idname"]), NA, args_v["anof_idname"])
rnum_idname<-ifelse(is.na(args_v["rnum_idname"]), "refseqid", args_v["rnum_idname"])
geneSym_name<-ifelse(is.na(args_v["geneSym_name"]), "gene_symbol", args_v["geneSym_name"])
geneType_anof=paste("../../club/21transcriptAno/",geno,".transcript.ano", sep="") #have gene type annotation
if(!is.na(args_v["geneType_anof"])){
	geneType_anof=ifelse(args_v["geneType_anof"] %in% c("NA","na"), NA, args_v["geneType_anof"])
}
geneType_ano_idName=ifelse(is.na(args_v["geneType_ano_idName"]), "name", args_v["geneType_ano_idName"])

#for what2do=formatGeneTb
formatGeneTbName<-ifelse(is.na(args_v["formatGeneTbName"]), "", args_v["formatGeneTbName"])
L2FC_scale_Max=ifelse(is.na(args_v["L2FC_scale_Max"]), 2, as.numeric(args_v["L2FC_scale_Max"]) )
filter_header=ifelse(is.na(args_v["filter_header"]), NA, args_v["filter_header"]) #used for filter the table when plotting the clustering plot
filter_text=ifelse(is.na(args_v["filter_text"]), NA, unlist(strsplit(as.character(args_v["filter_text"])," "))) #only plot rows with these filter_text values
max_geneList_num=ifelse(is.na(args_v["max_geneList_num"]), 100, as.numeric(args_v["max_geneList_num"]) )

comment_f=paste(out_file,"comment.tbl",sep=".")

define_parallel_fun(nCores=nCores)
print(paste("nCores=",nCores))
mkdir_if_not_exist(out_plot_root)

##common setup:
gexch_types <- c("up","dn","nc","na")
if(!("other_samples" %in% ls())){	other_samples<-NULL }
if(!is.na(args_v["all_sample_name"])){ 
	all_sample_name<- unlist(strsplit(as.character(args_v["all_sample_name"])," ")) 
}else if("sample_grps" %in% ls()){
	all_sample_name=unique(unlist(sample_grps));
}else{
	all_sample_name <- unique(c(ref_samples,test_samples,other_samples))
}
print(all_sample_name)

reads_num_names <- paste("num_",all_sample_name,sep="")
rpkm_names <- paste(ifelse(IfCal_RPM,"rpm_","rpkm_"), all_sample_name,sep="")
samplePairs=paste(test_samples,ref_samples,sep="_")
gchtp_names <- paste("GexType_",test_samples,"_",ref_samples,sep=""); gchtp_names
l2gch_names <- paste("L2FC_",test_samples,"_", ref_samples,sep="")
pval_names <- paste("SLog10P_",test_samples, "_", ref_samples, sep="")
adjPval_names <- paste("adjSLog10P_",test_samples, "_", ref_samples, sep="")
corr_p_cut <- p_cut/ifelse(if_correct_pval ,nrow(cmb_d),1) #Bonferroni corrected p value cut
isexpP_names <- paste("expP_",all_sample_name,sep="") #is expressed p value (poisson distribution)
isexp_names <- paste("exp_",all_sample_name,sep="")

if(pval_fun %in% c("DESeq2") ){ ##automatic define P-value type
	useP_type="padj"
}else{
	useP_type="pvalue"
}
if(!is.na(args_v["useP_type"])){ useP_type=args_v["useP_type"] }
if(useP_type %in% c("padj") ){ ##automatic define P-value type
	use_pval_names=adjPval_names
}else{
	use_pval_names=pval_names
}

min_TotalreadNum <- ifelse(is.na(args_v["min_TotalreadNum"]), length(reads_num_names), as.numeric(args_v["min_TotalreadNum"])); 


#for what2do=calCorWithPhenotype
# phenotypeName="tumorVol"
# phenotypeSamples=unlist(strsplit("Veh_1.1 Veh_1.2 Veh_1.3 Veh_1.7 Veh_1.8 PTC596_5mg_kg_2.1 PTC596_5mg_kg_2.4 PTC596_5mg_kg_2.9 PTC596_10mg_kg_4.1 PTC596_10mg_kg_4.3 PTC596_10mg_kg_4.9 PTC596_15mg_kg_6.4 PTC596_15mg_kg_6.8 PTC596_15mg_kg_6.10"," ")) 
# phenotypeVals=as.numeric(unlist(strsplit("1536 1797 1447 1856 2016 1744 1888 1712 961 950 1090 562 377 325"," ")) )
phenotypeName<-ifelse(is.na(args_v["phenotypeName"]), NA, args_v["phenotypeName"])
if(!is.na(args_v["phenotypeSamples"])){ phenotypeSamples=unlist(strsplit(as.character(args_v["phenotypeSamples"])," ")) }
if(!is.na(args_v["phenotypeVals"])){ phenotypeVals=as.numeric(unlist(strsplit(as.character(args_v["phenotypeVals"])," "))) }

figure_out_format=ifelse(capabilities()['png'],"png","pdf")
figure_size_fac=ifelse(figure_out_format=="png",1,1/72) #convert pixel in png to inch in pdf

################### RUN
if(any(what2do %in% c("all")) ){
	print(paste("read gene_rnum_f:",gene_rnum_f))
	gene_rnum_D <- read.table(gene_rnum_f, sep="\t", header=T, comment.char="", quote='')
	gene_rnum_D <- gene_rnum_D [,apply(!is.na(gene_rnum_D),2,any)]

	###combine with gene anotation
	if( !("g_ano_f" %in% ls()) ){g_ano_f<-NA}
	if(!is.na(g_ano_f)){
		print(paste("read g_ano_f: ",g_ano_f))
		g_ano_d <- read.table(g_ano_f, header=T, sep="\t", quote="", comment.char="")
		g_ano_d<-g_ano_d[setdiff(names(g_ano_d),c("exon_starts","exon_ends","refseqid.1"))]
		if("gene_Biotype" %in% names(g_ano_d) ){
			g_ano_d=g_ano_d[order(!is.na(g_ano_d$gene_Biotype) & g_ano_d$gene_Biotype %in% c("protein_coding"), decreasing=T ), ] #put protein coding genes on top
		}
		if( !is.na(anof_idname) & !is.na(rnum_idname) ){
			#cmb_d <- cbind(g_ano_d,gene_rnum_D[match(g_ano_d[,anof_idname], gene_rnum_D[,rnum_idname]),])
			cmb_d <- gene_rnum_D
			update_names=setdiff(names(g_ano_d),names(gene_rnum_D))
			cmb_d[,update_names]=g_ano_d[match(gene_rnum_D[,rnum_idname], g_ano_d[,anof_idname]), update_names]
		}else{
			if(nrow(g_ano_d)!=nrow(gene_rnum_D)){
				print("Error: row number in gene_rnum_D do not match g_ano_d!")
				quit()
			}else{
				cmb_d <- cbind(g_ano_d, gene_rnum_D[setdiff(names(gene_rnum_D),names(g_ano_d))])
			}
		}
	}else{
		cmb_d <- gene_rnum_D
	}
	rm ("gene_rnum_D","g_ano_d")
	save.image( paste(out_file,".img",sep="") )

	#combine reads number
	raw_num_names <- grep("num_",names(cmb_d),value=T)
	cmb_d[raw_num_names][is.na(cmb_d[raw_num_names])] <- 0
	for(i in 1:length(reads_num_names)){
		reads_num_name=reads_num_names[i]
		if( ! (reads_num_name %in% names(cmb_d)) ){
			if(all_sample_name[i] %in% names(comb_sample_l)){
					combine_frs <- paste('num_',comb_sample_l[[all_sample_name[i]]],sep='')
			}else{
				search_patt=paste("num_",".*",all_sample_name[i],sep="")
				combine_frs=grep(search_patt,names(cmb_d),value=T)
			}
			print (paste(c("combine read number for sample",reads_num_name, " from", combine_frs), collapse=" " ) )
			cmb_d[reads_num_name] <- rowSums(cmb_d[combine_frs], na.rm=T)
		}
	}
	cmb_d[reads_num_names][is.na(cmb_d[reads_num_names])] <- 0

	##remove redundancy in cmb_d
	#sort with read number; remove rows with duplicated rnum_idname
	cmb_d=cmb_d[order(rowSums(cmb_d[reads_num_names]),decreasing=T),]; nrow(cmb_d)
	if(!is.na(rnum_idname)){
		cmb_d=cmb_d[!duplicated(cmb_d[,rnum_idname]),]; nrow(cmb_d)
	}
	#remove transcripts with too small number of reads
	ifExpressed=rowSums(cmb_d[reads_num_names])>=min_TotalreadNum; sum(ifExpressed)
	cmb_d=cmb_d[ifExpressed,]; nrow(cmb_d)

	##calculate scaling factor if GexNorm_method=="TMM"
	if(GexNorm_method=="TMM"){
		library('edgeR')
		scaling_facs=calcNormFactors(cmb_d[,reads_num_names])
	}else{
		scaling_facs=rep(1,length(reads_num_names))
	}
	names(scaling_facs)=reads_num_names

	#calculate t_readsnum
	if(tlRnum_using=="colSum"){
		t_readsnum <- apply(cmb_d[reads_num_names],2,sum, na.rm=T); 
	}else if(tlRnum_using=="uniqueMapped"){
		t_readsnum <- tl_mapped_readsnum[all_sample_name]
		names(t_readsnum) <- reads_num_names
	}
	#calculate rpkm or rpm
	t_readsnum=round( t_readsnum*(scaling_facs) )
	t_readsnum
	if(gene_model=="club"){ 
		gene_len_var <- paste("_len",sep="")
		print (paste("gene_len_var=",gene_len_var),quote=F)
	}else if(gene_model=="repeat"){ 
		if("Trnum_cut" %in% ls()){ #filter low read number rows
			print (paste("original rows:", nrow(cmb_d)))
			cmb_d <- cmb_d[rowSums(cmb_d[reads_num_names],na.rm=T)>=Trnum_cut,]
			print (paste("after removing rows with total read number<",Trnum_cut, nrow(cmb_d), "rows remaining."))
		}
		cmb_d[,gene_len_var] <- abs(cmb_d$repTo-cmb_d$repFr)
	}
	if( is.null(cmb_d[1,rpkm_names[1]]) | IfCal_RPKM | IfCal_RPM){
		if(IfCal_RPM){
			cmb_d[rpkm_names] <- round(t( t(cmb_d[reads_num_names]) / (t_readsnum/1000000) ) , 3)
		}else{
			if(gene_len_var=="cds_len" & "transcript_len" %in% names(cmb_d)){
				if_change=!is.na(cmb_d$cds_len) & cmb_d$cds_len==0
				cmb_d$cds_len[if_change]=cmb_d$transcript_len[if_change]
			}
			cmb_d[rpkm_names] <- round(t( t(cmb_d[reads_num_names]/(cmb_d[,gene_len_var]/1000)) / (t_readsnum/1000000) ) , 3)
		}
	}

}

if(any(what2do %in% c("all")) ){
	if(pval_fun =="DESeq2"){
		library("DESeq2")
		library("BiocParallel")
		register(MulticoreParam(nCores))
		for(i in 1:length(test_samples)){
			samp_pair <- c(test_samples[i], ref_samples[i])
			names(samp_pair)=c("test","ctrl")

			gchtp_name <- gchtp_names[i]
			l2gch_name <- l2gch_names[i]
			print (samp_pair)
			reads_num_2names <- paste("num_",samp_pair,sep="")
			#build countData and colData used for DESeqDataSetFromMatrix
			condition=NULL
			DeSeq_samples=NULL
			for(sample_type in names(samp_pair) ){
				sample1=samp_pair[sample_type]
				if(sample1 %in% names(comb_sample_l) ){
					DeSeq_samples=c(DeSeq_samples, comb_sample_l[[sample1]])
					condition=c(condition, rep(sample_type, length(comb_sample_l[[sample1]])) )
				}else{
					DeSeq_samples=c(DeSeq_samples, sample1)
					condition=c(condition, sample_type )
				}
			}
			coldata=data.frame(condition=condition); rownames(coldata)=make.unique(DeSeq_samples)
			countData=cmb_d[paste("num_",DeSeq_samples,sep="")]
			if_no_counts=rowSums(countData)<1
			colnames(countData)=make.unique(DeSeq_samples)
			rownames(countData)=cmb_d[,rnum_idname]

			dds <- DESeqDataSetFromMatrix(countData = countData[!if_no_counts,], colData = coldata, design = ~ condition )
			dds <- DESeq(dds, parallel=TRUE)
			res <- results(dds)
			print(summary(res))
			save.image( paste(out_file,".img",sep="") )

			cmb_d[,l2gch_name]=NA
			cmb_d[!if_no_counts,l2gch_name] = res$log2FoldChange

			if(useP_type %in% c("padj") ){
				cmb_d[,adjPval_names[i]]=NA
				cmb_d[!if_no_counts,adjPval_names[i]]= -sign(res$log2FoldChange) * log10(res$padj)
			}else{
				cmb_d[,pval_names[i]]=NA
				cmb_d[!if_no_counts,pval_names[i]]= -sign(res$log2FoldChange) * log10(res$pvalue)
			}
		}
		cmb_d[l2gch_names][ is.na(cmb_d[use_pval_names]) | abs(cmb_d[use_pval_names])< abs(log10(setFC1_p_cut))  ] = 0
		write.csv(res, file="DeSeq2Res.csv", col.names=T, row.names=F, quote=F, append=T)
		rm("dds","res")


	}else{ #fisher's exact test
		if(!is.null(test_samples)){
			#fold change:
			cmb_d[l2gch_names] <- round(log2( cmb_d[paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),test_samples,sep="")] / cmb_d[paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),ref_samples,sep="")] ) ,3)

			#fisher's exact test/chi2 test,  for gene expression change of every pair
			for(i in 1:length(test_samples)){
				samp_pair <- c(test_samples[i], ref_samples[i])
				gchtp_name <- gchtp_names[i]
				l2gch_name <- l2gch_names[i]
				print (samp_pair)
				reads_num_2names <- paste("num_",samp_pair,sep="")
				numtb <- cbind(cmb_d[reads_num_2names], t(t_readsnum[reads_num_2names]-t(cmb_d[reads_num_2names])))
				pval_name <- pval_names[i]
				numtb=as.matrix(numtb)
				#cmb_d[pval_name] <- apply(numtb,1,paste("do",pval_fun,"test",sep="_") )

			 	if_do_test=rowSums(numtb[,1:2],na.rm=T)>=7; table(if_do_test)
			 	if(pval_fun !="none"){
				 	cmb_d[pval_name]=0
				 	cmb_d[if_do_test,pval_name] <- unlist(myApply((1:nrow(numtb))[if_do_test], function(row_i){
				 		do.call(paste("do",pval_fun,"test",sep="_"), list(count_vec=round(numtb[row_i,])) )
				 	 } ))
			 	}
			}
			cmb_d[l2gch_names][ is.na(cmb_d[pval_names]) | abs(cmb_d[pval_names])< abs(log10(setFC1_p_cut))  ] = 0
		}

	}




}


