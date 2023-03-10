
##additional code for 27ana_junc_map_info.R (used in interactive mode of R)
##1, study the distribution of read number, log2 ratio and P value for regulated AS events
##2, study splice site strength for regulated exon and their flanking exons
##3, study relationship between motif and splicing or gene expression regulation
##4, study 5'ss -4 and -3 position di-nucleotide frequency of up-regulated exons
##5, study 5'ss derived from transposable elements
##6.1, study splice site annotation vs. regulation (cryptic exon vs. annotated exon)
##6.2 study genes with cryptic exon inclusion and gene expression changes
##6.3, study basal PSI for different groups of exons (eg. iExon up regulated vs. annotated exon NC ... ...)
##7, output a table of junction read counts for GEO submission
##8, study FDR based on DMSO samples in HT-RNAseq


##1, study the distribution of read number, log2 ratio and P value for regulated AS events
#group by GAgt ss5 of up-regulated exon
#for ana_unit=="geneBlock"
regu_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.1.Ch0.26.tbl.regu.tbl",sep="")
regu_d=read.table(regu_f, header=T, sep="\t", comment.char="", quote="", stringsAsFactors=F)
use_P_names=adjpval_names

out_img_f=paste(regu_f,".density.pdf",sep="")
pdf(out_img_f, width=12, height=9)
layout(matrix(1:12,3))
for(i in 1:length(use_P_names)){
	ifsel=!is.na(regu_d[,ReguType_names[i]]) & regu_d[,ReguType_names[i]]=="UP" & !is.na(regu_d$region_ano) & regu_d$region_ano=="exon"; sum(ifsel)
	regu_d2=regu_d[ifsel, ]
	regu_d2$type="others"
	regu_d2$type[grep("GAgt",regu_d2$endSS_seq)]="GAgt"

	rawReadNum_names=paste("num_",unlist(sample_repl_l[sample_compare_matrix[,i]]), sep="")
	regu_d2$Avg_readNum=rowMeans(regu_d2[,rawReadNum_names])
	regu_d2$SS=log10( -log10(regu_d2[,use_P_names[i]]) )
	draw_ecdf_density(regu_d2, "Avg_readNum", "type", c("others","GAgt"), mycols=c("black","red"), 
		method = "density", xmin=0, xmax=150 )
	draw_ecdf_density(regu_d2, Log2Ratio_names[i], "type", c("others","GAgt"), mycols=c("black","red"), 
		method = "density", xmin=0, xmax=6 )
	draw_ecdf_density(regu_d2, "SS", "type", c("others","GAgt"), mycols=c("black","red"), xlab = paste("log10 (-log10(",use_P_names[i],"))"), 
		method = "density", xmin=0, xmax=3 )
	#draw_ecdf_density(regu_d2, use_P_names[i], "type", c("others","GAgt"), mycols=c("black","red"), 
	#	method = "density", xmin=0, xmax=0.1 )
}
dev.off()


##2, study splice site strength for regulated exon and their flanking exons
AS_reg_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl",sep="")
AS_reg_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.txt",sep="") #regulated only
AS_reg_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.allReguType.txt",sep="") #all exons with regulation type
all_gblock_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/exons.tbl",sep="") #used for additional annotation
out_reg_f=paste(AS_reg_f,".novelexon.ss_score.txt",sep="")
out_all_ss_f=paste(AS_reg_f,".all_intron_ss.txt",sep="")
geno="hg19"
all_ss5_score_f=paste("../../AS/ss_score/1all_ss5/GTss5.score.out",sep="")
transcript_ss_sco_f=paste("../../club/11ss_score/",geno,".refflat.redundant.ss.sco.txt",sep="")

all_ss5_score_d=read.table(all_ss5_score_f, header=T, sep="\t")
transcript_ss_sco_d=read.table(transcript_ss_sco_f, header=T, sep="\t", stringsAsFactors=F)
transcript_ss_sco_d$intron_id=paste(transcript_ss_sco_d$refseqid, ":E",transcript_ss_sco_d$introni,sep="")
transcript_ss_sco_d$intron_size=abs(transcript_ss_sco_d$ss5-transcript_ss_sco_d$ss3)
if_duplicated=duplicated(paste(transcript_ss_sco_d$contig,transcript_ss_sco_d$strand, transcript_ss_sco_d$ss5,transcript_ss_sco_d$ss3), sep=":"); table(if_duplicated)
transcript_ss_sco_d_nonRedundant=transcript_ss_sco_d[!if_duplicated,]

AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
required_headers=c("startSS_seq","endSS_seq")
if(!all(required_headers %in% names(AS_reg_d))){
	all_gblock_d=read.table(all_gblock_f, header=T, sep="\t", stringsAsFactors=F, quote="", comment.char="")
	update_headers=setdiff(required_headers, names(AS_reg_d))
	intersect(names(all_gblock_d),names(AS_reg_d))
	all_gblock_d$Gblock_id2=paste(all_gblock_d$contig, all_gblock_d$strand, all_gblock_d$start_pos, all_gblock_d$end_pos, sep=":")
	table(AS_reg_d$Gblock_id2 %in% all_gblock_d$Gblock_id2)
	AS_reg_d[update_headers]=all_gblock_d[match(AS_reg_d$Gblock_id2 ,all_gblock_d$Gblock_id2), update_headers]
}
table(AS_reg_d$region_ano)

#use ##6.3 defined AS_reg_d1
AS_reg_d2=AS_reg_d1

#study up-regulated novel exon, 
ifsel=AS_reg_d$region_ano=="exon" # & grepl("^n",AS_reg_d$endPosAno) & grepl("^n",AS_reg_d$startPosAno); table(ifsel) #usually not used
ifsel=AS_reg_d$region_ano=="exon" & grepl("^n",AS_reg_d$endPosAno) & grepl("^n",AS_reg_d$startPosAno); table(ifsel)
AS_reg_d2=AS_reg_d[ifsel,]
ifsel2=substr(AS_reg_d2$endSS_seq, 11, 12)=="gt" & substr(AS_reg_d2$startSS_seq, 49, 50)=="ag" & grepl("ss5$",AS_reg_d2$startPosAno); table(ifsel2)
AS_reg_d2=AS_reg_d2[ifsel2,]
ifsel3=AS_reg_d2$reg_len<200; table(ifsel3)
AS_reg_d2=AS_reg_d2[ifsel3,]

ss5_seqs=toupper(substr(AS_reg_d2$endSS_seq, 8, 16))
AS_reg_d2$ss5_sco=all_ss5_score_d$ss5_sco[match(ss5_seqs,all_ss5_score_d$ss5_seq)]
table(is.na(AS_reg_d2$ss5_sco))
ss3_seqs=toupper(substr(AS_reg_d2$startSS_seq, 31, 53))
AS_reg_d2$ss3_sco=cal_ss3_score(ss3_seqs, prog_dir="../../AS/ss_score/")
AS_reg_d2$if_psi_exon=grepl("^n",AS_reg_d2$endPosAno) & grepl("^n",AS_reg_d2$startPosAno) & grepl("ss5$",AS_reg_d2$startPosAno); table(AS_reg_d2$if_psi_exon)

transc_id_header<- ifelse(is.na(args_v["transc_id_header"]), "refseqid", args_v["transc_id_header"])
AS_reg_d2$intron_id=gsub("^n-d|ss5$|ss3$","",AS_reg_d2$startPosAno)
ifadd_transcript=!grepl(":",AS_reg_d2$intron_id); sum(ifadd_transcript)
AS_reg_d2$intron_id[ifadd_transcript]=paste(AS_reg_d2[ifadd_transcript,transc_id_header]  ,AS_reg_d2$intron_id[ifadd_transcript],sep=":")
AS_reg_d2=AS_reg_d2[AS_reg_d2$intron_id %in% transcript_ss_sco_d$intron_id, ]; nrow(AS_reg_d2)
AS_reg_d2$intron_i=as.numeric(gsub("^.*:E","",AS_reg_d2$intron_id))
table(AS_reg_d2$if_psi_exon, AS_reg_d2$intron_i==1)
AS_reg_d2=AS_reg_d2[!(!AS_reg_d2$if_psi_exon & AS_reg_d2$intron_i==1), ]
AS_reg_d2$ups_intron_i=AS_reg_d2$intron_i
AS_reg_d2$ups_intron_i[!AS_reg_d2$if_psi_exon]=AS_reg_d2$intron_i[!AS_reg_d2$if_psi_exon]-1
table(AS_reg_d2$ups_intron_i)
AS_reg_d2$ups_intron_id=paste0(sub("\\d+$","", AS_reg_d2$intron_id),AS_reg_d2$ups_intron_i)

AS_reg_d2[,c("ups_ss5","dns_ss3", "ups_ss5_sco","dns_ss3_sco","intron_size")]=NA
AS_reg_d2[AS_reg_d2$if_psi_exon,c("ups_ss5","dns_ss3", "ups_ss5_sco","dns_ss3_sco","intron_size")]=transcript_ss_sco_d[match(AS_reg_d2$intron_id[AS_reg_d2$if_psi_exon], transcript_ss_sco_d$intron_id), c("ss5","ss3", "ss5_sco","ss3_sco","intron_size")]
AS_reg_d2[!AS_reg_d2$if_psi_exon,c("dns_ss3", "dns_ss3_sco")]=transcript_ss_sco_d[match(AS_reg_d2$intron_id[!AS_reg_d2$if_psi_exon], transcript_ss_sco_d$intron_id), c("ss3", "ss3_sco")]
AS_reg_d2[!AS_reg_d2$if_psi_exon,c("ups_ss5", "ups_ss5_sco")]=transcript_ss_sco_d[match(AS_reg_d2$ups_intron_id[!AS_reg_d2$if_psi_exon], transcript_ss_sco_d$intron_id), c("ss5", "ss5_sco")]
AS_reg_d2[!AS_reg_d2$if_psi_exon,c("intron_size")]=transcript_ss_sco_d[match(AS_reg_d2$ups_intron_id[!AS_reg_d2$if_psi_exon], transcript_ss_sco_d$intron_id), c("intron_size")] + 
	transcript_ss_sco_d[match(AS_reg_d2$intron_id[!AS_reg_d2$if_psi_exon], transcript_ss_sco_d$intron_id), c("intron_size")]+ 
	AS_reg_d2$reg_len[!AS_reg_d2$if_psi_exon]

AS_reg_d2$iExon_ss5_m2m1=substr(AS_reg_d2$endSS_seq, 9, 10)
table(AS_reg_d2$iExon_ss5_m2m1, AS_reg_d2[,ReguType_names[1]])


##add iExon regulation data to transcript_ss_sco_d; can be used to study neighboring splice site sequence
add_regulation_names=ReguType_names
sort_regu_name=add_regulation_names[4]
transcript_ss_sco_d[,c(add_regulation_names,"iExon_ss5_m2m1")]=AS_reg_d2[match(transcript_ss_sco_d$intron_id, AS_reg_d2$intron_id), c(add_regulation_names,"iExon_ss5_m2m1")]
transcript_ss_sco_d[,c(add_regulation_names,"iExon_ss5_m2m1")][is.na(transcript_ss_sco_d[,c(add_regulation_names,"iExon_ss5_m2m1")])]="na"
transcript_ss_sco_d=transcript_ss_sco_d[order(match(transcript_ss_sco_d[,sort_regu_name],ReguTypes)), ]
if_duplicated=duplicated(paste(transcript_ss_sco_d$contig,transcript_ss_sco_d$strand, transcript_ss_sco_d$ss5,transcript_ss_sco_d$ss3), sep=":"); table(if_duplicated)
transcript_ss_sco_d_nonRedundant=transcript_ss_sco_d[!if_duplicated, ]
if_gt_ag_ss=substr(transcript_ss_sco_d_nonRedundant$ss5_seq, 4, 5)=="GT" & substr(transcript_ss_sco_d_nonRedundant$ss3_seq, 19, 20)=="AG"; table(if_gt_ag_ss)
table(if_gt_ag_ss, transcript_ss_sco_d_nonRedundant[,sort_regu_name])
transcript_ss_sco_d_nonRedundant=transcript_ss_sco_d_nonRedundant[if_gt_ag_ss,]
transcript_ss_sco_d_nonRedundant$upsExon_ss5_m2m1=substr(transcript_ss_sco_d_nonRedundant$ss5_seq, 2, 3)
table(transcript_ss_sco_d_nonRedundant$upsExon_ss5_m2m1, transcript_ss_sco_d_nonRedundant[,sort_regu_name], transcript_ss_sco_d_nonRedundant$iExon_ss5_m2m1=="GA")


##output AS_reg_d2
write.table(AS_reg_d2, file=out_reg_f, col.names=T, row.names=F, sep="\t", quote=F)
#output transcript_ss_sco_d_nonRedundant
write.table(transcript_ss_sco_d_nonRedundant, file=out_all_ss_f, col.names=T, row.names=F, sep="\t", quote=F)

##build reference parameters (annotated 5'ss, 3'ss score, intron/exon size)
ifsel=AS_reg_d$region_ano=="exon" & !grepl("^n",AS_reg_d$endPosAno) & !grepl("^n",AS_reg_d$startPosAno); table(ifsel) #annotated exons
Annotated_exon_d=AS_reg_d[ifsel,]
ifsel2=substr(Annotated_exon_d$endSS_seq, 11, 12)=="gt" & substr(Annotated_exon_d$startSS_seq, 49, 50)=="ag" & Annotated_exon_d$reg_len<200; table(ifsel2)
Annotated_exon_d=Annotated_exon_d[ifsel2,]
ss5_seqs=toupper(substr(Annotated_exon_d$endSS_seq, 8, 16))
Annotated_exon_d$ss5_sco=all_ss5_score_d$ss5_sco[match(ss5_seqs,all_ss5_score_d$ss5_seq)]
table(is.na(Annotated_exon_d$ss5_sco))
ss3_seqs=toupper(substr(Annotated_exon_d$startSS_seq, 31, 53))
Annotated_exon_d$ss3_sco=cal_ss3_score(ss3_seqs, prog_dir="../../AS/ss_score/")
##build ref_data_l for annotated intron sizes

ref_data_l=list(intron_size=transcript_ss_sco_d$intron_size [transcript_ss_sco_d$refseqid %in% Annotated_exon_d$refseqid] )

##study the upstream exon 5'ss -2 to -1 dinucleotide frequency
sel_ReguType_names=ReguType_names[c(2,4,6,8,9)]
sel_ReguType_names=ReguType_names
selReguTypes=c("UP","na","NC")
out_stats_f=paste(out_all_ss_f,".ss5_nt.stats.tbl",sep="")
write("5\'ss -2 to -1 position NT for upstream exon of iExon:", file=out_stats_f)
for (i in 1:length(sel_ReguType_names)){
	numtb=table(transcript_ss_sco_d_nonRedundant[,sel_ReguType_names[i]], transcript_ss_sco_d_nonRedundant$upsExon_ss5_m2m1)[selReguTypes, ]
	perctb=numtb/rowSums(numtb)
	write(sel_ReguType_names[i], file=out_stats_f, append=T)
	write.table(cbind(selReguTypes,numtb, selReguTypes, perctb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
}



#re-open
# AS_reg_d2=read.table(out_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
sco_names=c("ups_ss5_sco","ss3_sco","ss5_sco","dns_ss3_sco","intron_size")

##study the strength of ss for up-regulated exons: (for NVS_SM1, percentile)
ifsel=AS_reg_d2$ReguType_NVS_SM1_NVS_DMSO=="UP"
apply(AS_reg_d2[!ifsel,sco_names], 2, quantile, probs=0:20/20)
tb=apply(AS_reg_d2[ifsel,sco_names], 2, quantile, probs=0:20/20)
out_stats_f=paste(AS_reg_f,".ss_score.stats.txt",sep="")
write(paste("percentiles of ss score for up regulated exons (ReguType_NVS_SM1_NVS_DMSO)"), file=out_stats_f)
write.table(cbind(rownames(tb),tb), file=out_stats_f, quote=F, sep="\t", col.names=T, row.names=F, append=T)


##draw boxplot or CDF curve
#compare same splice site for different groups of exons

sco_names=c("ups_ss5_sco","ss3_sco","ss5_sco","dns_ss3_sco","intron_size","reg_len")

log_sco_names=c("intron_size","reg_len")
grp_vars=ReguType_names[c(2,4)]; grps=c("UP","annotated", "NC"); grp_cols=c("red", "black", "orange")
grp_vars="exon_type"; grps=c("Both:UP","Both:NC", "None:UP" ); grp_cols=c("blue","black", "red" ); ltys=1 ##from ##6.3
grp_vars="exon_type"; grps=c("Both:UP","Both:NC", "None:UP" ); grp_cols=c("black","gray", "black" ); ltys=c(3,2,1) ##from ##6.3

apply(AS_reg_d2[grp_vars],2,table)

out_img_f=paste(AS_reg_f,".exon.ss_score.UPvsNC.cdf.pdf",sep="")
pdf(out_img_f, width=11, height=7); layout(matrix(1:6,2,byrow=F))
for (sco_name in sco_names){
	annotated_sco_name=sub("ups_|dns_","",sco_name)
	for(grp_var in grp_vars){
		tmp_data=AS_reg_d2[, c(sco_name, grp_var)]
		# if(annotated_sco_name %in% names(Annotated_exon_d)){
		# 	tmp_data2=data.frame(Annotated_exon_d[,annotated_sco_name] )
		# }else{
		# 	tmp_data2=data.frame(ref_data_l[[annotated_sco_name]] )
		# }
		# names(tmp_data2)=sco_name
		# tmp_data2[,grp_var]="annotated"
		# tmp_data=rbind(tmp_data,tmp_data2, stringsAsFactors=F)
		median_vals=tapply(tmp_data[,sco_name], tmp_data[,grp_var], median, na.rm=T)[grps]
		draw_ecdf_density(tmp_data, sco_name, grp_var, grps, mycols=grp_cols, test_method = "wilcox.test", test_sample_pairs="adjacent",
			verticals=T, lwd=2, pch=NA, xlog=sco_name %in% log_sco_names, desc=paste(c("median=",round(median_vals,2)), collapse=" " ), ltys=ltys  ) #xmin=-5, xmax=15
		#draw median value lines
		if(sco_name %in% log_sco_names){ median_vals=log10(median_vals) }
		abline(v=median_vals, col=grp_cols, lty=ltys, lwd=2)
	}
}
dev.off()

out_img_f=paste(AS_reg_f,".novelexon.ss_score.UPvs.Annotated.cdf.pdf",sep="")
pdf(out_img_f, width=16, height=8); layout(matrix(1:8,2,byrow=T))
#compare to a common group of scores (eg. all ss in refseq)
ref_sco_l=list(all_ref_ss5=transcript_ss_sco_d_nonRedundant$ss5_sco,
	all_ref_ss3=transcript_ss_sco_d_nonRedundant$ss3_sco
)
compare_var_l=list(c("all_ref_ss5","ups_ss5_sco","ss5_sco"), c("all_ref_ss3","dns_ss3_sco","ss3_sco"))
var_cols=c("black","royalblue","red")
for(grp_var in grp_vars){
	for (sel_grp in grps){
		ifsel=!is.na(AS_reg_d2[,grp_var]) & AS_reg_d2[,grp_var]==sel_grp
		desc=paste(grp_var,"=",sel_grp)
		for(i in 1:length(compare_var_l)){
			compare_vars=compare_var_l[[i]]
			tmp_d=data.frame()
			for(compare_var in compare_vars){
				if(compare_var %in% names(ref_sco_l)){
					tmp_d1=data.frame(sco=ref_sco_l[[compare_var]], var=compare_var )
				}else{
					tmp_d1=data.frame(sco=AS_reg_d2[ifsel,compare_var], var=compare_var )
				}
				tmp_d=rbind(tmp_d,tmp_d1)
			}
			draw_ecdf_density(tmp_d, "sco", grpname="var", grpvals=compare_vars , mycols=var_cols, test_method = "wilcox.test", 
				verticals=T, lwd=2, pch=NA, xmin=-5, xmax=15, desc=desc, test_sample_pairs="first_oths")
		}
	}
}

dev.off()






##3, study relationship between motif and splicing or gene expression regulation
if(what2do %in% c("ano_motif_with_regu") ){
	#for motifs identified from DREME and scaned in all human genes, annotate them with overlapping 5'ss and associated regulation of the 5'ss
	#the output of this step will be used for cis elements analysis (address why many motifs are neither detected by RNA-seq nor show regulation in RNA-seq)
	motif_name="awGAGTArg"
	motif_name="motif_GAGTA"
	motif2gene_f=paste("../../motif_ana/Motif_search/12motif2gene/",study_name,"/",motif_name,".motif.tbl",sep="")
	motif2ss_f=paste("../../motif_ana/Motif_search/13motif2ref/",study_name,"/ss5_mp500.sense/site2ref.",motif_name,".map_detail.tbl",sep="")
	AS_reg_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl",sep="")
	gex_f=paste("8gexch_p/",study_name,"/refseqcds.DESeq2.format_tb/allGene.txt",sep="")
	geno="hg19"
	transcript_defi_f=paste("../../club/3add_gene_ano/",geno,".refflat.desc.nonRedundant",sep="")

	out_stats_f=paste(AS_reg_f,".",motif_name,".stats.tbl",sep="")
	motif_out_f=paste(AS_reg_f,".added2",motif_name,".tbl",sep="")
	selReguType_names=ReguType_names
	gchtp_names <- paste("GexType_",sample_pairs,sep=""); gchtp_names
	gl2gch_names <- paste("L2FC_",sample_pairs,sep="")
	gpval_names <- paste("adjSLog10P_",sample_pairs, sep="")

	transcript_defi_d=read.table(transcript_defi_f, header=T, sep="\t", stringsAsFactors=F, quote="")
	motif2gene_d=read.table(motif2gene_f, header=T, sep="\t", stringsAsFactors=F, quote="")
	motif2gene_d$motif_id=paste(motif2gene_d$chromosome, motif2gene_d$strand, motif2gene_d$motif_abs_start, motif2gene_d$motif_abs_end,sep=":")
	#motif2gene_d=motif2gene_d[!duplicated(motif2gene_d$motif_id), ]; nrow(motif2gene_d)
	transcript2motifnum=table(motif2gene_d$transcript)
	transcript_defi_d$geneMotif_num=0
	transcript_defi_d$geneMotif_num[match(names(transcript2motifnum), transcript_defi_d[,transc_id_header])]=transcript2motifnum
	quantile(transcript_defi_d$geneMotif_num)
	transcript_defi_d[transcript_defi_d$geneMotif_num==max(transcript_defi_d$geneMotif_num), ]
	#annotate with shortest motif to pA distance
	motif2gene_d$motif2pA_dist=(motif2gene_d$motif_abs_end-motif2gene_d$transc_end)*ifelse(motif2gene_d$strand=="-",-1,1); quantile(motif2gene_d$motif2pA_dist)
	motif2gene_d=motif2gene_d[order(motif2gene_d$motif2pA_dist,decreasing=T),]
	transcript_defi_d$minMotif_pADist=motif2gene_d$motif2pA_dist[match(transcript_defi_d[,transc_id_header], motif2gene_d$transcript)]
	transcript_defi_d$minMotif_pADist=abs(transcript_defi_d$minMotif_pADist)
	table(is.na(transcript_defi_d$minMotif_pADist))
	quantile(transcript_defi_d$minMotif_pADist,na.rm=T)
	#output transcript_defi_d
	out_gene_f=paste(motif2gene_f,".gene.txt",sep="")
	write.table(transcript_defi_d[setdiff(names(transcript_defi_d), c("exon_starts","exon_ends"))], file=out_gene_f, col.names=T, row.names=F, sep="\t", quote=F)



	motif2ss_d=read.table(motif2ss_f, header=T, sep="\t", stringsAsFactors=F, quote="")
	ifsel=motif2ss_d$rel_pos1<0 & motif2ss_d$rel_pos2>0; sum(ifsel)
	motif2ss_d=motif2ss_d[ifsel,]
	motif2ss_d$motif_id=paste(motif2ss_d$chromosome, motif2ss_d$strand, motif2ss_d$abs_pos1, motif2ss_d$abs_pos2,sep=":")
	table(motif2ss_d$motif_id %in% motif2gene_d$motif_id)

	AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="",comment.char="")
	table(AS_reg_d$region_ano)
	AS_reg_d=AS_reg_d[!is.na(AS_reg_d$strand),]
	#add 5'ss position
	AS_reg_d$ss5_pos=NA
	strand_s=ifelse(AS_reg_d$strand=="-",-1,1); table(strand_s)
	ifStartPos_ss5=AS_reg_d$startPosType=="ss5"
	AS_reg_d$ss5_pos[ifStartPos_ss5] = (AS_reg_d$start_pos - strand_s)  [ifStartPos_ss5]
	ifEndPos_ss5=AS_reg_d$endPosType=="ss5"; sum(ifEndPos_ss5)
	AS_reg_d$ss5_pos[ifEndPos_ss5] = AS_reg_d$end_pos [ifEndPos_ss5]
	table(is.na(AS_reg_d$ss5_pos), AS_reg_d$region_ano)
	AS_reg_d=AS_reg_d[!is.na(AS_reg_d$ss5_pos), ]
	AS_reg_d$ss5_id=paste(AS_reg_d$contig, AS_reg_d$strand, AS_reg_d$ss5_pos, sep=":")

	
	#add ss5 to motif2gene_d
	motif2gene_d$ss5_id=motif2ss_d$refSiteID[match(motif2gene_d$motif_id, motif2ss_d$motif_id)]; table(is.na(motif2gene_d$ss5_id))
	#add motif to AS_reg_d
	AS_reg_d$if_ss5_motif=ifelse(AS_reg_d$ss5_id %in% motif2ss_d$refSiteID,"Y","N"); table(AS_reg_d$if_ss5_motif)
	##output AS_reg_d (with if_ss5_motif information)
	write.table(AS_reg_d[setdiff(names(AS_reg_d),c("startSS_seq","endSS_seq"))], file=paste(AS_reg_f,".addSS5Motif.tbl",sep=""), col.names=T, row.names=F, sep="\t", quote=F)

	#study "5\'ss motif vs. AS regulation"
	write(paste("5\'ss motif vs. AS regulation",sep=""), file=out_stats_f, append=F)
	write(paste("motif2gene_f=",motif2gene_f,sep=""), file=out_stats_f, append=T)
	for(ReguType_name in ReguType_names){
		write(paste("\n",sep=""), file=out_stats_f, append=T)
		Gblock_anos2=intersect(Gblock_anos, AS_reg_d$region_ano)
		sum_tb=table( AS_reg_d$if_ss5_motif, AS_reg_d[,ReguType_name], AS_reg_d$region_ano ) [c("Y","N"),c("UP","NC","DN","na"),Gblock_anos2]
		SS_tb=data.frame(comparison=c("UP vs. NC", "UP vs. DN"))
		num_combine_tb=NULL
		for( region_ano in Gblock_anos2){
			sub_tb=sum_tb[,, region_ano]
			pval_UP_NC=do_fisher_test(sub_tb[,c("UP","NC")])
			pval_UP_DN=do_fisher_test(sub_tb[,c("UP","DN")])
			SS_tb[region_ano]= c(pval_UP_NC, pval_UP_DN)
			num_combine_tb=cbind(num_combine_tb, sub_tb[,c("UP","NC")])
			write(paste("ReguType_name=",ReguType_name, "; region_ano=",region_ano, "; Fisher\'s exact test P(-S*log10), UP/NC=", pval_UP_NC, " UP/DN=",pval_UP_DN,sep=""), file=out_stats_f, append=T)
			write.table(cbind(rownames(sub_tb),sub_tb), file=out_stats_f, append=T, col.names=T, row.names=F, sep="\t", quote=F)
		}
		num_combine_tb2=cbind(c("","","",rownames(num_combine_tb)) , 
			rbind(paste(ReguType_name,"vs. motif"), rep(Gblock_anos2,each=2), colnames(num_combine_tb), num_combine_tb))
		num_combine_tb2[1,-2]=""
		num_combine_tb2[2,c(T,F)]=""
		write(paste("Fisher\'s exact test P(-S*log10), S=1 for enrichment in with motif group and -1 otherwise:",sep=""), file=out_stats_f, append=T)
		write.table(num_combine_tb2, file=out_stats_f, append=T, col.names=F, row.names=F, sep="\t", quote=F)
		write.table(SS_tb, file=out_stats_f, append=T, col.names=T, row.names=F, sep="\t", quote=F)
	}


	#study "intron REMS vs. intron-retention" (optional)
	ifsel=AS_reg_d$region_ano %in% c("intron"); sum(ifsel)
	AS_reg_d1=AS_reg_d[ifsel, ] ; nrow(AS_reg_d1)
	AS_reg_d1=AS_reg_d1[ grep("^E",AS_reg_d1$startPosAno), ]; nrow(AS_reg_d1)
	AS_reg_d1=AS_reg_d1[ grep("^E",AS_reg_d1$endPosAno), ]; nrow(AS_reg_d1) #require whole intron are annotated (known)
	AS_reg_d1$upsExonID=paste(AS_reg_d1[, transc_id_header], AS_reg_d1$startPosAno,sep=":")
	motif2gene_d2=motif2gene_d[grep("^i", motif2gene_d$motif_PosAno), ]
	motif_PosAno_upsExon=sub("i","E",motif2gene_d2$motif_PosAno)
	Intron2MotifNum_v=table(paste(motif2gene_d2$transcript,motif_PosAno_upsExon,sep=":"))
	Intron2MotifNum_v=Intron2MotifNum_v[names(Intron2MotifNum_v) %in% AS_reg_d1$upsExonID]
	AS_reg_d1$IntronMotifNum=0
	AS_reg_d1$IntronMotifNum[match(names(Intron2MotifNum_v), as.character(AS_reg_d1$upsExonID)) ] = Intron2MotifNum_v
	IntronMotifNumBinTypes=c("0","1",">1")
	AS_reg_d1$IntronMotifNumBin=fromValue_break2Class(AS_reg_d1$IntronMotifNum, c(1,2), IntronMotifNumBinTypes ); table(AS_reg_d1$IntronMotifNumBin)

	AS_reg_d1$intron_size=abs(AS_reg_d1$start_pos-AS_reg_d1$end_pos)+1
	intron_size_binTypes=c("<300","mid",">5000")
	AS_reg_d1$intron_size_bin=fromValue_break2Class(AS_reg_d1$intron_size, c(300,5000), intron_size_binTypes ); table(AS_reg_d1$intron_size_bin) #roughly 20%, 60%, 20%
	#add log2ratio values
	in_cbd2=read.table(out_alltb_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	AS_reg_d1[,c(Log2Ratio_names,adjpval_names)]=in_cbd2[match( AS_reg_d1$Gblock_id,  in_cbd2$Gblock_id), c(Log2Ratio_names,adjpval_names)]
	##output AS_reg_d1 (with IntronMotifNum information)
	write.table(AS_reg_d1[setdiff(names(AS_reg_d1),c("startSS_seq","endSS_seq"))], file=paste(AS_reg_f,".IntronMotif.tbl",sep=""), col.names=T, row.names=F, sep="\t", quote=F)
	
	ifsel=AS_reg_d1$if_ss5_motif=="N"
	AS_reg_d2=AS_reg_d1[ifsel, ] #only introns with no motif in 5'ss
	compReguTypes=c("UP","NC")
	for(ReguType_name in selReguType_names){
		alltb=table(AS_reg_d2$IntronMotifNumBin, AS_reg_d2$intron_size_bin, AS_reg_d2[,ReguType_name]) [IntronMotifNumBinTypes,intron_size_binTypes,]
		numtb1=alltb[,, compReguTypes[1]]
		numtb2=alltb[,, compReguTypes[2]]
		perctb1=numtb1/sum(numtb1)
		perctb2=numtb2/sum(numtb2)
		delta_perctb=perctb1-perctb2
		
		col_labels=c("", rep(rep(compReguTypes,each=ncol(numtb1)),2),   rep(paste(compReguTypes[1],compReguTypes[2],sep="-"), ncol(numtb1)) )
		combineTb=cbind(IntronMotifNumBinTypes, numtb1,numtb2,  perctb1, perctb2,  delta_perctb)
		colnames(combineTb)[1]=""
		write(paste("matrix of IntronMotifNumBin (rows) x intron_size_bin (columns)","; ReguType_name=",ReguType_name, sep=""), file=out_stats_f, append=T)
		write.table(rbind(col_labels,colnames(combineTb),combineTb), file=out_stats_f, append=T, col.names=F, row.names=F, sep="\t", quote=F)
	}
	##draw box plot of log2ratio for introns with different length and motif number
	out_img_f=paste(AS_reg_f,".motif.intronlength.boxplot.pdf",sep="")
	pdf(out_img_f, width=18, height=12)
	layout(matrix(1:6,2,byrow=T))
	for(samplePair_id in match(selReguType_names, ReguType_names)){
		ReguType_name=ReguType_names[samplePair_id]
		log2ratioName=Log2Ratio_names[samplePair_id]
		pval_name=adjpval_names[samplePair_id]
		ifsel=AS_reg_d2[,ReguType_name] %in% c("UP","NC")
		ifsel=!is.na(AS_reg_d2[,pval_name]) & AS_reg_d2[,pval_name]<0.05
		AS_reg_d3=AS_reg_d2[ifsel, ]
		draw_2grped_boxplot(AS_reg_d3, "IntronMotifNumBin", IntronMotifNumBinTypes, "intron_size_bin", intron_size_binTypes[-1], plot_var_name=log2ratioName, 
			plot_var_snames=log2ratioName, mycolor=c("blue","cyan","red"), title=sample_pairs[samplePair_id], boxwex=0.4, ymin=-4, ymax=3 )
		
		draw_grped_boxplot(AS_reg_d3, "IntronMotifNumBin", IntronMotifNumBinTypes, log2ratioName, ymin=-4, ymax=3, spaceBwGrp=0,boxwex=1.3)
		draw_grped_boxplot(AS_reg_d3, "intron_size_bin", intron_size_binTypes, log2ratioName, ymin=-4, ymax=3, spaceBwGrp=0,boxwex=1.3)
	}
	dev.off()


	#add regulation and other annotations to motif2gene_d (output a file for cis elements analysis)
	update_headers=c("region_ano", selReguType_names)
	AS_reg_d1=AS_reg_d[order(match(AS_reg_d$region_ano,Gblock_anos)),] #AS_reg_d1 is redundant in terms of 5'ss, sorting the table based on region_ano allow certain priority  based on region_ano
	motif2gene_d[,update_headers]=AS_reg_d1[match(motif2gene_d$ss5_id, AS_reg_d1$ss5_id), update_headers]
	table(motif2gene_d[,selReguType_names[1]], motif2gene_d$region_ano)
	write.table(motif2gene_d, file=motif_out_f, col.names=T, row.names=F, sep="\t", quote=F)

	##study motif number and gene abundance change
	gex_d=read.table(gex_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F )
	gex_d$geneMotif_num= transcript_defi_d$geneMotif_num[match(gex_d$gene_symbol, transcript_defi_d$gene_symbol)]
	genesWith5ssMotif=unique(AS_reg_d$gene_symbol[AS_reg_d$if_ss5_motif=="Y"]); length(genesWith5ssMotif)
	gex_d$if_ss5_motif=ifelse(gex_d$gene_symbol %in% genesWith5ssMotif, "Y","N")
	table(gex_d$geneMotif_num>0, gex_d$if_ss5_motif) #double check 
	#bin geneMotif_num
	quantile(gex_d$geneMotif_num[gex_d$geneMotif_num>0], probs=0:10/10)
	GeneMotifNumBinTypes=c("0","1-3","4-29",">=30")
	geneMotifTypes=c(GeneMotifNumBinTypes,"ss5REMS","nSs5REMS","nSs5REMS.InFrame","nSs5REMS.FS","nSs5REMS.Stop")
	gex_d$geneMotifNumBin=fromValue_break2Class(gex_d$geneMotif_num, c(1,4,30), GeneMotifNumBinTypes ); table(gex_d$geneMotifNumBin)[GeneMotifNumBinTypes]
	gex_d$geneMotifType=as.character(gex_d$geneMotifNumBin)
	gex_d$geneMotifType[gex_d$if_ss5_motif=="Y"]="ss5REMS"; table(gex_d$geneMotifType)
	#anotate with novel exon ss5 motif and whether frame shift and stop codon in it
	if_novel_ss5Exon=AS_reg_d$region_ano=="exon" & grepl("^n",AS_reg_d$endPosAno); sum(if_novel_ss5Exon)
	if_frameShift=AS_reg_d$if_ss5_motif=="Y" & if_novel_ss5Exon & !is.na(AS_reg_d$exon_i_inFrame) & AS_reg_d$exon_i_inFrame=="N"; sum(if_frameShift)
	if_Inframe=AS_reg_d$if_ss5_motif=="Y" & if_novel_ss5Exon & !is.na(AS_reg_d$exon_i_inFrame) & AS_reg_d$exon_i_inFrame=="Y" & !is.na(AS_reg_d$exon_i_stopCodon) & AS_reg_d$exon_i_stopCodon=="N"; sum(if_Inframe)
	if_InframeStopCodon=AS_reg_d$if_ss5_motif=="Y" & if_novel_ss5Exon & !is.na(AS_reg_d$exon_i_inFrame) & AS_reg_d$exon_i_inFrame=="Y" & !is.na(AS_reg_d$exon_i_stopCodon) & AS_reg_d$exon_i_stopCodon=="Y"; sum(if_InframeStopCodon)
 	gex_d$geneMotifType[gex_d$if_ss5_motif=="Y" & gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_novel_ss5Exon]  ]="nSs5REMS";  table(gex_d$geneMotifType)
 	gex_d$geneMotifType[gex_d$if_ss5_motif=="Y" & gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_Inframe]  ]="nSs5REMS.InFrame";  table(gex_d$geneMotifType)
 	gex_d$geneMotifType[gex_d$if_ss5_motif=="Y" & gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_frameShift]  ]="nSs5REMS.FS";  table(gex_d$geneMotifType)
 	gex_d$geneMotifType[gex_d$if_ss5_motif=="Y" & gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_InframeStopCodon]  ]="nSs5REMS.Stop";  table(gex_d$geneMotifType)
 	table(gex_d$geneMotifType)[geneMotifTypes]
	#plotCols=c("blue","cyan","orange","red","purple","","","")
	#annotate with pA to motif distance (transcript_defi_d$minMotif_pADist)
	gex_d$minMotif_pADist=transcript_defi_d$minMotif_pADist[match(gex_d$gene_symbol, transcript_defi_d$gene_symbol)]
	quantile(gex_d$minMotif_pADist,na.rm=T)
	minMotif_pADist_bintypes=c("<0.5k:1nonss5REMS", "0.5-1k:1nonss5REMS", "1-2k:1nonss5REMS", ">=2k:1nonss5REMS","<0.5k","0.5-1k","1-2k",">=2k" )
	gex_d$minMotif_pADist_bin=fromValue_break2Class(gex_d$minMotif_pADist, c(500,1000,2000), minMotif_pADist_bintypes[5:8] ); 
	if_fewNonss5Motif=!is.na(gex_d$minMotif_pADist_bin) & gex_d$geneMotif_num ==1
	gex_d$minMotif_pADist_bin[if_fewNonss5Motif]=paste(gex_d$minMotif_pADist_bin[if_fewNonss5Motif], "1nonss5REMS",sep=":")
	table(gex_d$minMotif_pADist_bin)[minMotif_pADist_bintypes]
	
	#output gex_d
	out_gex_f=paste(gex_f,".addMotifInfo.txt",sep="")
	write.table(gex_d, file=out_gex_f, col.names=T, row.names=F, sep="\t", quote=F)

	plot_l=list(
		geneMotifNumBin=GeneMotifNumBinTypes,
		geneMotifType=geneMotifTypes,
		minMotif_pADist_bin=minMotif_pADist_bintypes
	)
	out_img_f=paste(AS_reg_f,".GeneMotif.Gex.boxplot.pdf",sep="")
	out_stats_f=paste(AS_reg_f,".GeneMotif.Gex.stats.txt",sep="")
	write("gene group numbers:", file=out_stats_f,)
	pdf(out_img_f, width=18, height=8); 
	for(grp_var in names(plot_l)){
		grps=plot_l[[grp_var]]
		tb=table(gex_d[,grp_var])[grps]
		write(paste(grp_var,"\tnum",sep=""), file=out_stats_f, append=T)
		write.table(cbind(names(tb),tb), file=out_stats_f, col.names=F, row.names=F, sep="\t", quote=F , append=T)
		plotCols=rainbow(length(grps))
		layout(matrix(1:2,1))
		for(gl2gch_name in gl2gch_names){
			draw_grped_boxplot(gex_d, grp_var, grps, gl2gch_name, ymin=-1.5, ymax=1, spaceBwGrp=0, 
				ylabel=gl2gch_name, xlabel="geneMotifNumBin"  , width=1.3, mycolor=plotCols[1:length(grps)],
				test_sample_pairs="first_oths", if_legend=gl2gch_name==gl2gch_names[length(gl2gch_names)])
		}
	}
	dev.off()

	out_img_f=paste(AS_reg_f,".GeneMotif.Gex.cdf.pdf",sep="")
	pdf(out_img_f, width=16, height=8); 
	for(grp_var in names(plot_l)){
		layout(matrix(1:2,1))
		grps=plot_l[[grp_var]]
		plotCols=rainbow(length(grps))
		draw_ecdf_density(gex_d, gl2gch_names, grp_var, grps, mycols=plotCols[1:length(grps)], 
			test_method = "wilcox.test", test_sample_pairs = "first_oths", verticals=T, pch=NA, lwd=2, xmin=-2, xmax=2, iflegend=T )

	}
	dev.off()

}


##4, study 5'ss -4 and -3 position di-nucleotide frequency of up-regulated exons
AS_reg_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl",sep="")
AS_reg_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.txt",sep="") #regulated only
AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
table(AS_reg_d$region_ano)
NTs=c("A","C","G","T")
AS_reg_d1=AS_reg_d[AS_reg_d$region_ano=="exon" & substr(AS_reg_d$endSS_seq,11,12)=="gt",]; nrow(AS_reg_d1)
AS_reg_d1$ss5m4=substr(AS_reg_d1$endSS_seq,7,7)
AS_reg_d1$ss5m3=substr(AS_reg_d1$endSS_seq,8,8)
AS_reg_d1$ss5m4m3=substr(AS_reg_d1$endSS_seq,7,8)
AS_reg_d1$ss5m2m1=substr(AS_reg_d1$endSS_seq,9,10)
diNTs=apply(expand.grid(NTs,NTs),1,paste,collapse=""); diNTs
sum_tb=NULL
selReguType_names=ReguType_names ###
for(ReguType_name in selReguType_names){
	for(ReguType in c("UP","NC") ){
			ifsel=AS_reg_d1[,ReguType_name]==ReguType
			AS_reg_d2=AS_reg_d1[ifsel,]
			if_GA=AS_reg_d2$ss5m2m1=="GA"
			tb1=table(factor(AS_reg_d2[ if_GA, "ss5m4"],levels=NTs), factor(AS_reg_d2[ if_GA, "ss5m3"],levels=NTs) ) [NTs,NTs]
			tb2=table(factor(AS_reg_d2[!if_GA, "ss5m4"],levels=NTs), factor(AS_reg_d2[!if_GA, "ss5m3"],levels=NTs) ) [NTs,NTs]
			chisq_res1=chisq.test(tb1)
			chisq_res2=chisq.test(tb2)
			m4Nums1=rowSums(tb1); m3Nums1=colSums(tb1)
			m4Nums2=rowSums(tb2); m3Nums2=colSums(tb2)
			
			sum_tb=rbind(sum_tb, c(m4Nums1,m4Nums1/sum(m4Nums1), m3Nums1, m3Nums1/sum(m3Nums1), 
				tb1, sum(tb1),chisq_res1$p.value,  (tb1-chisq_res1$expected)/sum(tb1), tb1/sum(tb1), chisq_res1$expected  ) )
			sum_tb=rbind(sum_tb, c(m4Nums2,m4Nums2/sum(m4Nums2), m3Nums2, m3Nums2/sum(m3Nums2), 
				tb2, sum(tb2), chisq_res2$p.value, (tb2-chisq_res2$expected)/sum(tb2), tb2/sum(tb2), chisq_res2$expected ) )

	}
}
colnames(sum_tb)=c(paste("m4",NTs,sep="."),paste("m4",NTs,"perc",sep="."), paste("m3",NTs,sep="."),paste("m3",NTs,"perc",sep="."),
	 diNTs,"total","chisq.P", paste(diNTs,"obs-exp.perc",sep="."), paste(diNTs,"perc",sep="."), paste(diNTs,"expected",sep=".") )
rownames(sum_tb)=paste(rep(selReguType_names,each=4), c("UP.GA","UP.non-GA","NC.GA","NC.non-GA"), sep="." )

out_stats_f=paste(AS_reg_f,".ss5m4m3diNT.stats.tbl",sep="")
write.table(cbind(rownames(sum_tb),sum_tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F)



##5, study 5'ss derived from transposable elements
ss5_TE_ano_f=paste("../../AS/ss_annotation/03anoSS_TE/",study_name,"/",geno,".TE.to.juncpos5.tbl", sep="")
repClasses=c("LINE","SINE","LTR","DNA")
ss5_TE_ano_d=read.table(ss5_TE_ano_f, header=T, sep="\t", quote="", stringsAsFactors=F)
ss5_TE_ano_d=ss5_TE_ano_d[order(ss5_TE_ano_d$TE_fullLen,decreasing=T),]
table(duplicated(ss5_TE_ano_d$ss_id))
TEinfo_headers=names(ss5_TE_ano_d)[-1]
ss5_TE_ano_d=ss5_TE_ano_d[ss5_TE_ano_d$repClass %in% repClasses,]

AS_reg_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl",sep="")
AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
table(AS_reg_d$region_ano)
AS_reg_d1=AS_reg_d[AS_reg_d$region_ano=="exon" & substr(AS_reg_d$endSS_seq,11,12)=="gt",]; nrow(AS_reg_d1)

AS_reg_d1$ss5_id=paste(AS_reg_d1$contig, AS_reg_d1$strand, AS_reg_d1$end_pos,sep=":")
table(AS_reg_d1$ss5_id %in% ss5_TE_ano_d$ss_id)
AS_reg_d1$ss5m2m1=substr(AS_reg_d1$endSS_seq,9,10)
AS_reg_d1[,TEinfo_headers]=ss5_TE_ano_d[match(AS_reg_d1$ss5_id, ss5_TE_ano_d$ss_id),TEinfo_headers]

repClass_l=tapply(ss5_TE_ano_d$repFamily, ss5_TE_ano_d$repClass, function(v){
	nums=sort(table(v),decreasing=T)
	setdiff( names(nums), grep("\\?",names(nums),value=T) )
} )[repClasses]
repClass_l

studyReguTypes=c("UP","NC")
exonGrps=c("GA;NC","GA;UP","nonGA;NC","nonGA;UP")

selReguType_names=ReguType_names[c(2,4,6,8,1,3,5,7,9)] ###
selReguType_names=ReguType_names[1] ###
cbsum_tb=data.frame()
for(ReguType_name in selReguType_names){
	print(ReguType_name)
	sum_tb=data.frame(ReguType_name=rep(ReguType_name,length(exonGrps)), exonGroup=exonGrps)
	ifsel=AS_reg_d1[,ReguType_name] %in% studyReguTypes
	AS_reg_d2=AS_reg_d1[ifsel,]
	exonGrp=paste(ifelse(AS_reg_d2$ss5m2m1=="GA","GA","nonGA"), AS_reg_d2[,ReguType_name], sep=";" )
	tb=table(exonGrp, !is.na(AS_reg_d2$repClass))[exonGrps,c("TRUE","FALSE")]
	tb[is.na(tb)]=0
	perc=tb[,1] / rowSums(tb)
	stderr= (perc*(1-perc)/rowSums(tb))^0.5
	sum_tb$TotalNum=rowSums(tb)
	sum_tb$num_TE=tb[,1]
	sum_tb$perc_TE=perc
	sum_tb$perc_TE.stderr=stderr
	sum_tb[paste("num",names(repClass_l),sep=".")]=table(exonGrp,AS_reg_d2$repClass)[exonGrps, names(repClass_l)]
	tb2=table(exonGrp,factor(AS_reg_d2$repFamily,levels=unlist(repClass_l)), AS_reg_d2$repClass)[exonGrps, ,names(repClass_l)]
	for(repClass in names(repClass_l) ){
		sum_tb[paste("num",repClass,repClass_l[[repClass]],sep=".")]=tb2[, repClass_l[[repClass]] , repClass]
	}
	cbsum_tb=rbind(cbsum_tb,sum_tb)
}
out_stats_f=paste(AS_reg_f,".ss5_matchTE.stats.tbl",sep="")
write.table(cbsum_tb, file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F)

##output a table of TE-derived exons (5'ss)
if_out=!is.na(AS_reg_d1$repClass)
out_data_f=paste(AS_reg_f,".ss5_matchTE.data.tbl",sep="")
write.table(AS_reg_d1[if_out,], file=out_data_f, sep="\t", quote=F, col.names=T, row.names=F)




##6.1, study splice site annotation vs. regulation (cryptic exon vs. annotated exon)
AS_reg_f=paste("27ana_juncreads/",study_name,"/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl",sep="")
AS_reg_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.allReguType.txt",sep="")
AS_reg_f=paste("../../dep_seq/othProject/",study_name,"/RNAseq_Splicing/06.exon_PSI/all_exon/formatted_tb/exons.delta_PSI20_Pfisher0.001.allReguType.txt",sep="")

AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")

table(AS_reg_d$region_ano)
out_stats_f=paste(AS_reg_f,".regulation_vs._annotation.stats.tbl",sep="")

ssSupp_headers=c("startSS_supp","endSS_supp")
ssSupp_types=c("Refseq","KnG_Ens","mRNA_EST","None")
#combine start and end position annotation to one
AS_reg_d$region_annotation="None"
if_start_annotated=AS_reg_d$startSS_supp %in% c(ssSupp_types[1:2]); table(if_start_annotated)
if_end_annotated=AS_reg_d$endSS_supp %in% c(ssSupp_types[1:2]); table(if_end_annotated)
AS_reg_d$region_annotation[if_start_annotated | if_end_annotated]="Either"
AS_reg_d$region_annotation[if_start_annotated & if_end_annotated]="Both"
table(AS_reg_d$region_annotation)
region_annotation_types=c("Both","Either","None")

write(paste(c("gene block regulation vs. annotation:\nAnnotated includes splice site in: ", ssSupp_types[1:2]),collapse=" "),  file=out_stats_f, append=F)
write(paste("Splicing regulation file=",AS_reg_f,collapse=" "),  file=out_stats_f, append=T)

region_anos=unique(AS_reg_d$region_ano) ; region_anos=region_anos[!is.na(region_anos)]
for(region_ano in region_anos){
	print(region_ano)
	ifsel=!is.na(AS_reg_d$region_ano) & AS_reg_d$region_ano==region_ano
	sum_tb=NULL
	headers1=NULL; headers2=NULL; headers3=NULL
	for(ssAno_header in "region_annotation"){
		for(ReguType_name in ReguType_names){
			tb=table( factor(AS_reg_d[ifsel, ssAno_header],levels=region_annotation_types),  factor(AS_reg_d[ifsel, ReguType_name],levels=ReguTypes))[region_annotation_types,ReguTypes]
			headers3=c(headers3,colnames(tb))
			headers2=c(headers2,c(ReguType_name, rep("", ncol(tb)-1) ) )
			sum_tb=cbind(sum_tb,tb)
		}
		headers1=c(headers1, c(ssAno_header, rep("", length(ReguType_names)*length(ReguTypes)-1 )  ) )
	}
	sum_tb2=rbind(headers1, headers2, headers3, sum_tb)
	sum_tb2=cbind(c("","","", region_annotation_types),sum_tb2)
	write(paste("\nRegion type=",region_ano),  file=out_stats_f, append=T)
	write.table(sum_tb2, file=out_stats_f, col.names=F, row.names=F, sep="\t", quote=F, append=T)
}


##6.2 study genes with cryptic exon inclusion and gene expression changes
table(AS_reg_d$region_annotation,AS_reg_d$region_ano)
out_AS_data_f=paste(AS_reg_f,".Exon.ssAnotation.txt",sep="")
dim(AS_reg_d)
remove_headers=NULL ; #c("startSS_seq","endSS_seq")
write.table(AS_reg_d[setdiff(names(AS_reg_d1), remove_headers)], file=out_AS_data_f, col.names=T, row.names=F, sep="\t", quote=F)


#open gene expression file
gex_f=paste("../../dep_seq/2gene_readsnum/8gexch_p/",study_name,"/refseqcds.DESeq2.format_tb/allGene.txt",sep="")
selReguType_names=ReguType_names
gchtp_names <- paste("GexType_",sample_pairs,sep=""); gchtp_names
gl2gch_names <- paste("L2FC_",sample_pairs,sep="")
gpval_names <- paste("adjSLog10P_",sample_pairs, sep="")

gex_d=read.table(gex_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F )
if_cryptic_exon=AS_reg_d$region_ano=="exon" & AS_reg_d$region_annotation %in% c("None","Either"); sum(if_cryptic_exon)
if_inframeNoPTC= if_cryptic_exon & !is.na(AS_reg_d$exon_i_inFrame) & AS_reg_d$exon_i_inFrame=="Y" & !is.na(AS_reg_d$exon_i_stopCodon) & AS_reg_d$exon_i_stopCodon=="N"; 
if_FS= if_cryptic_exon & !is.na(AS_reg_d$exon_i_inFrame) & AS_reg_d$exon_i_inFrame=="N"; 
table(if_cryptic_exon, if_FS)
if_PTC= if_cryptic_exon & !is.na(AS_reg_d$exon_i_stopCodon) & AS_reg_d$exon_i_stopCodon=="Y"; 
table(if_FS, if_PTC)

out_image_f=paste(AS_reg_f,".crypticExon_vs_Gex.pdf",sep="")
out_gex_f=paste(AS_reg_f,".crypticExon_vs_Gex.txt",sep="")
gene_cryptic_exon_types=c("Others","Inframe","FS","PTC","FS&PTC")
gene_cryptic_exon_types=c("Others","Inframe","FS|PTC")
gene_cryptic_exon_types=c("Others","FS|PTC")
pdf(out_image_f, width=12, height=8)
layout(matrix(1:6,2,byrow=F))
for(i in 1:length(ReguType_names)){
	ReguType_name=ReguType_names[i]
	gl2gch_name=gl2gch_names[i]
	if_exon_up=if_cryptic_exon & AS_reg_d[,ReguType_name]=="UP" & !is.na(AS_reg_d[,ReguType_name])
	sum(if_exon_up & if_FS)
	sum(if_exon_up & if_PTC)
	sum(if_exon_up & if_FS & if_PTC)
	gene_cryptic_exon_type_name=paste("activated_iExon_type.",sample_pairs[i],sep="")
	gex_d[,gene_cryptic_exon_type_name]="Others"
	gex_d[,gene_cryptic_exon_type_name][gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_exon_up & if_cryptic_exon & if_inframeNoPTC]]="Inframe"
	gex_d[,gene_cryptic_exon_type_name][gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_exon_up & if_FS]]="FS"
	gex_d[,gene_cryptic_exon_type_name][gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_exon_up & if_PTC]]="PTC"
	gex_d[,gene_cryptic_exon_type_name][gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_exon_up & if_FS & if_PTC]]="FS&PTC"
	gex_d[,gene_cryptic_exon_type_name][gex_d$gene_symbol %in% AS_reg_d$gene_symbol[if_exon_up & (if_FS | if_PTC)]]="FS|PTC"
	table(gex_d[,gene_cryptic_exon_type_name])[gene_cryptic_exon_types]
	draw_ecdf_density(gex_d,gl2gch_name,gene_cryptic_exon_type_name,gene_cryptic_exon_types,mycols=c("black","blue","red","green","orange"),
		desc="", xmin = -2, xmax = 1, test_method = "wilcox.test", test_sample_pairs = "first_oths", verticals=T, lwd=2, pch=19, cex=1, do.points=T)
	abline(v=0, lty=2, col="gray", lwd=2)
}
dev.off()

#output gene table
write.table(gex_d,file=out_gex_f, col.names=T, row.names=F, sep="\t", quote=F)



##6.3, study basal PSI for different groups of exons (eg. iExon up regulated vs. annotated exon NC ... ...)
AS_reg_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.allReguType.txt",sep="") #all exons with regulation type
all_gblock_f=paste("../../AS/5AS_deepSeq/04exon_PSI/",study_name,"/all_exon/exons.tbl",sep="") #used for additional annotation

AS_reg_f=paste("../../dep_seq/othProject/",study_name,"/RNAseq_Splicing/06.exon_PSI/all_exon/formatted_tb/exons.delta_PSI20_Pfisher0.001.allReguType.txt",sep="")
all_gblock_f=paste("../../dep_seq/othProject/",study_name,"/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl",sep="") #used for additional annotation

AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
basal_PSI_names="PSI_Ctrl_C3" # PSI_Ctrl_C3  PSI_DMSO
required_headers=c("startSS_seq","endSS_seq",basal_PSI_names)
if(!all(required_headers %in% names(AS_reg_d))){
	all_gblock_d=read.table(all_gblock_f, header=T, sep="\t", stringsAsFactors=F, quote="", comment.char="")
	update_headers=setdiff(required_headers, names(AS_reg_d))
	all_gblock_d$Gblock_id2=paste(all_gblock_d$contig, all_gblock_d$strand, all_gblock_d$start_pos, all_gblock_d$end_pos, sep=":")
	table(AS_reg_d$Gblock_id2 %in% all_gblock_d$Gblock_id2)
	AS_reg_d[update_headers]=all_gblock_d[match(AS_reg_d$Gblock_id2 ,all_gblock_d$Gblock_id2), update_headers]
}
table(AS_reg_d$region_ano)

#go to ##6.1, study splice site annotation vs. regulation (annotate exons with column region_annotation) .....
ifsel=AS_reg_d$region_ano=="exon" & substr(AS_reg_d$endSS_seq, 11, 12)=="gt" & substr(AS_reg_d$startSS_seq, 49, 50)=="ag" & AS_reg_d$reg_len<200; sum(ifsel)
AS_reg_d1=AS_reg_d[ifsel,]
sel_ReguType_names=ReguType_names[1]; sel_ReguType_names #ReguType_names[1:2] ReguType_names[1]
if_up_in_any=rowSums(!is.na(AS_reg_d1[sel_ReguType_names]) & AS_reg_d1[sel_ReguType_names]=="UP")>0; sum(if_up_in_any)
if_NC_in_all=apply(!is.na(AS_reg_d1[sel_ReguType_names]) & AS_reg_d1[sel_ReguType_names]=="NC", 1, all); sum(if_NC_in_all)
comb_regu_types=rep(NA, nrow(AS_reg_d1))
comb_regu_types[if_up_in_any]="UP"
comb_regu_types[if_NC_in_all]="NC"
table(comb_regu_types)
##define up-regulated in any, combine with region_annotation
AS_reg_d1$exon_type=paste(AS_reg_d1$region_annotation,comb_regu_types, sep=":"); table(AS_reg_d1$exon_type)
study_exon_types=c("Both:UP","Both:NC", "None:UP" )
study_exon_colors=c("blue","black", "red" ); ltys=1
study_exon_colors=c("black","gray", "black" ); ltys=c(3,2,1) #used for patent 
# Both:NA   Both:NC   Both:UP Either:NA Either:NC Either:UP   None:NA   None:NC      None:UP
#      1420     53055       203       104      5484        45       130      6467       205

out_img_f=paste(AS_reg_f,".Exon.basalPSI.density.pdf",sep="")
pdf(out_img_f, width=13, height=6)
layout(matrix(1:2, nrow=1))
for (method in c("density","ecdf")){
	for(basal_PSI_name in basal_PSI_names){
		median_vals=tapply(AS_reg_d1[,basal_PSI_name], AS_reg_d1[,"exon_type"], median, na.rm=T)[study_exon_types]
		#draw median value lines
		draw_ecdf_density(AS_reg_d1,basal_PSI_names, "exon_type", study_exon_types, mycols=study_exon_colors,
			method = method,  test_method = "wilcox.test", test_sample_pairs = "first_oths", iflegend=F, verticals=T, pch=NA, ltys=ltys, lwd=3,
			desc=paste(c("median=",median_vals),collapse=" "))
		abline(v=median_vals, col=study_exon_colors, lty=ltys, lwd=1.5)
	}
}
dev.off()

#output AS_reg_d1 (can be used for splice sites conservation score analysis)
out_AS_data_f=paste(AS_reg_f,".Exon.ssAnotation_Regulation.txt",sep="")
dim(AS_reg_d1)
remove_headers=c("startSS_seq","endSS_seq")
write.table(AS_reg_d1[setdiff(names(AS_reg_d1), remove_headers)], file=out_AS_data_f, col.names=T, row.names=F, sep="\t", quote=F)


##6.4 identify a set of expressed NMD-exons (exons that the inclusion of which will cause NMD of a gene)
load("/drive2/wli/analyze/dep_seq/othProject/CHDI_201908_iPSC_neuron_Cpds/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl.img")
#will load exon_defi_d
basal_PSI_cut=20

table(exon_defi_d$region_ano)
table(exon_defi_d$exon_i_inFrame)
table(exon_defi_d$exon_i_stopCodon)
sel_PSI_names=grep("DMSO",PSI_names,value=T)
avg_PSI_name="PSI_DMSO"
if_NMD_exons=exon_defi_d$region_ano %in% "exon" & (exon_defi_d$exon_i_inFrame %in% "N" | exon_defi_d$exon_i_stopCodon %in% "Y") & 
	( grepl("n-d",exon_defi_d$startPosAno) & grepl("n-d",exon_defi_d$endPosAno) ) ; table(if_NMD_exons); 
if_expressed=!is.na(exon_defi_d[,avg_PSI_name]) & exon_defi_d[,avg_PSI_name]>=basal_PSI_cut; table(if_expressed)
table(if_NMD_exons, if_expressed)

keep_headers=c("gene_symbol","coordinates","gene_id",transcript_id_name, "contig","strand","end_pos","start_pos", "reg_len","startPosAno", "endPosAno", 
	"startSS_seq", "endSS_seq", "exon_i_inFrame","exon_i_stopCodon","startSS_supp","endSS_supp", avg_PSI_name, sel_PSI_names,"gene_desc" )
NMD_exon_d=exon_defi_d[if_NMD_exons & if_expressed, keep_headers]; dim(NMD_exon_d)
out_NMD_exon_f=paste(out_root,"/",study_exon_type,"/formatted_tb/NMD_exons.basal_PSI",basal_PSI_cut,".tbl", sep="")
write.table(NMD_exon_d, file=out_NMD_exon_f, col.names=T, row.names=F, sep="\t", quote=F)



##7, output a table of junction read counts for GEO submission
setwd('/drive2/wli/analyze/Pipeline/RNAseq_Splicing')
load('/drive2/wli/analyze/dep_seq/othProject/2016_C3_828/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl.img')
keep_sample_ids=c(1:8); change_header_froms=c("intron_id", "PTC905"); change_header_tos=c("Junction_id", "HTT-C2") #2016_C3_828

setwd("/drive2/wli/analyze/dep_seq/2gene_readsnum")
load('/drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/201803_U1vars19/all_exon/exons.tbl.img')
keep_sample_ids=c(1,2,15,16); change_header_froms=c("intron_id", "v7CGA"); change_header_tos=c("Junction_id", "U1v-GA") #201803_U1vars19
out_root="/drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/201803_U1vars19/"

suppNum_names[keep_sample_ids]

read_num_f
read_num_d=read.table(read_num_f, header=T, sep="\t", stringsAsFactors=F, quote="")
read_num_d$intron_id=paste(read_num_d$contig, read_num_d$strand, read_num_d$juncpos5, read_num_d$juncpos3, sep=":")
read_num_d=read_num_d[!duplicated(read_num_d$intron_id),]


#remove all 0 counts rows
if_keep_row=rowSums(read_num_d[suppNum_names[keep_sample_ids]])>0; table(if_keep_row); table(is.na(if_keep_row))
keepheaders=c("gene_symbol", "intron_id", suppNum_names[keep_sample_ids]) ;  keepheaders
out_d=read_num_d[if_keep_row,keepheaders]
out_d=re_format_tb(out_d, ch_header_from = change_header_froms, ch_header_to = change_header_tos)
out_d=out_d[order(out_d[,1],out_d[,2]), ]
out_f_name=paste(out_root,"/",study_exon_type,"/formatted_tb/Junction.ReadCounts.GEO",".txt", sep=""); out_f_name
write.table(out_d,file=out_f_name, col.names=T, row.names=F, sep="\t", quote=F)
system(paste("gzip", out_f_name))



##8.1, study FDR based on DMSO samples in HT-RNAseq, find the best cutoff for significant events
study_name="HT_RNAseq_P2" #HT_RNAseq_P1 HT_RNAseq_P2
load(paste0("/cycvol1/NetApp_Bio_Data/ProcessedData/RNAseq/HT-RNAseq/",study_name,"/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl.img"))
Pval_setting=paste(paste(use_p_type,P_cut,sep=""), collapse="_")
AS_reg_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.delta_PSI",delta_PSI_cut,"_",Pval_setting, ".txt", sep="")
out_FDR_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.FDR.txt", sep="")

#AS_reg_d=read.table(AS_reg_f, header=T, sep="\t", stringsAsFactors=F, quote="")
AS_reg_d=exon_defi_d
AS_reg_d=cal_PSI_regu_type(AS_reg_d, use_p_type, delta_PSI_names, c(0.01,0.05), 2, ReguType_names, sample_pairs, zeroDeltaPSI_Pval=zeroDeltaPSI_Pval)
if_regulated=rowSums(AS_reg_d[ReguType_names]=="UP" | AS_reg_d[ReguType_names]=="DN")>0; table(if_regulated)
AS_reg_d=AS_reg_d[if_regulated, ]; nrow(AS_reg_d)

DMSO_comps=grepl("DP",delta_PSI_names); sum(DMSO_comps) #any DMSO well compared to DMSO control group
DMSOgrp_comps=sample_compare_matrix[1,] %in% sample_repl_l$DMSO; DMSOgrp_comps; sum(DMSOgrp_comps) #DMSO wells also used in DMSO control group

##8.1, calculate FDR in different combination of cutoff values
#-use_p_type "Pfisher Pttest" -P_cut "1e-7 0.05" -PSI_cal_rnumMin 20 -delta_PSI_cut 20
fisher_pval_names; ttest_pval_names; log2R_names; ReguType_names
cutoff_matrix=expand.grid(Pttest=c(0.05,0.01, 1e-3,1e-4,1e-5,1e-6),Pfisher=c(1e-7, 1e-10, 1e-20, 1e-30), dPSI=c(20,25,30) )
use_p_type=names(cutoff_matrix)[1:2]
cutoff_settings=apply(cutoff_matrix,1,function(vals){paste(paste0(names(cutoff_matrix),vals), collapse=".")})
NumSig_names=paste0("NumSig.",cutoff_settings) #number of significant events given a cutoff
FDR_names=paste0("FDR.",cutoff_settings) #FDR for a cutoff for a treatment (well)
FDR_tb=data.frame(test_sample=sample_compare_matrix[1,], if_DMSO=DMSO_comps, if_inControl_group=DMSOgrp_comps)

define_parallel_fun(nCores=workers)
define_parallel_fun(nCores=8)

FDR_tb[NumSig_names]=unlist(myApply(1:nrow(cutoff_matrix), function(cutoff_i){
	cutoff_setting_i=cutoff_settings[cutoff_i]
	AS_reg_d2=cal_PSI_regu_type(AS_reg_d, use_p_type, delta_PSI_names, unlist(cutoff_matrix[cutoff_i, use_p_type]), unlist(cutoff_matrix[cutoff_i, "dPSI"]), ReguType_names, sample_pairs, zeroDeltaPSI_Pval=zeroDeltaPSI_Pval)
	Regu_num_tb=apply(AS_reg_d2[ReguType_names],2,function(v){table(factor(v, levels=ReguTypes) )[ReguTypes] })
	colSums(Regu_num_tb[1:2,])
}))
Avg_FD=colMeans(FDR_tb[DMSO_comps & !DMSOgrp_comps , NumSig_names]) #average false discovery (number) given a cutoff
FDR_tb[FDR_names]=t(Avg_FD/t(FDR_tb[NumSig_names]))
FDR_tb[FDR_names][FDR_tb[FDR_names]>1]=1
FDR_tb[DMSO_comps,FDR_names]=1
FDR_tb$Avg_FDR=rowMeans(FDR_tb[FDR_names])
FDR_tb_sorted=FDR_tb[order(DMSO_comps, DMSOgrp_comps,FDR_tb$Avg_FDR, decreasing=T),]
write.table(FDR_tb_sorted, file=out_FDR_stats_f, quote=F, col.names=T, row.names=F, sep="\t")

##8.2, develop code to calculate P-value adjusted delta-PSI
cal_Padj_change<-function(x, Pvals, b, e){
	adj_indx=1-1/(1+10^(b*(log10(Pvals)-log10(e))))
	return(prod(c(x,adj_indx)))
}

cal_Padj_change(20,c(1e-7,1e-30), c(-1,-1), c(1e-7,1e-30) ) # =20*0.5*0.5=5
cal_Padj_change(20,c(1e-7,1e-30), c(-1,-1), c(1e-4,1e-20) ) # close to 20
cal_Padj_change(20,c(1,1), c(-1,-1), c(1e-3,1e-15) ) # close to 0

#need AS_reg_d, fisher_pval_names; ttest_pval_names; log2R_names; ReguType_names
adjDeltaPSI_names=sub("deltaPSI_","AdjDeltaPSI_",delta_PSI_names); adjDeltaPSI_names
cal_avg_fdr_using_adjDelta_PSI<-function(parameters=c(-0.5,-0.2,1e-4, 1e-20, 20),  AS_reg_d=AS_reg_d){
	#eg. b=c(-1,-1); e=c(1e-4, 1e-20); delta_PSI_cutoff=0.2
	b=parameters[1:2]; e=parameters[3:4]; delta_PSI_cutoff=parameters[5]
	AS_reg_d2=AS_reg_d
	AS_reg_d2[,adjDeltaPSI_names] = AS_reg_d2[,delta_PSI_names]
	AS_reg_d2[,adjDeltaPSI_names]=unlist(myApply(1:length(adjDeltaPSI_names), function(i){
		apply(AS_reg_d2[,c(delta_PSI_names[i],ttest_pval_names[i], fisher_pval_names[i])],1,function(v){cal_Padj_change(v[1], v[2:3], b, e)})
	}))
	AS_reg_d2=cal_PSI_regu_type( AS_reg_d2, use_p_type, adjDeltaPSI_names, c(1,1), delta_PSI_cutoff, ReguType_names, sample_pairs, zeroDeltaPSI_Pval=c(1,1) )
	Regu_num_tb=apply(AS_reg_d2[ReguType_names],2,function(v){table(factor(v, levels=ReguTypes) )[ReguTypes] })
	sigNums=colSums(Regu_num_tb[1:2,])+1 #1 is a pseudo-count to make sure FDR is not 0
	Avg_FD=mean(sigNums[DMSO_comps & !DMSOgrp_comps]) #average false discovery
	Avg_D=mean(sigNums[!DMSO_comps]) #average discovery
	FDRs=Avg_FD/sigNums
	FDRs[FDRs>1]=1
	FDRs[DMSO_comps]=1
	mean_FDR=mean(FDRs[!DMSO_comps],na.rm=T)
	total_FDR=(sum(sigNums[DMSO_comps & !DMSOgrp_comps])+1) / (sum(sigNums[!DMSO_comps])+1) #1 is a pseudo count
	print(c(parameters, Avg_FD, Avg_D, mean_FDR, total_FDR))
	return(total_FDR)
}
optim_res=optim(c(-0.5,-0.2,1e-4, 1e-20, 20), cal_avg_fdr_using_adjDelta_PSI, NULL, method = "L-BFGS-B", lower = c(-10,-10 ,1e-5, 1e-20, 1), upper = c(-0.01,-0.01, 0.05, 0.01,10),  AS_reg_d=AS_reg_d )
##b=c(-0.5,-0.2) #for both plate 1 and 2

para_matrix=expand.grid(b_pttest=-0.5, b_Pfisher=-0.2, e_pttest=c(0.05, 0.01, 1e-3,  1e-4,  1e-5), e_Pfisher=c(1e-2, 1e-3,1e-4, 1e-5, 1e-10, 1e-15), delta_PSI_cutoff=c(2, 5, 10, 15, 20) ) #, 10,15,20
FDR_d=transform(para_matrix)
FDR_d$avg_FDR=NA
for(i in 1:nrow(para_matrix)){
	FDR_d$avg_FDR[i]=cal_avg_fdr_using_adjDelta_PSI(unlist(para_matrix[i,]), AS_reg_d=AS_reg_d)
}
out_FDR_stats_f=paste(out_root,"/",study_exon_type,"/formatted_tb/exons.totalFDR_optimi_adjDeltaPSI.txt", sep="")
write.table(FDR_d, file=out_FDR_stats_f, quote=F, col.names=T, row.names=F, sep="\t")

