##3/21/2016 can combine PSI data

source("../SharedCodes/Rfunc.inc.R") 
source("05.inc.func.R")


study_set="CRC_shSF1"
input_root="/HPCTMP_NOBKUP/wl314/analyze/Projects/TJ/1804_Colorectal_cancer_splicing/"
output_root="/HPCTMP_NOBKUP/wl314/analyze/Projects/TJ/1804_Colorectal_cancer_splicing/Splicing_comb_tb/"
ana_meth="DEXSeq"; fext="/RNAseq_Splicing/05.DEXSeq/geneBlock/combine_d.tbl"
ana_meth="annoExonPSI"; fext="/RNAseq_Splicing/06.exon_PSI/annotated_exon/exons.tbl"
ana_meth="AllExonPSI"; fext="/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl"
ana_meth="junc5PSI"; fext="/RNAseq_Splicing/05.junc.Ttest/junc5/combine_d.tbl"
ana_meth="junc3PSI"; fext="/RNAseq_Splicing/05.junc.Ttest/junc3/combine_d.tbl"

file_setting_list=list(
 list(study="GSE50760_2014",species="human",input_f=paste(input_root, "GSE50760_2014",fext, sep=""), change_headerFroms=c("NC","CRC","MC"),change_headerTos=c("NC2014","CRC2014","MC2014"),  test_samples=unlist(strsplit("CRC2014 MC2014"," ")), ref_samples=unlist(strsplit("NC2014 NC2014"," ") ))
 ,list(study="GSE95132_CRC",species="human",input_f=paste(input_root, "GSE95132_CRC", fext, sep=""), change_headerFroms=c("NC","CRC"),change_headerTos=c("NC132","CRC132"),  test_samples=unlist(strsplit("CRC132"," ")), ref_samples=unlist(strsplit("NC132"," ") ))
 ,list(study="GSE104178_CRC",species="human",input_f=paste(input_root, "GSE104178_CRC",fext, sep=""), change_headerFroms=c("NM","CRC"),change_headerTos=c("NM178","CRC178"),  test_samples=unlist(strsplit("CRC178"," ")), ref_samples=unlist(strsplit("NM178"," ") ))
 #,list(study="Encode_shSF1",species="human",input_f=paste(input_root, "Encode_shSF1",fext, sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("shSF1.a.K562 shSF1.b.K562 shSF1.a.HepG2"," ")), ref_samples=unlist(strsplit("shC.a.K562 shC.b.K562 shC.a.HepG2"," ") ))
)



all_test_samples=unlist( sapply(file_setting_list, function(file_setting1){file_setting1$test_samples}) )
all_ref_samples=unlist( sapply(file_setting_list, function(file_setting1){file_setting1$ref_samples}) )
sample_compare_matrix=rbind(all_test_samples,all_ref_samples)
sample_pairs=paste(all_test_samples,all_ref_samples,sep="_")
adjpval_names=paste("Padj_", sample_pairs,sep="")
pval_names=paste("P_", sample_pairs,sep="")
ReguType_names=paste("ReguType_", sample_pairs,sep="")
ReguTypes=c("UP","DN","NC","na")
Log2Ratio_names=paste("Log2R_", sample_pairs,sep="")
delta_PSI_names=paste("deltaPSI_", sample_pairs,sep="")
fisher_pval_names=paste("Pfisher_", sample_pairs,sep="")
ttest_pval_names=paste("Pttest_", sample_pairs,sep="")
control_PSI_names= paste("PSI_",unique(all_ref_samples),sep="")

use_P_names=pval_names; use_change_names=Log2Ratio_names; start_head="start_pos"; end_header="end_pos" #for gene block table
use_P_names=fisher_pval_names; use_change_names=delta_PSI_names; start_head="start_pos"; end_header="end_pos" #for all exons PSI table using Fisher's exact test
use_P_names=fisher_pval_names; use_change_names=delta_PSI_names; start_head="ups_ss3"; end_header="ss5" #for annotated exon PSI table using Fisher's exact test
use_P_names=ttest_pval_names; use_change_names=delta_PSI_names; start_head="ups_ss3"; end_header="ss5" #for annotated exon PSI table using ttest
use_P_names=pval_names; use_change_names=delta_PSI_names; start_head="juncpos5"; end_header="juncpos3" #for junc5 or junc3 table

##1, combine all data:
####load and combine data using file_setting_list
remove_headers=c("cds_len","transcript_len","GeneUnit_id","unit_id","coordinates","endPos_id","startPos_id","region_id")
keep_headers=c("gene_symbol","Gblock_id","contig","strand",start_head,end_header,"skip_intron", "startPosType","startPosAno","endPosType","endPosAno","reg_len","region_ano","gene_id","ss5_Ano","ss3_Ano","startSS_supp","endSS_supp", "startSS_seq","endSS_seq","exon_i_inFrame","exon_i_stopCodon",
		"refseqid","gene_desc","ups_junc_num","dns_junc_num","skip_junc_num","introni","ss5_sco","ss3_sco")
comb_d=NULL
for(i in 1:length(file_setting_list)){ 
	file_setting=file_setting_list[[i]]
	inf=file_setting$input_f
	print(paste("open file:",inf))
	d=read.table(inf, header=T, sep="\t", stringsAsFactors=F, comment.char="", quote="")
	#if("ss5" %in% names(d) & !("end_pos" %in% names(d))){ names(d)[match("ss5", names(d))]="end_pos" }
	#if("ups_ss3" %in% names(d) & !("start_pos" %in% names(d))){ names(d)[match("ups_ss3", names(d))]="start_pos" }
	#if("ss3" %in% names(d) & !("start_pos" %in% names(d))){ names(d)[match("ss3", names(d))]="start_pos" }
	d$Gblock_id=paste(d$gene_symbol, ":",d$contig, ":", d$strand, ":", d[,start_head], "-", d[,end_header], sep="")
	if("skip_intron" %in% names(d)){
		d$Gblock_id=paste(d$Gblock_id,d$skip_intron,sep=":")
	}
	control_PSI_names1=paste("PSI_",unique(file_setting$ref_samples),sep="")

	d=re_format_tb(d, delete_headers=c(remove_headers), del_pattern="^num_|Num_",  ch_header_from=file_setting$change_headerFroms, ch_header_to=file_setting$change_headerTos)
	all_data_headers=grep("deltaPSI_|Pfisher_|Pttest_|Log2R_", names(d),value=T)
	used_data_heqaders=c( outer(c("deltaPSI_","Pfisher_","Pttest_","Log2R_"), paste(file_setting$test_samples, file_setting$ref_samples,sep="_"), paste, sep="") )
	discard_data_headers=setdiff(all_data_headers, used_data_heqaders)
	if(length(discard_data_headers)>0){ d=d[,setdiff(names(d),discard_data_headers)] } #remove unnesessary data

	if(i==1){
		keep_headers2=intersect(c(keep_headers,control_PSI_names1, use_change_names, use_P_names),names(d))
		comb_d=d[keep_headers2]
	}else{
		update_headers=intersect( c(control_PSI_names1,use_P_names,use_change_names), names(d) )
		# comb_d[,update_headers]=d[match(comb_d$Gblock_id, d$Gblock_id),update_headers]
		common_headers=intersect(names(comb_d),names(d))
		comb_d=merge(comb_d, d[,c( "Gblock_id",update_headers )], by="Gblock_id", all=T)
		if_update=is.na(comb_d$gene_symbol)
		if(any(if_update)){
			keep_headers2=setdiff(intersect(keep_headers, names(d)),"Gblock_id")
			comb_d[if_update, keep_headers2]=d[match(comb_d$Gblock_id[if_update], d$Gblock_id), keep_headers2]
		}
	}
	print(nrow(comb_d))
}
rm ("d")

#comb_d[comb_d$gene_symbol=="FOXM1" & comb_d$region_ano=="exon" & comb_d$reg_len==34,]
ifsel=rowSums( !is.na(comb_d[,c(use_P_names,use_change_names)]) ) >0; table(ifsel)
comb_d=comb_d[ifsel, ]

if(all(c("contig","strand",start_head,end_header) %in% names(comb_d))){
	comb_d$Coordinates=paste(comb_d$contig,":",ifelse(comb_d$strand=="+",comb_d[,start_head], comb_d[,end_header]), "-",ifelse(comb_d$strand=="+",comb_d[,end_header], comb_d[,start_head]), sep="" )
}else{
	comb_d$Coordinates=sub("^[^:]*:","",comb_d$Gblock_id, perl=T)
	comb_d$Coordinates=sub(":\\+:|:\\-:",":",comb_d$Coordinates, perl=T)
}

#2, calculate ReguType_names
P_cut_top=0.001; change_cut=log2(2) 
P_cut_top=0.001; change_cut=40
P_cut_top=1e-7; change_cut=40
for(i in 1:ncol(sample_compare_matrix)){
	comb_d[,ReguType_names[i]]=cal_regu_type(data=comb_d, p_name=use_P_names[i], change_name=use_change_names[i], p_cut=P_cut_top, change_cut=change_cut )
}
##re-organize orders of columns
front_common_headers=c("gene_symbol", "Coordinates", setdiff(intersect(keep_headers, names(comb_d)), c("gene_symbol","gene_desc") ) )
comb_d=re_format_tb(comb_d, front_headers=c(front_common_headers,control_PSI_names, use_P_names,use_change_names,ReguType_names,"gene_desc") )

#output all data to a file (the file may be huge)
out_combine_f= paste(output_root, study_set,"/",study_set,".Gblock.Combine.",ana_meth,"All.txt" ,sep="");
mkdir_if_not_exist(out_combine_f)
write.table(comb_d, file=out_combine_f, col.names=T,row.names=F, sep="\t", quote=F, na="")


##output a file of Gblocks regulated in at least one condition:
sel_gex_ids=1:length(ReguType_names)
ifreg <- rowSums( !is.na(comb_d[ReguType_names[sel_gex_ids]]) & (comb_d[ReguType_names[sel_gex_ids]]=="UP" | comb_d[ReguType_names[sel_gex_ids]]=="DN")  )>0
sum(ifreg); 
out_format_f= paste(output_root, study_set,"/Gblock.Combine.",ana_meth,".P",P_cut_top,".Ch",round(change_cut,2),".txt" ,sep="");
mkdir_if_not_exist(out_format_f)
write.table(comb_d[ifreg,], file=out_format_f, col.names=T,row.names=F, sep="\t", quote=F)

##output a subset of exons belong to a defined gene set (eg. expression down-regulated and unannotated exons)
set_name="downRegulated_iExons"
sel_genes=unlist(strsplit("Tomm70a Mipep Memo1 Trappc3 Pcca Hsp90aa1 Rps6ka5 Ecd Adam19 Agrn Akap9 Lss Tmem57 Fkbp11 Slit3 Angel2 Srek1 Dcaf8 Ahcyl1 Slc4a8 Fbxl22 Grb14 Dynll1 Spon2 Stau2 Dok5 Lynx1 Gm6307 Adamts8 Lrtm1 Myh2 Nrep 2310001H17Rik Ptcd3 Mb Akr1e1 Ddah1 Mmaa Sirt3 Ggh Cpeb1 Akr1b10 Ramp1 Map2k6 Sypl2 Efr3b Cdc37l1 P2rx3 Zfp142 Dok4 Bhlhe41 Oxct1 Zfp385b Cab39l Ergic3 Arntl Grip2 Zfp454 Dusp10 Uqcr10 Fam195a Atp5g1 Rps6kl1 Hspa4l Dhrs1 Aplnr Nudt8 Ptrhd1 Ndufa12 Zfp750 Fam136a Tead1 Bcl6b Oxld1 Perm1 Gng5 Fam65b Snai3 Snapc3 Ndufab1 Pole3 Chchd10 9630001P10Rik 2310015D24Rik Dusp23 Plch2 Gamt Deb1 Ndufa6 Relb B230312C02Rik Cacna2d3 Aurka Phyh"," "))
if_iExon=comb_d$region_ano %in% "exon" & !(comb_d$startSS_supp %in% c("Refseq") & comb_d$endSS_supp %in% c("Refseq")); table(if_iExon)
sel_gex_ids=1:length(ReguType_names)
if_up <- rowSums( !is.na(comb_d[ReguType_names[sel_gex_ids]]) & comb_d[ReguType_names[sel_gex_ids]]=="UP"  )>0
if_output=comb_d$gene_symbol %in% sel_genes & if_iExon & if_up; table(if_output)
out_subset_f= paste(output_root, study_set,"/Gblock.Combine.",ana_meth,".P",P_cut_top,".Ch",round(change_cut,2),".",set_name, ".txt" ,sep="");
write.table(comb_d[if_output,], file=out_subset_f, col.names=T,row.names=F, sep="\t", quote=F)



#stats of regulation type
library("fmsb")
out_stats_f=paste(out_format_f,".stats.tbl",sep="")
tb=as.matrix( apply(comb_d[ReguType_names], 2, function(v){ table(factor(v,levels=ReguTypes))} ) [ReguTypes,] , nrow=length(ReguTypes) )
colnames(tb)=ReguType_names
print(tb)
write(paste("Regulation"),  file=out_stats_f)
write.table(cbind(ReguTypes,tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)

if( !("region_ano" %in% names(comb_d)) ){ comb_d$region_ano=ana_meth }
region_anos=unique(comb_d$region_ano) ; region_anos=region_anos[!is.na(region_anos)]
for(region_ano in region_anos){
	print(region_ano)
	ifsel=!is.na(comb_d$region_ano) & comb_d$region_ano==region_ano
	tb=as.matrix( apply(as.matrix(comb_d[ifsel, ReguType_names]), 2, function(v){ table(factor(v,levels=ReguTypes))} ) [ReguTypes,] , nrow=length(ReguTypes) )
	colnames(tb)=ReguType_names
	write(paste("\nRegion annotation=",region_ano),  file=out_stats_f, append=T)
	write.table(cbind(ReguTypes,tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
}
type_headers1=ReguType_names[11]  ;  type_headers2=ReguType_names[1:10]
type_headers1=ReguType_names;  type_headers2=ReguType_names
if(length(ReguType_names)>1){
for(region_ano in region_anos){
	print(region_ano)
	ifsel=!is.na(comb_d$region_ano) & comb_d$region_ano==region_ano
	compare_types_among_vars(out_stats_f, comb_d[ifsel,], type_headers1=type_headers1, Types=ReguTypes, type_headers2=type_headers2, 
		prefix=paste("\ncross table: region_ano=",region_ano) )
}}



##3, calculate overlap (Venn diagram) for activated exons (not used)
ReguType_names
P_cut_bottom=0.1; minLog2ratio_cut=log2(1.2); ajacentLog2ratio_cut=log2(2); maxLog2ratio_cut=log2(8)
P_cut_bottom=0.001; minLog2ratio_cut=log2(2); ajacentLog2ratio_cut=log2(2); maxLog2ratio_cut=log2(2)
P_cut_bottom=0.05; P_cut_top=1e-5; minLog2ratio_cut=log2(1.5); ajacentLog2ratio_cut=log2(4); maxLog2ratio_cut=log2(10)
sel_regu_type="UP"
sel_region_type=c("exon","exon_a5ss","exon_a3ss","exon_part","intron","intron_a5ss","intron_a3ss","intron_part")
sel_region_type=c("exon")

sel_ids=c(1,2); sampleSet="NVS_HTT"
sel_ids=c(1,3,6,8); sampleSet="NVS_PTC_2RO"
sel_ids=c(5:8); sampleSet="2RO_2concentration"
sel_ids=c(6,8,14); sampleSet=paste("247_067_NVS_normalCell",".P",P_cut_bottom,".Ch",round(ajacentLog2ratio_cut,2), sep="")
sel_ids=c(1,2,5,6); sampleSet=paste("247_twoHsCells_TwoConc",".P",P_cut_bottom,".Ch",round(ajacentLog2ratio_cut,2), sep="")
sel_ids=c(3,4,7,8); sampleSet=paste("067_twoHsCells_TwoConc",".P",P_cut_bottom,".Ch",round(ajacentLog2ratio_cut,2), sep="")
sel_ids=c(6,9); sampleSet=paste("247_pA_vs_rz",".P",P_cut_bottom,".Ch",round(ajacentLog2ratio_cut,2), sep="")

ifsel=comb_d$region_ano %in% sel_region_type; sum(ifsel)
comb_d1=comb_d[ifsel, ]
comb_d1=comb_d1[!duplicated(comb_d1$Gblock_id), ]
#remove NA type exons
if_notNA=rowSums(is.na(comb_d1[,use_change_names[sel_ids]]) | is.na(comb_d1[,use_P_names[sel_ids]]))==0; table(if_notNA)
comb_d1=comb_d1[if_notNA,]
##only study GA-type exon
if_GA=substring(comb_d1$endSS_seq,9,10)=="GA";
comb_d1=comb_d1[if_GA,]
if_regulatedInOne=rowSums(comb_d1[,use_change_names[sel_ids]]>maxLog2ratio_cut & comb_d1[,use_P_names[sel_ids]]<P_cut_top)>0; sum(if_regulatedInOne)
comb_d1=comb_d1[if_regulatedInOne,]

# for(i in sel_ids){ #if >minLog2ratio_cut study this exon
# 	comb_d1[,ReguType_names[i]]=cal_regu_type(data=comb_d1, p_name=use_P_names[i], change_name=use_change_names[i], p_cut=P_cut_bottom, change_cut=minLog2ratio_cut )
# }
# if_regu_1=rowSums(comb_d1[,ReguType_names[sel_ids]]==sel_regu_type)>0 ; sum(if_regu_1)
# comb_d1=comb_d1[if_regu_1,]
for(i in sel_ids){ #if >maxLog2ratio_cut set to UP 
	comb_d1[,ReguType_names[i]]=cal_regu_type(data=comb_d1, p_name=use_P_names[i], change_name=use_change_names[i], p_cut=P_cut_top, change_cut=maxLog2ratio_cut )
}

Log2Ratio_tb=comb_d1[use_change_names[sel_ids]]
Log2Ratio_tb[ is.na(comb_d1[use_P_names[sel_ids]]) | comb_d1[use_P_names[sel_ids]]>=P_cut_bottom | comb_d1[use_change_names[sel_ids]]<=minLog2ratio_cut ] =NA
#comb_d1[,ReguType_names[sel_ids]]="NC"
#every row must have at least one UP
#ifMax_FC= Log2Ratio_tb==apply(Log2Ratio_tb,1,max,na.rm=T)
#ifMax_FC[is.na(ifMax_FC)]=F
#comb_d1[,ReguType_names[sel_ids]][ifMax_FC]="UP"
for(i in 1:(length(sel_ids)-1)){
	Log2Ratio_tb2=Log2Ratio_tb
	Log2Ratio_tb2[ comb_d1[,ReguType_names[sel_ids]]!=sel_regu_type ]=NA
	newLog2ratio_cut=apply(Log2Ratio_tb2,1,min,na.rm=T)-ajacentLog2ratio_cut
	if_UP=!is.na(Log2Ratio_tb) & Log2Ratio_tb>newLog2ratio_cut
	comb_d1[,ReguType_names[sel_ids]] [if_UP] ="UP"
}
table(unlist(comb_d1[,ReguType_names[sel_ids]]))

comb_d1$combine_regu_code=apply( comb_d1[,ReguType_names[sel_ids]]==sel_regu_type, 1, paste, collapse="_" )
comb_d1$combine_regu_code=gsub("TRUE","T",comb_d1$combine_regu_code)
comb_d1$combine_regu_code=gsub("FALSE|NA","F",comb_d1$combine_regu_code)
tb=table(comb_d1$combine_regu_code); tb
table(comb_d1$combine_regu_code,comb_d1$region_ano); 
delete_headers=setdiff( grep("Padj_|P_|Log2R_|ReguType_",names(comb_d1),value=T), c(ReguType_names[sel_ids], use_P_names[sel_ids], use_change_names[sel_ids]) )
comb_d1=re_format_tb(comb_d1, delete_headers=delete_headers,back_headers=c(ReguType_names[sel_ids], use_P_names[sel_ids], use_change_names[sel_ids],"combine_regu_code","gene_desc") )


out_f1=paste(output_root,study_set,"/",sampleSet,"/Gblock.regulated.txt", sep="")
mkdir_if_not_exist(out_f1)
write.table(comb_d1, file=out_f1, col.names=T, row.names=F, sep="\t", quote=F)

venn_d=comb_d1[c("Gblock_id","gene_symbol","region_ano")]
venn_d[sample_pairs[sel_ids]]=as.numeric(comb_d1[,ReguType_names[sel_ids]]==sel_regu_type)
out_region_types=sel_region_type
for(out_region_type in out_region_types){
	venn_out_f1=paste(output_root,study_set,"/",sampleSet,"/Gblock.regulated.",out_region_type,".venn.txt", sep="")
	subvenn_d1=venn_d[venn_d$region_ano==out_region_type,]
	subvenn_gene_d=data.frame(gene_symbol=unique(subvenn_d1$gene_symbol))
	for(ReguType_name in sample_pairs[sel_ids] ){
		subvenn_gene_d[,ReguType_name]=tapply(subvenn_d1[,ReguType_name], subvenn_d1$gene_symbol, sum, na.rm=T)[subvenn_gene_d$gene_symbol]
	}
	write.table(subvenn_d1, file=venn_out_f1, col.names=T, row.names=F, sep="\t", quote=F)
	write.table(subvenn_gene_d, file=paste(venn_out_f1,".gene",sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}

##4, output sequence for sequence Logo
if_presel=substring(comb_d1$endSS_seq,9,10)=="GA"; table(if_presel)
combine_regu_codes=unique(comb_d1$combine_regu_code)
for(out_region_type in out_region_types){
	for(combine_regu_code in combine_regu_codes){
		sub_d=comb_d1[if_presel & comb_d1$region_ano==out_region_type & comb_d1$combine_regu_code==combine_regu_code,]
		if(nrow(sub_d)>0){
			for(seq_head in c("startSS_seq","endSS_seq")){
				out_fa_f=paste(output_root,study_set,"/",sampleSet,"/fasta/",out_region_type,"/",seq_head,".",combine_regu_code,".fa", sep="")
				out_img_f=paste(output_root,study_set,"/",sampleSet,"/weblogo/",out_region_type,".",seq_head ,"/",combine_regu_code,".jpeg", sep="")
				mkdir_if_not_exist(out_fa_f)
				mkdir_if_not_exist(out_img_f)
				sub_d2=sub_d[!is.na(sub_d[,seq_head]) & sub_d[,seq_head]!="",]
				if(nrow(sub_d2)>0){
					if(nchar(sub_d2[1,seq_head])>30 ){
						substring(sub_d2[,seq_head], 1, 30)=""
					}
					fa_text=paste(">",sub_d2$Gblock_id,"\n",sub_d[,seq_head],"\n", sep="")
					print(out_fa_f)
					
					write(fa_text, file=out_fa_f, ncolumns=1)
					cmd=paste("weblogo -f", out_fa_f, "-o", out_img_f, 
						"-F jpeg -A dna -s large --errorbars NO -U probability --title \"", nrow(sub_d2), out_region_type, combine_regu_code,seq_head   ,
						"\"", "--first-index",  -(nchar(sub_d2[1,seq_head])-10-1) )
					print (cmd)
					system(cmd)
				}
			}
		}

	}
}


##5, output a table for GO analysis or IPA analysis (each gene only keep one row, keep the most significant exon-type splicing regulation)
select_gblock_types=c("exon","exon_a5ss","exon_a3ss")
comb_d1=comb_d[comb_d$region_ano %in% select_gblock_types,]; nrow(comb_d1)
comb_gene_d=data.frame(gene_symbol=unique(comb_d1$gene_symbol))
for(i in 1:length(use_P_names)){
	update_names=c(use_P_names[i], use_change_names[i]); print(update_names)
	comb_d1=comb_d1[order(comb_d1[,use_P_names[i]], -abs(comb_d1[,use_change_names[i]])),]
	comb_gene_d[, update_names] = comb_d1[match(comb_gene_d$gene_symbol, comb_d1$gene_symbol) , update_names ]
}
#output

out_f1=paste(output_root,"gene.mostSig.Gblock.txt", sep="")
write.table(comb_gene_d, file=out_f1, col.names=T, row.names=F, sep="\t", quote=F)

