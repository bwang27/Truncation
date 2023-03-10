#2021/02/18 add min_read_cutoff to allow filter output clusters by minimal read counts

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../../../Pipeline/SharedCodes/Rfunc.inc.R")
define_parallel_fun(nCores=6)

study_name="2020-06-Research7_CMC5_AAVs_Seq"
samples=unlist(strsplit(as.character("SS-AAV9GTX-CMV-ReelinR3p6V5-03 SS-AAV9GTX-CMV-MECP2-03 SS-AAV9-WT-UbC-UBE3A-V5-01 SS-AAV9-WT.CMV-ReelinR3p6V5-01 SS-AAV9-WT.CMV-CDKL5-01 ssAAV9-WT.CMV-IgK-TATk-GFP-CDKL5-01 SS-AAV9-WT.CMV-ReelinR3p6V2-01 Vigene_vector_052718 Applied_Viromics_high_full_vector_073019 Paragon_demo_batch_vector Paragon_CsCl_band_two_vector Paragon_shake_flask_vigene_pRepCap")," "))
if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }
project_root=study_name
if(!is.na(args_v["project_root"])){ project_root= as.character(args_v["project_root"]) }
if(!is.na(args_v["samples"])){ samples= unlist(strsplit(as.character(args_v["samples"])," "))}

bed_in_dir=paste0(project_root, "/07_readCluster/")
out_prefix=paste0(project_root, "/08_ClusterAnnotation/hostCellDNA_anno")

if(!is.na(args_v["bed_in_dir"])){ bed_in_dir= as.character(args_v["bed_in_dir"]) }
if(!is.na(args_v["out_prefix"])){ out_prefix= as.character(args_v["out_prefix"]) }
min_read_cutoff=ifelse(is.na(args_v["min_read_cutoff"]), 1, as.numeric(args_v["min_read_cutoff"]) )

genomicAnnotation_regions = c("Promoter", "5' UTR","3' UTR", "Exon",  "Intron",  "Downstream", "Distal Intergenic")
cancer_gene_census_f="/drive2/wli/data/cancer/COSMIC/202009_download/cancer_gene_census.csv"
CGC_headers=c("Role.in.Cancer","Mutation.Types","Translocation.Partner")

mkdir_if_not_exist(out_prefix)
save.image(paste0(out_prefix, ".R.data.image"))

library("ChIPseeker")
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(org.Hs.eg.db)
library("openxlsx")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

bed_fs=paste0(bed_in_dir, samples, ".readMap.cluster.bed")
names(bed_fs)=samples
peakAnno_l=list()
for(i in 1:length(samples)){
	sample=samples[i]
	print(sample)
	peakAnno_l[[sample]]=annotatePeak(bed_fs[i], tssRegion=c(-3000, 1000), TxDb=txdb, level="gene", annoDb="org.Hs.eg.db", 
		genomicAnnotationPriority = c("Promoter", "5UTR","3UTR", "Exon",  "Intron",  "Downstream", "Intergenic"),
		overlap = "all")
}
names(peakAnno_l)

##create peakAnno_df_l (list of data.frames )
cluInfo_fs=paste0(bed_in_dir, samples, ".readMap.cluster.txt")
names(cluInfo_fs)=samples
sample_short_names=abbreviate(samples, method = "both.sides", minlength = 20 ) # left.kept  both.sides
any(duplicated(sample_short_names)) #make sure it's False

##open cancer_gene_census_f and load oncogene information
cancer_gene_census_d=read.csv(cancer_gene_census_f)

peakAnno_df_l=list()
for(i in 1:length(samples)){
	sample=samples[i]
	pval_name=paste0("Pval_",sample)
	Num_name=paste0("Num_",sample)
	clu_info_f=cluInfo_fs[sample]
	clu_info_d=read.table(clu_info_f, sep="\t", header=T, quote="", comment.char="", stringsAsFactors=F, check.names=F)
	if(sample %in% "AllSamples"){
		all_pval_names=grep("Pval_", names(clu_info_d), value=T)
		all_num_names=grep("Num_", names(clu_info_d), value=T)
		clu_info_d[,pval_name]=apply(clu_info_d[,all_pval_names, drop=F],1,function(p){ 10^mean(log10(p),na.rm=T) })
		clu_info_d[,Num_name]=apply(clu_info_d[,all_num_names, drop=F],1,sum)
	}
	peakAnno_d=as.data.frame(peakAnno_l[[sample]])
	peakAnno_d=re_format_tb(peakAnno_d, delete_headers=c("V6","V5","strand","width"), ch_header_from=c("V4","seqnames","start","end"),ch_header_to=c("name","rname","clu_from","clu_to"))
	peakAnno_d$geneStrand=ifelse(peakAnno_d$geneStrand==1, "+","-")
	peakAnno_d$annotation_region=sub(" \\(.*","",peakAnno_d$annotation)
	peakAnno_d=merge(clu_info_d, peakAnno_d, all=T)
	peakAnno_d=peakAnno_d[peakAnno_d[,Num_name]>=min_read_cutoff, ]

	peakAnno_d[,c("cancer_gene_census_tier",CGC_headers)]=cancer_gene_census_d[match(peakAnno_d$SYMBOL , cancer_gene_census_d$Gene.Symbol), c("Tier", CGC_headers)]
	peakAnno_df_l[[sample]]=peakAnno_d
	peakAnno_d=peakAnno_d[order(peakAnno_d[, pval_name]),]
	print(paste("output", sample))
	out_ano_f=paste0(out_prefix,".",sample,".xlsx")
	wb <- createWorkbook()
	addWorksheet(wb, sheetName = sample_short_names[i])
	freezePane(wb, sheet = sample_short_names[i], firstRow = TRUE, firstCol = F)
	writeDataTable(wb, sheet = sample_short_names[i], x = peakAnno_d, colNames = T, rowNames = F,  firstColumn=T)
	headerStyle <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", border = "Bottom", textRotation=90)
	setColWidths(wb, sheet = sample_short_names[i], cols =1:ncol(peakAnno_d), widths = 6)
	addStyle(wb, sample_short_names[i], style = headerStyle, rows = 1, cols=1:ncol(peakAnno_d), gridExpand = F, stack=T)
	setColWidths(wb, sheet = sample_short_names[i], cols = match(c("name","annotation","GENENAME"), names(peakAnno_d) ), widths = 20)
	conditionalFormatting(wb, sheet = sample_short_names[i], cols=match(pval_name,names(peakAnno_d)), rows=1+1:(nrow(peakAnno_d)), type = "colourScale", rule=c(0,0.1), style = c("orange","white"))
	saveWorkbook(wb, out_ano_f, overwrite = TRUE)
}


##analyze annotation region numbers
anno_region_sum_d=NULL
for(i in 1:length(samples)){
	Num_name=paste0("Num_",samples[i])
	peakAnno_d=peakAnno_df_l[[samples[i]]]
	peakAnno_d$annotation_region=factor(peakAnno_d$annotation_region, levels=genomicAnnotation_regions)
	if_oncogene=grepl("oncogene",peakAnno_d$Role.in.Cancer) ; sum(if_oncogene)
	peakAnno_d_onco=peakAnno_d[if_oncogene,]
	number_d=rbind(
		table(peakAnno_d$annotation_region)[genomicAnnotation_regions],
		table(peakAnno_d_onco$annotation_region )[genomicAnnotation_regions],
		tapply(peakAnno_d[,Num_name], peakAnno_d$annotation_region, sum)[genomicAnnotation_regions],
		tapply(peakAnno_d_onco[,Num_name], peakAnno_d_onco$annotation_region, sum)[genomicAnnotation_regions]
	)
	number_d[is.na(number_d)]=0
	sum_d=data.frame(sample=samples[i], GeneType=c("any","oncogene"), Count=rep(c("cluster#","read#"),each=2) )
	sum_d[genomicAnnotation_regions]=number_d
	if(i==1){anno_region_sum_d=sum_d}else{anno_region_sum_d=merge(anno_region_sum_d,sum_d, all=T)}
}
anno_region_sum_d=anno_region_sum_d[order(anno_region_sum_d$Count, anno_region_sum_d$GeneType),]
out_region_sum_f=paste0(out_prefix, ".summary.txt")
mkdir_if_not_exist(out_region_sum_f)
write.table(anno_region_sum_d, file=out_region_sum_f, col.names=T, row.names=F, sep="\t", quote=F)
wb <- createWorkbook()
addWorksheet(wb, sheetName = "Anno_regions.stats")
writeDataTable(wb, sheet = "Anno_regions.stats", x = anno_region_sum_d, colNames = T, rowNames = F,  firstColumn=T)
out_summary_f=paste0(out_prefix, ".summary.xlsx")
saveWorkbook(wb, out_summary_f, overwrite = TRUE)


out_img_f=paste0(out_prefix, ".barplot.pdf")
mkdir_if_not_exist(out_img_f)
pdf(out_img_f, width=7, height=7)
plotAnnoBar(peakAnno_l)
dev.off()

save.image(paste0(out_prefix, "R.data.image"))
