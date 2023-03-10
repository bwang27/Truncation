#generate read coverage plot using sashimi-plot.py
#2021/03/02 add function to plot each individual sample and plasmid
#2022/09/06 WL: changed the way to plot the genome coverage from sashimi-plot.py to a R library plotgardener to read and plot the bigwig files (plotgardener, installed on R4.2 /home/wli/packages/R-devel/bin/)
#2022/09/07 WL: need to install R library GenomicFeatures, #will use the makeTxDbFromGFF function

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
if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }

sample_info_f="/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/2020-06-Research7_CMC5_AAVs_Seq/Samples.xlsx"
if(!is.na(args_v["sample_info_f"])){ sample_info_f= as.character(args_v["sample_info_f"]) }

mapping_data_root=paste0("/home/nas02/ProcessedData/Nanopore_seq/",study_name,"/minimap2/" )
if(!is.na(args_v["mapping_data_root"])){ mapping_data_root= as.character(args_v["mapping_data_root"]) }

out_folder=ifelse(is.na(args_v["out_folder"]), paste0(study_name, "/04_readCov_plot/"), args_v["out_folder"])
GOI_fa_dir=ifelse(is.na(args_v["GOI_fa_dir"]), "../../../summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences", args_v["GOI_fa_dir"])
index_folder=ifelse(is.na(args_v["index_folder"]), "target_index", args_v["index_folder"])
#bedpe_folder="/home/Research_Archive/ProcessedData/Nanopore_seq/2022-04-19_AAV_SEQ/02b_ChimericRead/visu/"
bedpe_folder=ifelse(is.na(args_v["bedpe_folder"]), paste0(study_name, "/02b_ChimericRead/visu/"), args_v["bedpe_folder"])

##default plot parameters (in inches)
plot_width=8
archPlot_width=10
signalPlot_height=1.5
archPlot_height=1
textPlot_height=0.2
left_margin=0.5
page_top_bottom_margin=0.1
xaxis_height=0.4
feature_box_height=0.1

library("openxlsx")
#library("GenomicFeatures")
library("plotgardener")
#library("OrganismDbi")
#library("GO.db")
#library("IRanges")
#library("org.Hs.eg.db")
#library("TxDb.Hsapiens.UCSC.hg19.knownGene")
chimeric_types=list(
	Inversion=c("Inversion+","Inversion-"), 
	Deletion=c("Deletion+","Deletion-"),
	Insertion=c("Insertion+","Insertion-"),
	Reversion=c("Reversion+","Reversion-"),
	fusion="fusion")
archLinecolors=c("black","red","blue"); names(archLinecolors)=c("*","+","-")
axis_range2bestLabel<-function(maxVal,tickNum=3){ #this function is mainly for counts values (such as read counts, the lower range is always 0)
	if(maxVal<10){
		label_axis_at=c(0,maxVal)
	}else{
		binsize=round(maxVal/tickNum)
		binsize_numDigs=nchar(binsize)
		binsize2=round(binsize/10^(binsize_numDigs-1))*10^(binsize_numDigs-1)
		label_axis_at=seq(0,maxVal, by=binsize2)
	}
	label_axis_at
}

modify_fusion_bedpe<-function(bedpe_d, keepChro, ChroSize){ #modify fusion bedpe for visualization only
	bedpe_d=bedpe_d[bedpe_d[,1] %in% keepChro | bedpe_d[,4] %in% keepChro, ]
	if_to_GOI=bedpe_d[,4] %in% keepChro #from others to GOI
	bedpe_d2=bedpe_d
	bedpe_d2[if_to_GOI,c(1,2,3,9)]=bedpe_d[if_to_GOI,c(4,5,6,10)]
	bedpe_d2[if_to_GOI,c(4,5,6,10)]=bedpe_d[if_to_GOI,c(1,2,3,9)]
	bedpe_d2[,4]=bedpe_d2[,1]
	bedpe_d2[,10]=bedpe_d2[,9]
	if_plus=bedpe_d2[,9]=="+"
	bedpe_d2[if_plus,6]=ChroSize+(ChroSize-bedpe_d2[if_plus,3])
	bedpe_d2[if_plus,5]=bedpe_d2[if_plus,6]-1
	if_minus=bedpe_d2[,9]=="-"
	bedpe_d2[if_minus,5:6]= bedpe_d2[if_minus,2:3]
	bedpe_d2[if_minus,2:3]= -bedpe_d2[if_minus,c(3,2)]
	bedpe_d2
}

sample_info_d=read.xlsx(sample_info_f)
sample_info_d$index_folder=apply(sample_info_d[,c("pHelper","pRepCap","GOI")],1,paste,collapse=".")
sample_info_d$fa_sizes_f=paste0(index_folder,"/",sample_info_d$index_folder,"/", sample_info_d$index_folder,".fa.sizes")
sample_info_d$Ad5="Human_Ad5_first5k"
sample_info_d$Ad5[is.na(sample_info_d$pHelper)]=NA
sample_info_d$bigwig_f=paste0(mapping_data_root, "/",sample_info_d[,1],".bw")
mkdir_if_not_exist(out_folder)


##1, read reference chromosome/plasmid sizes
chr_sizes_d=NULL
for(i in 1:nrow(sample_info_d)){
	print(sample_info_d$fa_sizes_f[i])
	d1=read.table(sample_info_d$fa_sizes_f[i], header=F, sep="\t", quote="",stringsAsFactors=F)
	names(d1)=c("chr","size")
	if(i==1){
		chr_sizes_d=d1
	}else{
		chr_sizes_d=merge(chr_sizes_d,d1,all=T)
	}
}
save.image(paste0(out_folder,"R.data.image"))


#for each plasmid, find all samples used that vector and plot genome coverage using sashimi-plot.py (old way)
# for (plasmid_type in c("pHelper","pRepCap","GOI","Ad5")){
# 
# 	plasmids=setdiff( unique(sample_info_d[,plasmid_type]),c("NA","na",NA) )
# 	for(plasmid1 in plasmids){
# 		sample_info_d_pla1=sample_info_d[sample_info_d[,plasmid_type] %in% plasmid1, ]
# 		
# 		plot_sample_l=c(list(1:nrow(sample_info_d_pla1)), 1:nrow(sample_info_d_pla1) )
# 		names(plot_sample_l)=c("All", unlist(sample_info_d_pla1[,1]) )
# 		for(s in names(plot_sample_l)){
# 			sample_info_d1=sample_info_d_pla1[plot_sample_l[[s]], ]
# 			gtf_subfolder=ifelse(plasmid_type %in% "GOI",sample_info_d1[1,"GOI_folder"]  ,sub("^p","",plasmid_type))
# 			gtf_f=paste0(GOI_fa_dir,"/", gtf_subfolder,"/", sample_info_d1[1,plasmid_type], ".gtf")
# 			gtf_d=read.table(gtf_f, header=F, sep="\t", quote="")
# 			feature_num=nrow(gtf_d)
# 			bam_list_d=data.frame(sample=sample_info_d1[,1], bam_f=paste0(mapping_data_root,"/minimap2/",sample_info_d1[,1] ,".sorted.bam"), id=1:nrow(sample_info_d1))
# 			palette_d=data.frame(color=sample_info_d1$Color)
# 			bam_list_f=paste0(out_folder,"/",s,"/", plasmid1,".bam_list.tsv")
# 			palette_f=paste0(out_folder,"/",s,"/", plasmid1,".palette.tsv")
# 			mkdir_if_not_exist(bam_list_f)
# 			write.table(bam_list_d, file=bam_list_f, col.names=F, row.names=F, sep="\t", quote=F)
# 			write.table(palette_d, file=palette_f, col.names=F, row.names=F, sep="\t", quote=F)
# 			plasmid_size=chr_sizes_d[match(plasmid1, chr_sizes_d[,1]),2]
# 			plot_region=paste0(plasmid1,":1-",plasmid_size)
# 			out_prefix=paste0(out_folder,"/",s,"/",plot_region,".readCoverage")
# 			sashimi_cmd=paste("sashimi-plot.py -b ",bam_list_f," --gtf ",gtf_f," -c ",plot_region," -o ",out_prefix,
# 				" --height 1.4 --ann-height ", feature_num/15*2 ," --palette ",palette_f," --width 12 --base-size 12 --color-factor 3 --out-format pdf")
# 			system(sashimi_cmd)
# 			convert_cmd=paste0("convert  -density 150 ",out_prefix,".pdf ",out_prefix,".png")
# 			system(convert_cmd)
# 
# 		}
# 	}
# }

#for each plasmid, find all samples used that vector and plot genome coverage using plotgardener (WL, 2022/09)
##2, read all bigwig files:
readCov_l=list()
for(i in 1:nrow(sample_info_d)){
	print(sample_info_d$bigwig_f[i])
	bigwig_d=readBigwig(file = sample_info_d$bigwig_f[i]) #all converted to 1-based coordinates
	readCov_l[[sample_info_d[i,1]]]=bigwig_d
}

##3, make read coverage plots
for (plasmid_type in c("pHelper","pRepCap","GOI","Ad5")){
	plasmids=setdiff( unique(sample_info_d[,plasmid_type]),c("NA","na",NA) )
	for(plasmid1 in plasmids){
		plasmid_size=chr_sizes_d[match(plasmid1, chr_sizes_d[,1]),2]
		sample_info_d_pla1=sample_info_d[sample_info_d[,plasmid_type] %in% plasmid1, ]

		# gtf_subfolder=ifelse(plasmid_type %in% "GOI",sample_info_d_pla1[1,"GOI_folder"]  ,sub("^p","",plasmid_type))
		# gtf_f=paste0(GOI_fa_dir,"/", gtf_subfolder,"/", sample_info_d_pla1[1,plasmid_type], ".gtf")
		# gtf_d=read.table(gtf_f, header=F, sep="\t", quote="")
		# gtf_d[,9]=gsub("\\'","p",gtf_d[,9])
		bed_subfolder=ifelse(plasmid_type %in% "GOI",sample_info_d_pla1[1,"GOI_folder"]  ,sub("^p","",plasmid_type))
		bed_f=paste0(GOI_fa_dir,"/", bed_subfolder,"/", sample_info_d_pla1[1,plasmid_type], ".bed")
		bed_d=read.table(bed_f, header=F, sep="\t", quote="", stringsAsFactors=F)
		names(bed_d)=c("chromosome","start","end","name","score","strand")
		if_GOI=grepl("GOI", bed_d[,4])
		if(any(if_GOI)){
			bed_d[if_GOI,4]=paste0(bed_d[if_GOI,4],"(",bed_d[if_GOI,3]-bed_d[if_GOI,2],")")
		}
		
		feature_num=nrow(bed_d)
		feature_plot_total_height=feature_num*feature_box_height
		#chrominfo=data.frame(chrom = plasmid1, length=plasmid_size, is_circular=TRUE)
		#FeatureTxDB=makeTxDbFromGFF(gtf_f, chrominfo=chrominfo, organism=NA, taxonomyId=NA) #, "Homo sapiens" format="gtf"
		#odb <- makeOrganismDbFromTxDb(txdb=FeatureTxDB)
		#asGTF(FeatureTxDB)
		#asGFF(TxDb.Hsapiens.UCSC.hg19.knownGene)
		#plotAssembly=assembly("AAVgenome", TxDb=FeatureTxDB, OrgDb = NA, gene.id.column="ENTREZID", display.column="ENTREZID") #, OrgDb=odb, "org.Hs.eg.db" , , 
		plot_sample_l=c(list(1:nrow(sample_info_d_pla1)), 1:nrow(sample_info_d_pla1) )
		names(plot_sample_l)=c("All", unlist(sample_info_d_pla1[,1]) )
		
		for(s in names(plot_sample_l)){
			sample_info_d1=sample_info_d_pla1[plot_sample_l[[s]], ]
 			plot_region=paste0(plasmid1,":1-",plasmid_size)
 			out_prefix=paste0(out_folder,"/",s,"/",plot_region,".readCoverage")
			#out_prefix=paste0(out_folder,"/",s,"/",plasmid1,".readCoverage")
			mkdir_if_not_exist(out_prefix)
			if(s !="All" & length(plot_sample_l[["All"]])==1){
				if(plot_sample_l[["All"]] == plot_sample_l[[s]]){ #if the current sample is the same as All samples, copy the image files from All folder
					for(fileext in c("pdf","png")){
						cp_cmd=paste0("cp -f \'", out_folder,"/All/",plot_region,".readCoverage.",fileext,"\' \'",out_prefix,".",fileext,"\'")
						system(cp_cmd)
					}
				}
			}else{
				plot_total_height=page_top_bottom_margin*2+(signalPlot_height+textPlot_height*2)*nrow(sample_info_d1)+xaxis_height+feature_plot_total_height
				
				pdf(paste0(out_prefix,".pdf"), width=plot_width, height=plot_total_height)
				pageCreate(width = plot_width, height = plot_total_height, default.units = "inches", showGuides=F)
				for(rowi in 1:nrow(sample_info_d1)){
					readCov_d1=readCov_l[[sample_info_d1[rowi,1]]]
					plotText(label=sample_info_d1[rowi,1], fontcolor=sample_info_d1$Color[rowi], 
						x=(plot_width-left_margin)/2, y=page_top_bottom_margin+(signalPlot_height+textPlot_height*2)*(rowi-1),
						just=c("center","top") )
					ymax=ifelse(any(readCov_d1[,1]==plasmid1), max(readCov_d1[readCov_d1[,1]==plasmid1, 6]),  1)
					label_y_at=axis_range2bestLabel(ymax,tickNum=3)
					label_y_at_label=label_y_at
					if(max(label_y_at)>=1000){label_y_at_label=format(label_y_at, big.mark=",")}
					signal1plot <- plotSignal( data = readCov_d1, chrom=plasmid1, chromstart=1, chromend=plasmid_size, linecolor=sample_info_d1$Color[rowi],
						scale=T, cex = 0.8, x = left_margin, 
						y = page_top_bottom_margin+(signalPlot_height+textPlot_height*2)*(rowi-1) + textPlot_height, 
						width = plot_width-left_margin, height = signalPlot_height, 
						just = c("left", "top"), lwd=2, default.units = "inches")
					annoYaxis(plot = signal1plot, at=label_y_at, label=label_y_at_label, fontsize = 7, axisLine=T)
				}
				#plot chromosome coordinates:
				annoGenomeLabel( plot = signal1plot, scale = "bp", x = left_margin, y = page_top_bottom_margin+(signalPlot_height+textPlot_height*2)*nrow(sample_info_d1), 
					sequence=F, at=c(1,(1:as.integer(plasmid_size/1000))*1000), just = c("left", "top"), default.units = "inches")

				#plot gene features
				for(i in 1:nrow(bed_d)){
					boxLeft=left_margin+(bed_d[i,2]+1)/plasmid_size*(plot_width-left_margin)
					boxWidth=(bed_d[i,3]-bed_d[i,2])/plasmid_size*(plot_width-left_margin)
					y=page_top_bottom_margin+(signalPlot_height+textPlot_height*2)*nrow(sample_info_d1)+xaxis_height+(i-1)*feature_box_height
					plotRect(x=boxLeft, width=boxWidth,
						y=y, height=feature_box_height-0.01,
						just = c("left", "top"), lwd=NA, fill=ifelse(bed_d[i,6]=="-","blue","red"), alpha=ifelse(grepl("GOI", bed_d[i,4]),1,0.5)
					)
					labelOnLeft=boxLeft>plot_width-boxLeft-boxWidth #label on left or right based on the size of open space
					plotText(bed_d[i,4], fontsize=7, 
						x=ifelse(labelOnLeft, boxLeft, boxLeft+boxWidth),
						y=y, just=c(ifelse(labelOnLeft, "right","left"),"top")
					)
				}

				dev.off()
				convert_cmd=paste0("convert  -density 150 ",out_prefix,".pdf ",out_prefix,".png")
				system(convert_cmd)
			}
			
		}
		
	}
}

##4, plot the bedpe files for chimeric read analysis:

for (plasmid_type in c("pHelper","pRepCap","GOI","Ad5") ){
	plasmids=setdiff( unique(sample_info_d[,plasmid_type]),c("NA","na",NA) )
	for(plasmid1 in plasmids){
		plasmid_size=chr_sizes_d[match(plasmid1, chr_sizes_d[,1]),2]
		sample_info_d_pla1=sample_info_d[sample_info_d[,plasmid_type] %in% plasmid1, ]

		bed_subfolder=ifelse(plasmid_type %in% "GOI",sample_info_d_pla1[1,"GOI_folder"]  ,sub("^p","",plasmid_type))
		bed_f=paste0(GOI_fa_dir,"/", bed_subfolder,"/", sample_info_d_pla1[1,plasmid_type], ".bed")
		bed_d=read.table(bed_f, header=F, sep="\t", quote="", stringsAsFactors=F)
		names(bed_d)=c("chromosome","start","end","name","score","strand")
		if_GOI=grepl("GOI", bed_d[,4])
		if(any(if_GOI)){
			bed_d[if_GOI,4]=paste0(bed_d[if_GOI,4],"(",bed_d[if_GOI,3]-bed_d[if_GOI,2],")")
		}
		feature_num=nrow(bed_d)
		feature_plot_total_height=feature_num*feature_box_height
		arch_plot_total_height=(archPlot_height+textPlot_height*2)*length(chimeric_types)
		
		for(rowi in 1:nrow(sample_info_d_pla1)){
			s=sample_info_d_pla1[rowi,1] #sample name
 			plot_region=paste0(plasmid1,":1-",plasmid_size)
 			out_prefix=paste0(out_folder,"/",s,"/",plot_region,".chimeric"); print(out_prefix)
			mkdir_if_not_exist(out_prefix)

			plot_total_height=page_top_bottom_margin*2+(signalPlot_height+textPlot_height*2)+arch_plot_total_height+xaxis_height+feature_plot_total_height
			
			pdf(paste0(out_prefix,".pdf"), width=archPlot_width, height=plot_total_height)
			pageCreate(width = archPlot_width, height = plot_total_height, default.units = "inches", showGuides=F)
				#read coverage plot 
				readCov_d1=readCov_l[[s]]
				plotText(label=s, fontcolor=sample_info_d_pla1$Color[rowi], 
					x=(archPlot_width-left_margin)/2, y=page_top_bottom_margin,
					just=c("center","top") )
				ymax=ifelse(any(readCov_d1[,1]==plasmid1), max(readCov_d1[readCov_d1[,1]==plasmid1, 6]),  1)
				label_y_at=axis_range2bestLabel(ymax,tickNum=3)
				label_y_at_label=label_y_at
				if(max(label_y_at)>=1000){label_y_at_label=format(label_y_at, big.mark=",")}
				signal1plot <- plotSignal( data = readCov_d1, chrom=plasmid1, chromstart=1, chromend=plasmid_size, linecolor=sample_info_d_pla1$Color[rowi],
					scale=T, cex = 0.8, x = left_margin, 
					y = page_top_bottom_margin+ textPlot_height, 
					width = archPlot_width-left_margin, height = signalPlot_height, 
					just = c("left", "top"), lwd=2, default.units = "inches")
				annoYaxis(plot = signal1plot, at=label_y_at, label=label_y_at_label, fontsize = 7, axisLine=T)

				#plotbedpe files for chimeric reads (Arch plot)
				for(chimeric_typei in 1:length(chimeric_types)){
					chimeric_type=names(chimeric_types)[chimeric_typei]
					bedpe_fs=paste0(bedpe_folder, "/", s, ".", chimeric_types[[chimeric_type]], ".bedpe") #bedpe files for the sample
					archSubStrandPlot_height=archPlot_height/2
					y_adjust=0
					if(chimeric_type=="fusion"){
						if(file.exists(bedpe_fs)){
							bedpe_d=read.table(bedpe_fs, header=F, sep="\t", quote="", stringsAsFactors=F)
							bedpe_d2=modify_fusion_bedpe(bedpe_d, plasmid1, plasmid_size)
						}else{
							print(paste0("Error/Warning: bedpe_f ",bedpe_fs," not found!"))
							bedpe_d2=NULL
						}
						bedpe_dl=list(bedpe_d2[bedpe_d2[,9]=="+",], bedpe_d2[bedpe_d2[,9]=="-",])
					}else{
						bedpe_dl=list()
						for(filei in 1:length(bedpe_fs)){
							bedpe_f=bedpe_fs[filei]
							if(file.exists(bedpe_f)){
								bedpe_d=read.table(bedpe_f, header=F, sep="\t", quote="", stringsAsFactors=F)
								bedpe_d=bedpe_d[bedpe_d[,1] %in% plasmid1 | bedpe_d[,4] %in% plasmid1,]
							}else{
								bedpe_d=NULL
								print(paste0("Error/Warning: bedpe_f ",bedpe_f," not found!"))
							}
							bedpe_dl[[filei]]=bedpe_d
						}
					}
					names(bedpe_dl)=c("+","-")
					chimeric_read_counts=sapply(bedpe_dl, function(bedpe_d){sum(bedpe_d[,8],na.rm=T)})

					plotText(label=paste0(chimeric_type," (read counts=",sum(chimeric_read_counts), ": " ,
							paste( paste(names(chimeric_read_counts),chimeric_read_counts,sep=""), collapse=" "), ")" ), 
						fontcolor="black", 
						x=(archPlot_width-left_margin)/2, 
						y=page_top_bottom_margin+(signalPlot_height+textPlot_height*2)+(archPlot_height+textPlot_height*2)*(chimeric_typei-1),
						just=c("center","top") )

					for(filei in 1:2){
						y_adjust=ifelse(filei==2, archSubStrandPlot_height, 0)
						bedpe_d=bedpe_dl[[filei]]
						ymax=max(c(bedpe_d[,8],1)) #max read counts
						strand=names(bedpe_dl)[filei]
						ArchPlot=plotPairsArches(bedpe_d, chrom=plasmid1, chromstart=1, chromend=plasmid_size, 
							flip=strand=="-", archHeight="V8", fill=archLinecolors[strand], range=c(0, ymax),
							linecolor="fill", alpha=1, x=left_margin, 
							y=page_top_bottom_margin+(signalPlot_height+textPlot_height*2)+(archPlot_height+textPlot_height*2)*(chimeric_typei-1) + textPlot_height+y_adjust,
							width = archPlot_width-left_margin, height = archSubStrandPlot_height, just=c("left", "top"), default.units="inches"
						)
						annoYaxis(plot = ArchPlot, at=c(0, ymax),  fontsize = 7, axisLine=T)
					}
					
				}

			#plot chromosome coordinates:
			annoGenomeLabel( plot = signal1plot, scale = "bp", x = left_margin, y = page_top_bottom_margin+(signalPlot_height+textPlot_height*2)+arch_plot_total_height, 
				sequence=F, at=c(1,(1:as.integer(plasmid_size/1000))*1000), just = c("left", "top"), default.units = "inches")
			#plot gene features
			for(i in 1:nrow(bed_d)){
				boxLeft=left_margin+(bed_d[i,2]+1)/plasmid_size*(archPlot_width-left_margin)
				boxWidth=(bed_d[i,3]-bed_d[i,2])/plasmid_size*(archPlot_width-left_margin)
				y=page_top_bottom_margin+(signalPlot_height+textPlot_height*2)+arch_plot_total_height+xaxis_height+(i-1)*feature_box_height
				plotRect(x=boxLeft, width=boxWidth,
					y=y, height=feature_box_height-0.01,
					just = c("left", "top"), lwd=NA, fill=ifelse(bed_d[i,6]=="-","blue","red"), alpha=ifelse(grepl("GOI", bed_d[i,4]),1,0.5)
				)
				labelOnLeft=boxLeft>archPlot_width-boxLeft-boxWidth #label on left or right based on the size of open space
				plotText(bed_d[i,4], fontsize=7, 
					x=ifelse(labelOnLeft, boxLeft, boxLeft+boxWidth),
					y=y, just=c(ifelse(labelOnLeft, "right","left"),"top")
				)
			}
			dev.off()
			convert_cmd=paste0("convert  -density 150 ",out_prefix,".pdf ",out_prefix,".png")
			system(convert_cmd)
			
		}
		
	}
}
