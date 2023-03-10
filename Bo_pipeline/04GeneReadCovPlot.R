#generate read coverage plot using sashimi-plot.py
#2021/03/02 add function to plot each individual sample and plasmid

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("/drive2/wli/analyze/Pipeline/SharedCodes/Rfunc.inc.R")

study_name="2020-06-Research7_CMC5_AAVs_Seq"
if(!is.na(args_v["study_name"])){ study_name= as.character(args_v["study_name"]) }

sample_info_f="/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/2020-06-Research7_CMC5_AAVs_Seq/Samples.xlsx"
if(!is.na(args_v["sample_info_f"])){ sample_info_f= as.character(args_v["sample_info_f"]) }

mapping_data_root=paste0("/home/nas02/ProcessedData/Nanopore_seq/",study_name )
if(!is.na(args_v["mapping_data_root"])){ mapping_data_root= as.character(args_v["mapping_data_root"]) }

out_folder=ifelse(is.na(args_v["out_folder"]), paste0(study_name, "/04_readCov_plot/"), args_v["out_folder"])
GOI_fa_dir=ifelse(is.na(args_v["GOI_fa_dir"]), "/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences", args_v["GOI_fa_dir"])
index_folder=ifelse(is.na(args_v["index_folder"]), "target_index", args_v["index_folder"])


library("openxlsx")
sample_info_d=read.xlsx(sample_info_f)
sample_info_d$index_folder=apply(sample_info_d[,c("pHelper","pRepCap","GOI")],1,paste,collapse=".")
sample_info_d$fa_sizes_f=paste0(index_folder,"/",sample_info_d$index_folder,"/", sample_info_d$index_folder,".fa.sizes")
sample_info_d$Ad5="Human_Ad5_first5k"
sample_info_d$Ad5[is.na(sample_info_d$pHelper)]=NA
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

#for each plasmid, find all samples used that vector and plot genome coverage
for (plasmid_type in c("pHelper","pRepCap","GOI","Ad5")){

	plasmids=setdiff( unique(sample_info_d[,plasmid_type]),c("NA","na",NA) )
	for(plasmid1 in plasmids){
		sample_info_d_pla1=sample_info_d[sample_info_d[,plasmid_type] %in% plasmid1, ]
		
		plot_sample_l=c(list(1:nrow(sample_info_d_pla1)), 1:nrow(sample_info_d_pla1) )
		names(plot_sample_l)=c("All", unlist(sample_info_d_pla1[,1]) )
		for(s in names(plot_sample_l)){
			sample_info_d1=sample_info_d_pla1[plot_sample_l[[s]], ]
			gtf_subfolder=ifelse(plasmid_type %in% "GOI",sample_info_d1[1,"GOI_folder"]  ,sub("^p","",plasmid_type))
			gtf_f=paste0(GOI_fa_dir,"/", gtf_subfolder,"/", sample_info_d1[1,plasmid_type], ".gtf")
			gtf_d=read.table(gtf_f, header=F, sep="\t", quote="")
			feature_num=nrow(gtf_d)
			bam_list_d=data.frame(sample=sample_info_d1[,1], bam_f=paste0(mapping_data_root,"/minimap2/",sample_info_d1[,1] ,".sorted.bam"), id=1:nrow(sample_info_d1))
			palette_d=data.frame(color=sample_info_d1$Color)
			bam_list_f=paste0(out_folder,"/",s,"/", plasmid1,".bam_list.tsv")
			palette_f=paste0(out_folder,"/",s,"/", plasmid1,".palette.tsv")
			mkdir_if_not_exist(bam_list_f)
			write.table(bam_list_d, file=bam_list_f, col.names=F, row.names=F, sep="\t", quote=F)
			write.table(palette_d, file=palette_f, col.names=F, row.names=F, sep="\t", quote=F)
			plasmid_size=chr_sizes_d[match(plasmid1, chr_sizes_d[,1]),2]
			plot_region=paste0(plasmid1,":1-",plasmid_size)
			out_prefix=paste0(out_folder,"/",s,"/",plot_region,".readCoverage")
			sashimi_cmd=paste("sashimi-plot.py -b ",bam_list_f," --gtf ",gtf_f," -c ",plot_region," -o ",out_prefix,
				" --height 1.4 --ann-height ", feature_num/15*2 ," --palette ",palette_f," --width 12 --base-size 12 --color-factor 3 --out-format pdf")
			system(sashimi_cmd)
			convert_cmd=paste0("convert  -density 150 ",out_prefix,".pdf ",out_prefix,".png")
			system(convert_cmd)

		}
	}
}
