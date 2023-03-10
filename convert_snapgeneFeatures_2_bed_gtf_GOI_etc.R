##2020/05/28: using the feature table exported from Snapgene, reformat to bed and gtf format, also generate a file to define gene-of-interest (GOI) region.
##2020/06/24: added function to output GOI sequences
##2021/03/01: added function to report additional features (eg. KanR, AmpR etc) to GOI_info_d
##2021/05/18: added function to output a refflat format file (for RNAseq analysis)

setwd("/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences")
library("openxlsx")
library("Biostrings") #has readDNAStringSet function

sel_plasmid_names=NULL
sel_plasmid_names=c("pscAAV.MeP426-hMECP2_e1.RDH1pA-shift1190","pAAV-UBE3A-KanR-3-shift1400","pscAAV.CMV-hMECP2_e1-shift1190","pscAAV.MeP426-hMECP2_e1-shift1190","pAAV-FXN026-AmpR-shift1383")
sel_plasmid_names=c("pAAV-FXN026-AmpR-shift1383")
sel_plasmid_names=c("pAAV.hSyn-hCDKL5_1_optimized4human-WPRE_shift","pAAV-GTX-KanR","pALD-X80")
sel_plasmid_names=c("pAAV-PTC001-eGFP_v081920_shift")
sel_plasmid_names=c("pAAV-PTC002-eGFP_shift","pAAV-PTC004-eGFP_shift","pAAV-PTC005-eGFP_shift")
sel_plasmid_names=c("pDG-KanR","pAAV-hDDC-KanR-shift")
sel_plasmid_names=c("pAAV-GFP-shift","pAAV-PTC001-temp-eGFP-v081920-shift","pAAV-PTC003-eGFP-shift","PTC-020-eGFP-filler-v102320-shift")
sel_plasmid_names=c("lambda_J02459")
sel_plasmid_names=c("FA-pAAV-FXN026-KanR-Shift1400")
sel_plasmid_names=NULL
sel_plasmid_names="Human_Ad5_first5k" #use this to define specific plasmids to analyze, set NULL to analyze all available ones
sel_plasmid_names=c("ID23_scAAV-SCA3-4E10x2")
sel_plasmid_names=c("NEUROD1_cds","pAAV9-WT","pHelper","Human_Ad5_first5k")
sel_plasmid_names=c("pDG-KanR")
sel_plasmid_names=c("pAAV-hDDC-KanR-shift")
sel_plasmid_names=c("pcDNA3.1-Hygro-E2A-rE4-VA-KanR","pcDNA3.1-Hygro-RC9","pcDNA3.1-Hygro-RC9-SIE4_WL")
sel_plasmid_names=unlist(strsplit("pAAV-GFP pAAV-GFP-ON3T11pA7-UU-CSTF3 pAAV-GFP-ON3T11pA7-hGHpA pAAV-GFP-ON3T11pA7"," "))
sel_plasmid_names=c("ID23_scAAV-SCA3-4E10x2")
sel_plasmid_names=c("scAAV-SCA3-E10x2_CAGdelIR3")
sel_plasmid_names=c("scAAV-SCA3-E10x2_hPGKp")
sel_plasmid_names=c("pHelper-KanR")
sel_plasmid_names=c("scAAV-SCA3-SR3-E10E8_hPGKp_optCDS")
sel_plasmid_names=c("pAAV-CMV-R3p6-v2-notag-KanRVERIFIED_shift")
sel_plasmid_names=c("pscAAV-hARC-V5-BMI1_TSS1_shift")
sel_plasmid_names=c("pscAAV-RK-hybrid-intron-hRP2-BGHpA-KanR_WT_shift")
sel_plasmid_names=c("Stuffer-U7hUsh2aEX13-46_shift")
sel_plasmid_names=c("Stuffer-U7hUsh2aEX13-48_shift")
sel_plasmid_names=c("scAAV-SCA3-E10x2_CAGdel7")
sel_plasmid_names=c("pAAV-RC2")
sel_plasmid_names=c("pAdDeltaF6")
sel_plasmid_names=c("pAAV.hSyn-mCDKL5_1-optM-V5-WPRE_shift")
sel_plasmid_names=c("pSTUB-1_shift")
sel_plasmid_names=c("pAAV-PTC001-eGFP-temp_shift")
sel_plasmid_names=c("pscAAV.CMV-Intron-GFP_shift")
sel_plasmid_names=c("pscAAV.CMV-GFP-WPRE_shift")
sel_plasmid_names=c("pAAV.CMV-GFP-WPRE_shift")
sel_plasmid_names=c("pAAV.CMV-GFP-T2A-Nluc-WPRE_shift")
sel_plasmid_names=c("pAAV-sc-MeP426-miniMECP2-V5-reg2-RDH1PA_shift")
sel_plasmid_names=c("pscAAV-RK-hybrid-intron-hRP2-BGHpA-KanR_WT_shift")
sel_plasmid_names=c("pAAV.hSyn-hCDKL5_1_optimized_for_human-V5-WPRE_shifted")
sel_plasmid_names=c("pscAAV.MeP426-MECP2_e1-7-2_shift")
sel_plasmid_names=c("pscAAV.MeP426-MECP2_e1-7-3_shift")
sel_plasmid_names=c("pAAV-sc-MeP426-miniMECP2-V5-reg2-RDH1PA_shift")
sel_plasmid_names=c("pAAV-Rep2rh10-KanR")
sel_plasmid_names=c("Stuffer-U7mUsh2aEX12-48_shift")
sel_plasmid_names=c("pscAAV-hARC-V5_shift")
sel_plasmid_names=c("pAAV-FXN026-KanR_CMC_shift")
#sel_plasmid_names=c("Cap9")
#sel_plasmid_names=c("pAAV-Reelin-R3p6-v2-v5-FINAL_shift")
#sel_plasmid_names=c("pHelper-KanR_CB")
#sel_plasmid_names=c("pAAV9-WT-KanR")

Out_addi_feature_patt="NeoR|KanR|AmpR|E1[AB]" #text pattern of additional features (eg. KanR, AmpR etc) to be output to GOI_info_d
master_tb_f="plasmid_files_metaData.xlsx"
master_tb_d=read.xlsx(master_tb_f, sheet=1)

processed_files=NULL
for(rowi in 1:nrow(master_tb_d)){
	feature_f=paste0(master_tb_d[rowi,"Folder"],"/Features from ",master_tb_d[rowi,"plasmid.name"],".txt")
	if(file.exists(feature_f) & !is.na(master_tb_d[rowi,"plasmid.name"]) & ! ( feature_f %in% processed_files) ){
		if(!is.null(sel_plasmid_names) & !(master_tb_d[rowi,"plasmid.name"] %in% sel_plasmid_names) ){
			print (paste("Plasmid not defined in sel_plasmid_names, not processed:", feature_f))
		}else{
			processed_files=c(processed_files, feature_f)
			print (paste0("Process file ",length(processed_files),": ", feature_f))
			feature_d=read.table(feature_f, header=F, sep="\t", quote="", comment.char="", stringsAsFactors=F, strip.white=T)
			names(feature_d)=c("name","position","length","strand","type")
			feature_d$name=gsub(" ",".",feature_d$name)
			feature_d$name=make.unique(feature_d$name)
			feature_d$position=gsub(",","",feature_d$position)
			feature_d$strand[feature_d$strand %in% c("=>","==")]="+"
			feature_d$strand[feature_d$strand %in% "<="]="-"
			feature_d$from=as.numeric(sub("\\.\\..*$","",feature_d$position ))
			feature_d$to=as.numeric(sub("^.*\\.\\.","",feature_d$position ))
			feature_d=feature_d[,setdiff(names(feature_d),c("position","type") )]
			feature_d$chromosome=master_tb_d[rowi,"plasmid.name"]
			
			out_GOI_info_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".GOI_info.txt")
			GOI_info_d=setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("chromosome","from","to", "name"))
			##add "GOI" feature for GOI vector
			if(master_tb_d[rowi,"plasmid.category"] == "GOI"){
				if_ITR_row=grepl("ITR",feature_d$name)
				sel_ITR_2rows=range((1:nrow(feature_d))[if_ITR_row])
				if(sum(if_ITR_row)>=2 & length(sel_ITR_2rows)==2){
					GOI_from_to=c( feature_d[sel_ITR_2rows[1],"from"], feature_d[sel_ITR_2rows[2],"to"] )
					feature_d=rbind(feature_d[sel_ITR_2rows[2],], feature_d)
					feature_d[1,"name"]="GOI_full"
					feature_d[1,c("from","to","length")]=c(GOI_from_to, GOI_from_to[2]-GOI_from_to[1]+1)
					GOI_info_d=rbind(GOI_info_d, feature_d[1,c("chromosome","from","to", "name")])
					
					print(paste("output file:", out_GOI_info_f))
					
					
					##extract GOI sequence and output .fa file
					plasmid_fa_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".fa")
					plasmid_GOIfa_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".GOI.fa")
					plasmid_stringSet=readDNAStringSet(plasmid_fa_f)
					plasmid_fa=toString(plasmid_stringSet[1])
					GOI_seq=subseq(plasmid_fa, start=GOI_from_to[1], end=GOI_from_to[2])
					out_fa_string=c(paste0(">",master_tb_d[rowi,"plasmid.name"],".GOI"), GOI_seq)
					write(out_fa_string, file=plasmid_GOIfa_f, ncol=1)
				}else{
					print(paste("Error, no enough/correct ITR annotation found for ", feature_f));
				}
			}
			#add additional features (match Out_addi_feature_patt) to GOI_info_d
			ifout_feat=grepl(Out_addi_feature_patt, feature_d$name, ignore.case=T)
			if(any(ifout_feat)){
				GOI_info_d=rbind(GOI_info_d, feature_d[ifout_feat,c("chromosome","from","to", "name")])
			}
			if(nrow(GOI_info_d)>0){
				write.table(GOI_info_d, file=out_GOI_info_f, sep="\t", col.names=T, row.names=F, quote=F)
			}

			feature_d$from_m1=feature_d$from-1
			feature_d$score_bed=1
			feature_d$score_gtf="."
			feature_d$frame="."
			feature_d$source="plasmid"
			feature_d$feature="exon"
			feature_d$group=paste0("gene_id \"",feature_d$name,"\"; transcript_id \"",feature_d$name,"\";")
			feature_d$exon_num=1
			feature_d$gene_id=0
			feature_d$gene_desc=paste(feature_d$chromosome,feature_d$name, sep=".")
			feature_d$gene_Biotype="mRNA" #arbitrary value
			feature_d$size=feature_d$to -feature_d$from_m1
			feature_d$too_small=feature_d$size<50

			#output bed and gtf files:
			out_bed_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".bed")
			out_gtf_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".gtf")
			out_refflat_f=paste0(master_tb_d[rowi,"Folder"],"/",master_tb_d[rowi,"plasmid.name"],".refflat")
			out_bed_d=feature_d[,c("chromosome","from_m1","to", "name", "score_bed","strand")]
			out_gtf_d=feature_d[,c("chromosome","source","feature","from","to", "score_gtf","strand","frame","group")]
			out_refflat_d=feature_d[,c("gene_desc", "gene_desc", "chromosome","strand","from_m1","to","from_m1","to","exon_num","from_m1","to","gene_id",  "gene_desc", "gene_Biotype")]
			#gene_symbol     refseqid        contig  strand  transc_start    transc_end      cds_start       cds_end exon_num        exon_starts     exon_ends       gene_id gene_desc       gene_Biotype
			names(out_refflat_d)=unlist(strsplit("gene_symbol refseqid contig strand transc_start transc_end cds_start cds_end exon_num exon_starts exon_ends gene_id gene_desc gene_Biotype"," "))

			print(paste("output file:", out_bed_f))
			print(paste("output file:", out_gtf_f))
			print(paste("output file:", out_refflat_f))
			write.table(out_bed_d, file=out_bed_f, sep="\t", col.names=F, row.names=F, quote=F)
			write.table(out_gtf_d, file=out_gtf_f, sep="\t", col.names=F, row.names=F, quote=F)
			write.table(out_refflat_d[!feature_d$too_small,], file=out_refflat_f, sep="\t", col.names=T, row.names=F, quote=F)
		}

	}else{
		print (paste("Not found file or already processed:", feature_f))
	}
}

print("Processed files:")
print(processed_files)


##create an .fa and .bed file for all vectors #to create IGV genome file
out_prefix="All_Plasmids_Genome/All_Plasmids_Genome."
system("mkdir All_Plasmids_Genome")
allbed_files=list.files(pattern = "*/*.bed$", recursive=T)
all_fa_files=list.files(pattern = "*/*.fa$", recursive=T)

rm_vectors="All_Plasmids_Genome|lambda|\\.GOI\\.fa$"
all_fa_files=all_fa_files[! grepl(rm_vectors, all_fa_files) ]; all_fa_files
allbed_files=allbed_files[! grepl(rm_vectors, allbed_files) ]; allbed_files
out_bed_f=paste0(out_prefix,"bed")
out_fa_f=paste0(out_prefix,"fa")

merge_bed_cmd=paste(c("cat ",allbed_files, ">", out_bed_f), collapse=" ")
merge_fa_cmd=paste(c("cat ",all_fa_files, ">", out_fa_f), collapse=" ")
system(merge_bed_cmd)
system(merge_fa_cmd)
