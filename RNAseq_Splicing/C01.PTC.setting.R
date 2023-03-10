

setwd("/drive2/wli/analyze/Pipeline/RNAseq_Splicing")


study_set="SplicingCpds_U1vars"
input_root="/drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/"
input_root2="/drive2/wli/analyze/dep_seq/othProject/public_data/"
output_root="/drive2/wli/analyze/dep_seq/othProject/compound_splicing/1Splicing_comb_tb/"
ana_meth="DEXSeq"; fext="/RNAseq_Splicing/05.DEXSeq/geneBlock/combine_d.tbl"
ana_meth="AllExonPSI"; fext="/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl"
file_setting_list=list(
 list(study="201803_U1vars19",species="human",input_f=paste(input_root, "201803_U1vars19","/all_exon/exons.tbl" , sep=""), change_headerFroms="mock",change_headerTos="mock1803",  test_samples=unlist(strsplit("v1CTT v2CAT v3CGT v4CCT v5CTA v6CAA v7CGA v8CCA v9CTC v10CAC v11CGC v12CCC v13CTG v14CGG v15CCG v16GAATgtaaA v17AGAc v19AGAa"," ")), ref_samples=paste(unlist(strsplit("mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock"," ") ),"1803",sep="") )
 ,list(study="2017_U1AGA",species="human",input_f=paste(input_root, "2017_U1AGA","/all_exon.v19/exons.tbl" , sep=""), change_headerFroms="mock",change_headerTos="mock1711",  test_samples=unlist(strsplit("U1GA1ug U1GA4ug U1WT1ug U1WT4ug U1GA1ug.10nM905"," ")), ref_samples=unlist(strsplit("mock1711 mock1711 mock1711 mock1711 mock1711"," ") ))
 ,list(study="201712_U1var",species="human",input_f=paste(input_root, "201712_U1var","/all_exon.v19/exons.tbl" , sep=""), change_headerFroms=c("U1v13","U1v14","mock"), change_headerTos=c("U1v13AT","U1v14AT","mock1712"),  test_samples=unlist(strsplit("U1WT U1v13AT U1v14AT"," ")), ref_samples=unlist(strsplit("mock1712 mock1712 mock1712"," ") ))
 ,list(study="PTC201801_159_MRC5",species="human",input_f=paste(input_root, "PTC201801_159_MRC5/all_exon/exons.tbl" , sep=""), change_headerFroms=c("DMSO"),change_headerTos="DMSO_MRC5",  test_samples=unlist(strsplit("P159.10nM P159.300nM P159.600nM"," ")), ref_samples=unlist(strsplit("DMSO_MRC5 DMSO_MRC5 DMSO_MRC5"," ") ))
 ,list(study="novartis2015",species="human",input_f=paste(input_root, "novartis2015/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("NVS_SM1"," ")), ref_samples=unlist(strsplit("NVS_DMSO"," ") ))
 ,list(study="2016_C3_828",species="human",input_f=paste(input_root, "2016_C3_828/all_exon/exons.tbl" , sep=""), change_headerFroms=c("DMSO"), change_headerTos=c("DMSO_SHSY5Y"),  test_samples=unlist(strsplit("PTC905.24nM PTC905.100nM SMN.C3.300nM SMN.C3.3uM PTC828.500nM PTC023.500nM"," ")), ref_samples=unlist(strsplit("DMSO_SHSY5Y DMSO_SHSY5Y DMSO_SHSY5Y DMSO_SHSY5Y DMSO_SHSY5Y DMSO_SHSY5Y"," ") ))
 ,list(study="2016_RO247_RO067_hs",species="human",input_f=paste(input_root, "2016_RO247_RO067_hs/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL, change_headerTos=NULL,  test_samples=unlist(strsplit("PNN247300.pA PNN2473uM.pA PNN06730.pA PNN0671uM.pA HDF247300.pA HDF2473uM.pA HDF06730.pA HDF0671uM.pA HDF2473uM.rz"," ")), ref_samples=unlist(strsplit("PNNDMSO.pA PNNDMSO.pA PNNDMSO.pA PNNDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.pA HDFDMSO.rz"," ") ))
 ,list(study="NC2018_PK4C9_Roche",species="human",input_f=paste(input_root2, "NC2018_PK4C9_Roche/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl" , sep=""), change_headerFroms=c("DMSO"), change_headerTos=c("DMSO_GM03813C"),  test_samples=unlist(strsplit("PK4C9"," ")), ref_samples=unlist(strsplit("DMSO_GM03813C"," ") ))
)


study_set="FD_Cpds2"
input_root="/drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/"
input_root2="/drive2/wli/analyze/dep_seq/othProject/2018_kinetin_MGH/RNAseq_Splicing/06.exon_PSI/"
output_root="/drive2/wli/analyze/dep_seq/othProject/compound_splicing/1Splicing_comb_tb/"
ana_meth="AllExonPSI"; 
file_setting_list=list(
 list(study="2017_FD_Cpds",species="human",input_f=paste(input_root, "2017_FD_Cpds/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("p5403.1000 p7671.30 p7671.300 p7744.300 p7750.30 p7808.30 p7814.300"," ")), ref_samples=unlist(strsplit("DMSO DMSO DMSO DMSO DMSO DMSO DMSO"," ") ))
 ,list(study="2015_Kinetin",species="human",input_f=paste(input_root, "2015_Kinetin/all_exon/exons.tbl" , sep=""), change_headerFroms=c("FD","CR","CK","CD"), change_headerTos=c("",".R",".K",".D"),  test_samples=unlist(strsplit("C.R30 P.K100 P.R1 P.R30"," ")), ref_samples=unlist(strsplit("C.D P.D P.D P.D"," ") ))
 ,list(study="2018_kinetin_MGH",species="human",input_f=paste(input_root2, "all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("kin C15477"," ")), ref_samples=unlist(strsplit("ctrl ctrl"," ") ))
)


study_set="U1var_w905_wo905"
input_root="/drive2/wli/analyze/AS/5AS_deepSeq/04exon_PSI/"
output_root="/drive2/wli/analyze/dep_seq/othProject/compound_splicing/1Splicing_comb_tb/"
ana_meth="AllExonPSI"; fext="/RNAseq_Splicing/06.exon_PSI/all_exon/exons.tbl"
file_setting_list=list(
 list(study="201806_U1v_w905",species="human",input_f=paste("/drive2/wli/analyze/dep_seq/othProject/201806_U1v_w905/RNAseq_Splicing/06.exon_PSI","/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("mock.905 WT.905 v1CTT.905 v2CAT.905 v3CGT.905 v4CCT.905 v5CTA.905 v6CAA.905 v7CGA.905 v8CCA.905 v9CTC.905 v10CAC.905 v11CGC.905 v12CCC.905 v13CTG.905 v14CGG.905 v15CCG.905"," ")), ref_samples=unlist(strsplit("mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO"," ") ) )
 ,list(study="201803_U1vars19",species="human",input_f=paste(input_root, "201803_U1vars19","/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("v1CTT v2CAT v3CGT v4CCT v5CTA v6CAA v7CGA v8CCA v9CTC v10CAC v11CGC v12CCC v13CTG v14CGG v15CCG"," ")), ref_samples=unlist(strsplit("mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock"," ") ) )
)
##consolidated exon table recalculated PSI etc
file_setting_list=list(
 list(study="201806_U1v_w905",species="human",input_f=paste("/drive2/wli/analyze/dep_seq/othProject/201806_U1v_w905/RNAseq_Splicing/06.exon_PSI/consolidated.U1var_w905_wo905/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("mock.905 WT.905 v1CTT.905 v2CAT.905 v3CGT.905 v4CCT.905 v5CTA.905 v6CAA.905 v7CGA.905 v8CCA.905 v9CTC.905 v10CAC.905 v11CGC.905 v12CCC.905 v13CTG.905 v14CGG.905 v15CCG.905"," ")), ref_samples=unlist(strsplit("mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO mock.DMSO"," ") ) )
 ,list(study="201803_U1vars19",species="human",input_f=paste(input_root, "201803_U1vars19","/consolidated.U1var_w905_wo905/all_exon/exons.tbl" , sep=""), change_headerFroms=NULL,change_headerTos=NULL,  test_samples=unlist(strsplit("v1CTT v2CAT v3CGT v4CCT v5CTA v6CAA v7CGA v8CCA v9CTC v10CAC v11CGC v12CCC v13CTG v14CGG v15CCG"," ")), ref_samples=unlist(strsplit("mock mock mock mock mock mock mock mock mock mock mock mock mock mock mock"," ") ) )
)

##for U1var_w905_wo905 select exons activated by 905+U1var, but not by 905 only or U1v only
#set P_cut_top=0.001; change_cut=40
use_change_names
ReguType_names
comb_d$if_905_not_up=comb_d[,use_change_names[1]] <10; table(comb_d$if_905_not_up)
if_U1v_w905_up=!is.na(comb_d[,use_change_names[3:17]]) & !is.na(comb_d[,use_change_names[18:32]]) & comb_d[,ReguType_names[3:17]] == "UP" & ( comb_d[,use_change_names[3:17]]-comb_d[,use_change_names[18:32]]>20 ) #exons activated only when both 905 and U1v was added
if_U1v_w905_up=if_U1v_w905_up & comb_d$if_905_not_up
apply(if_U1v_w905_up,2,table)
add_headers=sub("ReguType_","IfFurtherUp_",ReguType_names[3:17]); add_headers
comb_d[,add_headers]=if_U1v_w905_up


##combine with dominant disease information
OMIM_genemapf="/drive2/wli/analyze/summary/GeneSet/OMIM/genemap2.txt"
OMIM_genemapd=read.table(OMIM_genemapf, sep="\t", header=T, comment.char="#", quote="", stringsAsFactors=F)
OMIM_genemapd$if_dominant=grepl("dominant",OMIM_genemapd$Phenotypes, ignore.case=T); sum(OMIM_genemapd$if_dominant)
OMIM_genemapd2=OMIM_genemapd[OMIM_genemapd$if_dominant,]

if_up = rowSums(!is.na(comb_d[,ReguType_names]) & comb_d[,ReguType_names] == "UP")>0; sum(if_up)
if_iExon=comb_d$region_ano=="exon" & !(comb_d$startSS_supp %in% c("Refseq","KnG_Ens"))  & !(comb_d$endSS_supp %in% c("Refseq","KnG_Ens")) ; sum(if_iExon)
sum(if_up & if_iExon)
out_d=comb_d[if_up & if_iExon, ]; nrow(out_d)
out_d$dominant_disease=OMIM_genemapd2$Phenotypes[match(out_d$gene_id, OMIM_genemapd2$Entrez.Gene.ID, incomparables=F, nomatch="")]
out_d$dominant_disease[is.na(out_d$gene_id)| is.na(out_d$dominant_disease)]=""

out_f=paste(output_root, study_set,"/Gblock.Combine.",ana_meth,".P",P_cut_top,".Ch",round(change_cut,2),".iExon.UP.DominantDisease.txt" ,sep="");
write.table(out_d, file=out_f, col.names=T,row.names=F, sep="\t", quote=F)
