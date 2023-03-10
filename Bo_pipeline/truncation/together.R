install.packages("ggplot2",repos = "http://cran.us.r-project.org")
#install.packages("ggbreak",repos = "http://cran.us.r-project.org")
#install.packages("ggthemes",repos = "http://cran.us.r-project.org")
#install.packages("ggpubr",repos = "http://cran.us.r-project.org")
library(ggplot2)
library(ggbreak)
library(ggthemes)
library(ggpubr)
theme_set(theme_pubr())

data1=read.table("/home/Research_Archive/ProcessedData/Nanopore_seq/2022-04-19_AAV_SEQ/truncation/ID23/5cutoff/bin/plot_bin_irreads_only/final.after.bin",header=T,sep="\t")
#pdf("IR_reads.plus.pdf")
reads_coverage_plus=data1[,3]
position=data1[,1]
p1=ggplot(data1,aes(x=position, y=reads_coverage_plus))+geom_line(aes(y=reads_coverage_plus,colour="IR_read (+)"))+scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=seq(0,5000,by=200))+ scale_color_manual( values = c("IR_read (+)" = "red"))+theme(axis.text.x = element_text(angle=45,hjust=1,size=8))+labs(x='',y="reads_coverage") #+theme(axis.title.x = element_text(size = rel(1), angle = 00))

reads_coverage_minus=data1[,4]
p2=ggplot(data1,aes(x=position, y=reads_coverage_minus))+geom_line(aes(y=reads_coverage_minus,colour="IR_read (-)"))+scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=seq(0,5000,by=200))+ scale_color_manual( values = c("IR_read (-)" = "blue"))+theme(axis.text.x = element_text(angle=45,hjust=1,size=8))+labs(x='',y="reads_coverage") 


data2=read.table("/home/Research_Archive/ProcessedData/Nanopore_seq/2022-04-19_AAV_SEQ/truncation/ID23/5cutoff/bin/plot_bin_irreads_only/final.after.bin",header=T,sep="\t")
#reads_coverage_plus=data1[,5]
position2=data2[,1]


reads_terminal_plus1=data2[,5]
p3=ggplot(data2,aes(x=position2, y=reads_terminal_plus1))+geom_line(aes(y=reads_terminal_plus1,colour="%read_terminal (+)"))+scale_y_continuous(limits=c(0,50))+scale_x_continuous(breaks=seq(0,5000,by=200))+ scale_color_manual( values = c("%read_terminal (+)" = "red"))+theme(axis.text.x = element_text(angle=45,hjust=1,size=8))+labs(x='',y="reads_coverage") 


reads_terminal_minus1=data2[,6]
p4=ggplot(data2,aes(x=position2, y=reads_terminal_minus1))+geom_line(aes(y=reads_terminal_minus1,colour="%read_terminal (-)"))+scale_y_continuous(limits=c(0,50))+scale_x_continuous(breaks=seq(0,5000,by=200))+ scale_color_manual( values = c("%read_terminal (-)" = "blue"))+theme(axis.text.x = element_text(angle=45,hjust=1,size=8))+labs(x='',y="reads_coverage") 

#g2 <- ggplotGrob(p2)
data3=read.table("/home/Research_Archive/ProcessedData/Nanopore_seq/2022-04-19_AAV_SEQ/truncation/ID23/5cutoff/bin/plot_bin_irreads_only/final.after.bin",header=T,sep="\t")
#ggplot(data1,aes(x=Start_position,y=Reads_coverage))+geom_line(aes(y="IR_read_+",colour="blue"))+scale_x_continuous(breaks = seq(0,4808,by=500))+scale_y_continuous(breaks = seq(0,5,by=1))
reads_coverage=data3[,2]
position3=data3[,1]
p5=ggplot(data3,aes(x=position3, y=reads_coverage))+geom_line(aes(y=reads_coverage,colour="Number_of_total_reads"))+scale_y_continuous(limits=c(0,200000))+scale_x_continuous(breaks=seq(0,5000,by=200))+ scale_color_manual( values = c("Number_of_total_reads" = "grey"))+theme(axis.text.x = element_text(angle=45,hjust=1,size=8))+labs(x="position",y="reads_coverage") #+theme(axis.title.x = element_text(size = rel(1), angle = 00))

ID23_20220419binironly=ggarrange(p1,p2,p3,p4,p5, ncol=1,nrow=5, heights = c(10,10,10,10,10),widths=c(30,30,30,30,30), align = "v",legend = "right")

ggsave("20220419_ID23_binironly.png",ID23_20220419binironly)
#ID23_scAAV-SCA3-4E10x2_AAV
#multiplot4(p1,p2,p3,p4,p5,cols=1)

dev.off()
