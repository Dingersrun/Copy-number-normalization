#plotting the genome-wide copy number variation
load("/fml/chones/data/Dingwen/G4_BlmSyn/Novaseq/WT_BS_ChIP/08CallingCNV/ATACseq_gDNA/02bin50kb/D_K_CNV.before_segment_filtering/Fib_BSvsWT.CNR.replor2024.Rdata")
library(tidyverse)
#read in the cumulative start and end of each chromosome in the order of 1,2,3,4...X,Y and chr M
accum <- read.table("/fml/chones/data/Dingwen/G4_BlmSyn/Novaseq/WT_BS_ChIP/03BamCoverage/hg38.chrM_end.accum",sep="\t",header = F)
colnames(accum)<-c("chr","start","end","cumstart","cumend")

#read in the centromeric regions
centrm <- read.table("/fml/chones/data/Dingwen/G4_BlmSyn/Novaseq/WT_BS_ChIP/03BamCoverage/hg38.centromers.UCSC.bed",sep="\t",header = F)
colnames(centrm)<-c("chr","cent_start","cent_end","cent_name") 
#calculate the cumulative coordinates for centromeric regions
centrm<-left_join(centrm,accum,by="chr")
centrm<-centrm%>%transmute(chr=chr,start=cent_start,end=cent_end,cent_cumstart=cent_start+cumstart-1,Cumend=cent_end+cumstart-1)

#read in the CNV file for each 50kb windows, later polotted as yellow dots
mybin<-read.table("/fml/chones/data/Dingwen/G4_BlmSyn/Novaseq/WT_BS_ChIP/08CallingCNV/ATACseq_gDNA/02bin50kb/BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.cnr",sep="\t",header=T)
colnames(mybin)[1]<-"chr"
mybin<-left_join(mybin, accum, by = "chr")
mybin<-mybin%>%transmute(chr=chr,start=start.x,end=end.x,log2=log2,bin_cumstart=start.x+cumstart-1,bin_cumend=end.x+cumstart-1)

#read in the DNA segment files, later plotted as lines
myseg<-read.table("BS_D_ref2K_segment.CopyNumberRatio.bed",sep="\t",header=F)
colnames(myseg)<-c("chr","seg_start","seg_end","CNR")
myseg<-left_join(myseg,accum,by="chr")
myseg<-myseg%>%transmute(chr=chr,start=seg_start,end=seg_end,log2=log2,seg_umstart=seg_start+cumstart-1,seg_cumend=seg_end+cumstart-1)

#extract the intervals
#myseg<-myseg[,c(1,2,3,5)]
#colnames(myseg)[1]<-"chr";
#myseg<-left_join(myseg,accum,by="chr")
#myseg<-myseg%>%transmute(chr=chr,start=start.x,end=end.x,log2=log2,Cumstart=start.x+cumstart-1,Cumend=end.x+cumstart-1)

#plot(mybin$Cumstart,2^mybin$log2, col="grey",pch=20,cex=0.6,ylim=c(-10,10))
#segments(myseg[,5],2^myseg[,4],myseg[,6],2^myseg[,4],col="blue",lty=1,lwd=3)
#rect(centrm$Cumstart,-10,centrm$Cumend,-8,color = "orange")

#plot with ggplot
p_genome<-ggplot(data=mybin)+theme_classic() + 
    #background to  highlight the interval for each chromosome
    geom_rect(data=accum,aes(xmin=cumstart,xmax=cumend,ymin=min(2^(mybin$log2))+0.05,ymax=max(2^(mybin$log2))-0.05),fill=c(rep(c("grey90","white"),12),"white")) +
    #dots of each 50kb bins
    geom_point(aes(x=round((bin_cumstart+bin_cumend)/2),y=2^log2),color="orange",alpha=0.5,size=0.5,shape=20) + 
    #lines of CNV
    geom_segment(data=myseg,aes(x=seg_cumstart,y=CNR,xend=seg_cumend,yend=CNR),color="blue",lwd=1) + 
    #centromeres
    geom_vline(xintercept=c(centrm$cent_cumstart,centrm$Cumend),color="grey",lty=2,alpha=0.5) + 
    labs(x="Position",y="copy number ratio") + 
    scale_x_continuous(breaks=c((accum$cumstart+accum$cumend)/2)[1:24],labels=accum$chr[1:24])
 p_genome + theme(axis.text.x=element_text(angle=45),axis.text = element_text(size = 14,color = "black"),axis.title = element_text(size = 14))
p_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV")
p_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV")+coord_cartesian(ylim=c(0,3))

p_genome2<-ggplot(data=mybin%>%filter(chr!="chrY"))+theme_classic() + 
    #background to  highlight the interval for each chromosome
    geom_rect(data=accum%>%filter(chr!="chrY"),aes(xmin=cumstart,xmax=cumend,ymin=min(2^(mybin$log2))+0.05,ymax=max(2^(mybin$log2))-0.05),fill=c(rep(c("grey90","white"),12))) +
    #dots of each 50kb bins
    geom_point(aes(x=round((bin_cumstart+bin_cumend)/2),y=2^log2),color="orange",alpha=0.5,size=0.5,shape=20) + 
    #lines of CNV
    geom_segment(data=myseg%>%filter(chr!="chrY"),aes(x=seg_cumstart,y=CNR,xend=seg_cumend,yend=CNR),color="blue",lwd=1) + 
    #centromeres
    geom_vline(xintercept=c((centrm%>%filter(chr!="chrY"))$cent_cumstart,(centrm%>%filter(chr!="chrY"))$Cumend),color="grey",lty=2,alpha=0.5) + 
    labs(x="Position (Mb)",y="copy number ratio") + 
    scale_x_continuous(breaks=c((accum$cumstart+accum$cumend)/2)[1:23],labels=accum$chr[1:23])

p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV")
p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV")+coord_cartesian(ylim=c(0,3.6))
p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV")+coord_cartesian(ylim=c(0,3))


#plot the CNV for each chromosome

#plot stuffs for each of hg38 chromosmoe

hg38chr<-c( paste0("chr",c(seq(1:22),"X","Y")))
list_plot<-list()
index <- 1
for (mychr in hg38chr) {

p<-ggplot(data = mybin %>% filter(chr == mychr)) + theme_classic() + 
  labs(x="Position (Mb)",y="Copy number ratio",title=paste0("Fib-BS vs Fib-WT:",mychr)) + 
  geom_rect(data = centrm%>% filter(chr == mychr),
                aes(
                  xmin = start/1000000,
                  xmax = end/1000000,
                  ymin = 0.05,
                  ymax = round(2^max((mybin %>% filter(chr == mychr))$log2)+0.5)-0.5),
                fill = "grey90") +
  geom_point(
    aes(x = start/1000000, y = 2 ^ log2),
    color = "orange",
    alpha = 0.5,
    size = 0.5,
    shape = 16)  + 
  geom_segment(
    data = myseg%>% filter(chr == mychr),
    aes(x = start/1000000,y = CNR,xend = end/1000000,yend = CNR),
    color = "blue",
    lwd = 1
  ) + ylim(0,round(2^max((mybin %>% filter(chr == mychr))$log2)+0.5)) +xlim(0,((accum%>% filter(chr == mychr))[,3])/1000000)

list_plot[[index]] <- p
index=index+1

pdf(paste0("CNR_FibBSvsWT_replot",mychr,".pdf"),width=6,height=2.5)
print(p)
dev.off()
}

pdf("CNR_FibBSvsWT_replot_each_chromosome.pdf",width=6,height=2.5)
for (index in 1:24) {
print(list_plot[[index]])
}
dev.off()

pdf("Fib_bin50kb_chr12modified_2024_replot.pdf",width=6,height=3)
p_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")

p_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")+coord_cartesian(ylim=c(0,3.6))

p_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")+coord_cartesian(ylim=c(0,3))

p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")
p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")+coord_cartesian(ylim=c(0,3.6))

p_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Copy number ratio")+coord_cartesian(ylim=c(0,3))
dev.off()

p2_genome<-ggplot(data=mybin)+theme_classic() + 
    #background to  highlight the interval for each chromosome
    geom_rect(data=accum,aes(xmin=cumstart,xmax=cumend,ymin=min(mybin$log2)+0.05,ymax=max(mybin$log2)-0.05),fill=c(rep(c("grey90","white"),12),"white")) +
    #dots of each 50kb bins
    geom_point(aes(x=round((bin_cumstart+bin_cumend)/2),y=log2),color="orange",alpha=0.5,size=0.5,shape=20) + 
    #lines of CNV
    geom_segment(data=myseg,aes(x=seg_cumstart,y=log2(CNR),xend=seg_cumend,yend=log2(CNR)),color="blue",lwd=1) + 
    #centromeres
    geom_vline(xintercept=c(centrm$cent_cumstart,centrm$Cumend),color="grey",lty=2,alpha=0.5) + 
    labs(x="Position",y="log2(copy number ratio)") + 
    scale_x_continuous(breaks=c((accum$cumstart+accum$cumend)/2)[1:24],labels=accum$chr[1:24])

p2_genome2<-ggplot(data=mybin%>%filter(chr!="chrY"))+theme_classic() + 
    #background to  highlight the interval for each chromosome
    geom_rect(data=accum%>%filter(chr!="chrY"),aes(xmin=cumstart,xmax=cumend,ymin=min(mybin$log2)+0.05,ymax=max(mybin$log2)-0.05),fill=c(rep(c("grey90","white"),12))) +
    #dots of each 50kb bins
    geom_point(aes(x=round((bin_cumstart+bin_cumend)/2),y=log2),color="orange",alpha=0.5,size=0.5,shape=20) + 
    #lines of CNV
    geom_segment(data=myseg%>%filter(chr!="chrY"),aes(x=seg_cumstart,y=log2(CNR),xend=seg_cumend,yend=log2(CNR)),color="blue",lwd=1) + 
    #centromeres
    geom_vline(xintercept=c((centrm%>%filter(chr!="chrY"))$cent_cumstart,(centrm%>%filter(chr!="chrY"))$Cumend),color="grey",lty=2,alpha=0.5) + 
    labs(x="Position (Mb)",y="copy number ratio") + 
    scale_x_continuous(breaks=c((accum$cumstart+accum$cumend)/2)[1:23],labels=accum$chr[1:23])

pdf("Fib_bin50kb_chr12modified_log2CNR.pdf",width=6,height=3)
p2_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Log2(copy number ratio)")

p2_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Log2(copy number ratio)")+coord_cartesian(ylim=c(-2,2))

p2_genome + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Log2(copy number ratio)")+coord_cartesian(ylim=c(-3,3))

p2_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Log2(copy number ratio)")
p2_genome2 + theme(axis.text.x=element_text(angle=90),axis.text = element_text(size = 8,color = "black"),axis.title = element_text(size = 8))+labs(title="Fib BS vsWT, whole genome CNV",y="Log2(copy number ratio)")+coord_cartesian(ylim=c(-2,2))

dev.off()