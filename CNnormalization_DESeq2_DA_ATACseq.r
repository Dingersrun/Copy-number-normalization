#2025-Jan
#usage: Rscript sampletable_readcount_with_header_path scaling-facot_for_peaks prefix_for_plot
#remember to modify the ref level accordingly before running the script.
#an example: 
#Rscript /fml/chones/data/Marek/NovaSeqS4_1/bcl_delivery/DS_noMis/Fib_CNcorrection/05CN_DESeq/CNnormalization_DESeq2_DA_ATACseq.R /fml/chones/data/Marek/NovaSeqS4_1/bcl_delivery/DS_noMis/Fib_CNcorrection/05CN_DESeq/Fib.sample.info /fml/chones/data/Marek/NovaSeqS4_1/bcl_delivery/DS_noMis/Fib_CNcorrection/05CN_DESeq/ATAC_Fib_D_K_Merge.ATAC_D_K.20Meach.HTseq.withheader.NochrXY.readcount /fml/chones/data/Marek/NovaSeqS4_1/bcl_delivery/DS_noMis/Fib_CNcorrection/05CN_DESeq/ATAC_D_K.peak_scaling_factor.CopywriteR.bed FibBSvsWT.CopyWriter_DESeq2

#silence the screen output when loading R packages, not suppressing other output messages from other commands
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

#load packages
library(tidyverse)
library(DESeq2)

#read the input for R script
args <- commandArgs(trailingOnly = TRUE)

#sample info: args[1] # the contrast level in the the column named as condisiton
#count matrix: args[2]
#Copy number ratio file: args[3], 6 column: "chr","start","end","Peak_name","Segment","CNR"
#name for the plot: args[4]

#name for the plot
plotname <-args[4]
message("prefix of the plots is ", plotname)
#save the screen output into
sink(paste0(plotname,"DESeq2.R.out"))

# reading in the input files
sampletable <- read.table(args[1],sep="\t",header =TRUE)
sampletable$condition<-as.factor(sampletable$condition)

message("reading in sampletable ...")
sampletable

message("reading in count matrix")
mycount <- read.table(args[2], header=TRUE,sep="\t")
head(mycount)
summary(mycount)
#Copy number ratio file
#myCNR<-read.table("ATAC_D_K.peak_scaling_factor.QDNAseq.bed",header=F, sep="\t")
myCNR<-read.table(args[3],header=F, sep="\t")
colnames(myCNR)<-c("chr","start","end","Peak","Seg","CNR")

message("apply the copy number correction")
# -------apply the copy number correction ot the peaks------------
mycount_cor<-left_join(mycount,myCNR[,c(4,6)],by="Peak")%>%na.omit()

#CNR=BS/WT
#CNR>=1, we scale down the readcount in BS sample; when CNR<1, we scale down readcount in WT sample
mycount_cor <- mycount_cor %>%
  mutate(across(starts_with("D"), ~ ifelse(CNR >= 1, round(. /CNR + 0.5), .)),
         across(starts_with("K"), ~ ifelse(CNR < 1, round(. *CNR + 0.5), .)))

#------------runing DESeq2 analysis-----------------
rownames(mycount_cor)<-mycount_cor$Peak
dds<-DESeqDataSetFromMatrix(countData=mycount_cor[,c(2:7)], colData=sampletable, design = ~condition)

#define the reference level
dds$condition <- relevel(dds$condition, ref = "WTFibroblasts")
dds$condition

#filtering the peaks with low read count, the total number of reads from all samples less than 40
keep <- rowSums(counts(dds)) >= 40
dds<-dds[keep,]
#perform normalization (with the effective library size) and differential analysis
message("DA analysis")
dds<-DESeq(dds)

#resultsNames(dds)
head(results(dds, alpha=0.05))
summary(results(dds, alpha=0.05))
summary(results(dds, alpha=0.01))
summary(results(dds, alpha=0.001))

#transform the results
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
resLFCOrdered <- resLFC[order(resLFC$pvalue),]
DEresult<-as.data.frame(resLFCOrdered)
DEresult$Peak<-row.names(DEresult); row.names(DEresult)<-NULL
#---------------save the results-------------
message("save the DA results")
message(paste0("saving results for ",plotname))
write.table(DEresult[,c(6,1:5)],file=paste0(plotname,".resLFCOrdered.table"),quote=F,col.names=T,row.names=F,sep="\t")

#colnames(DEresult)<-c(paste0(myname,"_baseMean"), paste0(myname,"_log2FC"), paste0(myname,"_lfcSE"), paste0(myname,"pval"), paste0(myname,"padj"),"Peak")
#write.table(DEresult,file=paste0("DESeq2",plotname,".resLFCOrdered.table"),quote=F,col.names=T,row.names=F,sep="\t")

message("plotting...")
#-----------produce and save MA and volcano plots------------
#get the range of log2FC
myy<-max(round(abs(DEresult$log2FoldChange)+0.5))

#vocalno
p1<-ggplot(DEresult) +
	geom_point(aes(x=log2FoldChange,y=(-log10(padj)),color=(padj<0.05)),size=1,shape=16,alpha=0.5) +
	xlim(-myy,myy) +theme_classic() + 
	scale_color_manual(values=c("grey","blue")) + 
	theme(axis.text=element_text(size=10,colour = "black"), axis.title=element_text(size=10)) + 
	ggtitle(paste0("Fib BSvsWT, ATACseq \n",plotname)) + 
	labs(x=paste0("CNV_log2FC ",plotname))

#MAplot
p2<-ggplot(DEresult) + 
 geom_point(aes(x=baseMean, y=log2FoldChange, color=(padj<0.05)),size=1,shape=16,alpha=0.5)+ 
	scale_x_log10() + ylim(-myy,myy) +theme_classic() + 
	scale_color_manual(values=c("grey","blue")) + 
	theme(axis.text=element_text(size=10,colour = "black"),axis.title=element_text(size=10)) + 
	ggtitle(paste0("Fib BSvsWT, ATACseq \n",plotname)) + 
	labs(x="Normalized readcount")

pdf(paste0("FibBSvsWT.",plotname,"volcano_MAplot.pdf"),width=4,height=4)
p1
p1+ theme(legend.position="Null")
p2
p2+ theme(legend.position="Null")
dev.off()

width_px <- 6 * 600 / 2.54  # width in pixels
height_px <- 6 * 600 / 2.54  # height in pixels

png(filename = paste0(plotname,"1.png"), width = width_px, height = height_px, res = 600)
p1+ theme(legend.position="Null")
dev.off()  # Close the PNG device

# Save p2
png(filename = paste0(plotname,"2.png"), width = width_px, height = height_px, res = 600)
p2+ theme(legend.position="Null")
dev.off()
sink()

message("save the R enviroment")
save.image(paste0(plotname,".Rdata"))
