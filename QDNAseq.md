Copy number normalization using QDNAseq to chracterize CNR

### Identify CNV in D and K cell lines; using Fib-BS (D) as an example
```
#change the bamfile name for each sampel accordingly
#identify copy number variation with QDNAseq. By default it uses the samples from 1000 genome proejct as the control.

#remotes::install_github("asntech/QDNAseq.hg38@main")
library(QDNAseq.hg38) #yesh, it works, I dont have to build the bin annotation file myself
library(QDNAseq)
library(tidyverse)

bins <- getBinAnnotations(binSize=50, genome="hg38")
readCounts <- binReadCounts(bins, bamfiles="BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam")
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
#plot(copyNumbersSmooth) —> plot saved

exportBins(copyNumbersSmooth, file="Fib_BS_D_50kb.txt")
exportBins(copyNumbersSmooth, file="Fib_BS_D_50.bed", format="bed")

#segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersCalled <- callBins(copyNumbersSegmented)

# Assuming 'copyNumbersCalled' is your QDNAseqCopyNumbers object and the data extraction has been done properly

# Combine all data into a single data frame
mydata <- data.frame(
  copynumber = copyNumbersCalled@assayData$copynumber,
  calls = copyNumbersCalled@assayData$calls,
  probamp = copyNumbersCalled@assayData$probamp,
  probdloss = copyNumbersCalled@assayData$probdloss,
  probgain = copyNumbersCalled@assayData$probgain,
  probloss = copyNumbersCalled@assayData$probloss,
  probnorm = copyNumbersCalled@assayData$probnorm,
  segmented = copyNumbersCalled@assayData$segmented,
  row.names = rownames(copyNumbersCalled@assayData$copynumber)  # assuming all elements share the same row names
)
# Setting column names manually to ensure they are correct
colnames(mydata) <- c("copynumber", "calls", "probamp", "probdloss", "probgain", "probloss", "probnorm", "segmented")

mydata$bed<-rownames(mydata)

# Display the corrected data frame with column names and make it in the bed format
head(mydata)

mydata<- mydata%>% separate(bed, into = c("chr", "interval"), sep = ":") %>% separate(interval, into = c("start", "end"), sep = "-", convert = TRUE)
mydata<-mydata[,c(9,10,11,1:8)]
head(mydata)
rownames(mydata)<-NULL
mydata$chr<-paste0("chr",mydata$chr)
write.table(mydata,"Fib_BS_CN.QDNAseq.50kbinterval.table",row.names=F,col.names=T,quote=F,sep="\t")

#merge adjacent intervals to recover the bed intervals for segments
#testing the command
myseg<-mydata%>% na.omit() %>% group_by(chr, segmented) %>% 
summarize(SegStart = min(start), SegEnd= max(end), .groups = "drop") %>% 
arrange(chr, SegStart)

myseg<-myseg[,c(1,3,4,2)]
myseg$chr<-paste0("chr",myseg$chr)
write.table(myseg,"Fib_BS_CN.QDNAseq.Segment.table",row.names=F,col.names=T,quote=F,sep="\t")
```

### get the copy number from each sample and calculate the copy number ratio of BS/WT
#note that copy number is not log2 transformed.
```r
library(tidyverse)
myBS<-read.table("Fib_BS_CN.QDNAseq.50kbinterval.table",header=T, sep="\t")
myBS<-myBS[,c(1,2,3,4,5,11)]
colnames(myBS)[4:6]<-paste0("BS_",colnames(myBS)[4:6])

myWT<-read.table("Fib_WT_CN.QDNAseq.50kbinterval.table",header=T, sep="\t")
myWT<-myWT[,c(1,2,3,4,5,11)]
colnames(myWT)[4:6]<-paste0("WT_",colnames(myWT)[4:6])

myCNR<-cbind(myWT,myBS[,c(4,5,6)])%>%mutate(CNR_interval=BS_copynumber/WT_copynumber, CNR_segment=BS_segmented/WT_segmented)

write.table(na.omit(myCNR)[,c(1,2,3,6,9,11)],"Fib_BSvsWT.QDNAseq.CNR_interval_withheader.txt",col.names=T,row.names=F,quote=F,sep="\t")

save.image("Fib_BSvsWT.QDNAseq.CNR_segment.RData") 
```
### final CNR
collapse the adjacent intervals with the same segment values and extract the segment information.
using the Copy number ratio calculated from log2 converted CN for copy number normalization

```bash
cat Fib_BSvsWT.QDNAseq.CNR_interval_withheader.txt|awk 'BEGIN {FS=OFS="\t"} {if (NR == 1) {current_chr = $1; current_start = $2; current_end = $3; current_WT = $4; current_BS = $5; current_CNR = $6; next} if ($1 == current_chr && $4 == current_WT && $5 == current_BS && $6 == current_CNR) {current_end = $3} else {print current_chr, current_start, current_end, current_WT, current_BS, current_CNR; current_chr = $1; current_start = $2; current_end = $3; current_WT = $4; current_BS = $5; current_CNR = $6}} END {print current_chr, current_start, current_end, current_WT, current_BS, current_CNR}' > Fib_BSvsWT.QDNAseq.CNR_segment_withheader.txt

tail -n+2 Fib_BSvsWT.QDNAseq.CNR_segment_withheader.txt|awk -v OFS='\t' '{print $1,$2,$3,$4}' >Fib_BSvsWT.QDNAseq.CNR_segment.bed

/fml/chones/data/Marek/NovaSeqS4_1/bcl_delivery/DS_noMis/Fib_CNcorrection/02QDNAseq/Fib_BSvsWT.QDNAseq.CNR_segment.bed
```
Then assign the scaling factors to peaks, perform copy number normalization and differential analysis with DEseq2 using https://github.com/Dingersrun/Copy-number-normalization/edit/main/QDNAseq.md#:~:text=CNnormalization_DESeq2_DA_ATACseq


