#### Running CopywriteR in R to detect the copy number for each sample
1, `CopywriteR` is the `R` package for copy number detection used in `cobra` containerized workflow.
Here we use it for benchmarking our proposed pipeline. 
2, `CopywriteR` by default used the data from 1000 genome project as the control
3, performing differential analysis with `DESeq2`

##### Running CopyWriter in `R`
```
samples=c("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/03HMcan_diff//04BSFib_D_gDNAbam/BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam","/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/03HMcan_diff/05WTFib_K_gDNAbam/WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam")
controls=rep("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/03HMcan_diff/05WTFib_K_gDNAbam/WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam",2)  
sample.control <- data.frame(samples=samples,controls=controls)
bp.param <- SnowParam(workers = 10, type = "SOCK")
CopywriteR(sample.control = sample.control, destination.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter"),reference.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter/hg38_50kb_chr/"),bp.param = bp.param)
#this step takes long
plotCNA(destination.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter"))
```
##### Extract the copy number from the output 
```
#The first one is interesting to us. as it is the CNR...
log2.BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam.vs.log2.WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam

#the segment data is saved in the output segment.Rdata
#load it into R and check the strucuture of the object
#IMPORTANT,   the segment copy number ratio is log2 transformed!!
load("segment.Rdata")
test<-segment.CNA.object$output
myseg<-test%>%filter(ID=="log2.BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam.vs.log2.WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam")
myseg<-myseg[,-1]%>%transmute(chr=paste0("chr",chrom),start=round(0.5+loc.start),end=rounf(loc.end),num.mark=num.mark,log2segment=seg.mean, CNR_segment=2^seg.mean)
write.table(myseg,"CopywriteR.FibBSvsWT.CNR_segment.table",col.names=T,row.names=F,quote=F,sep="\t")
save.image("segment_BSvsWT_log2CNR.Rdata")

interval<-segment.CNA.object$data
mybin<-interval[,c(1,2,3)]%>%transmute(chr=paste0("chr",chrom),start=round(maploc+0.5),ends=round(maploc+50000),log2CN=log2.BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam.vs.log2.WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam, CN_interval=2^log2.BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam.vs.log2.WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam)

write.table(mybin,"CopywriteR.FibBSvsWT.CNR_50kb_bins.table",col.names=T,row.names=F,quote=F,sep="\t")
save.image("Intervals_BSvsWT_log2CNR.Rdata")
```
#In bash to assign the scaling factors to peaks
```
cat CopywriteR.FibBSvsWT.CNR_segment.table|awk -v OFS='\t' '{print $1,$2,$3,$6}'|sort -k1,1 -k2,2n > CopywriteR.FibBSvsWT.CNR_segment.bed
```
Then perform copy number normalization and differential analysis
