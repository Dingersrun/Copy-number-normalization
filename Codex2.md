#### Running Codex2 in R to detect copy number for each sample
1, Codex2 does not support the mode of using samples as negative control if there is only 1 control/normal sample in `R`  
2, Codex2 performs the best when one has matched tumor and normal/control samples in `bash`
3, performing differential anlaysis 

##### Run Codex2 in `r`
```r
install.packages('devtools')
devtools::install_github("yuchaojiang/CODEX2/package")
library(CODEX2)

library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38

dirPath<-"/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/Bam_gDNA" 
bamFile <- list.files(dirPath, pattern = '*.bam$') #make sure that index files for .bam are also int he same directory
bamdir <- file.path(dirPath, bamFile)
sampname <- substr(bamFile,1,2) #1,2 means taking the first 2 letters in the filename as the samole name
bedFile <- file.path(dirPath, "chr1.50kb.bed")
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX2_chr1")
#bedFile <- file.path(dirPath, "hg38.norm.sizes.50kb.bed")
#"/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/Bam_gDNA/hg38.norm.sizes.50kb.bed"
#bambedObj <- getbambed(bamdir = bamdir,sampname = sampname, bedFile=bedFile,projectname = "Fib_CNVdetection")

bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

gc <- getgc(ref, genome = genome)
mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))  
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
summary(Y[,1:2])

qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 4000),length_thresh = c(29990, 51000), mapp_thresh = 0.9,gc_thresh = c(20, 80)) #if no exons wre excluded at this step, the prprogram will complain

#Excluded 3669 exons due to extreme coverage.
#Excluded 10 exons due to extreme exonic length.
#Excluded 3935 exons due to extreme mappability.
#Excluded 2366 exons due to extreme GC content.
#After taking union, excluded 7531 out of 57509 exons in QC.

Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc

Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)


normObj.null <- normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1, N = N)
# Computing normalization with k = 1 latent factors ...
# Iteration 1	beta diff =0.0403	f(GC) diff =4.13e-10
# Iteration 2	beta diff =0.000386	f(GC) diff =8.37e-14
# Iteration 3	beta diff =1.91e-07	f(GC) diff =7.6e-18
# Stop at Iteration 3.
# AIC1 = 365688415.854
# BIC1 = 365213013.372
# RSS1 = 0
# 	AIC1 = 365696209.6
# 	BIC1 = 365220775.581
# 	RSS1 = 0

Yhat.null <- normObj.null$Yhat
AIC.null <- normObj.null$AIC; BIC.null <- normObj.null$BIC
RSS.null <- normObj.null$RSS

#testing the command
chr <- "chr22" # This can be run for one chromosome or multiple chromosomes
chr.index <- which(seqnames(ref_qc)==chr)
chr1normObj.null <- normalize_null(Y_qc = Y_qc[chr.index,],
                               gc_qc = gc_qc[chr.index],
                               K = 1, N = N)

test22<- segmentCBS(Y_qc[chr.index,], 
                            Yhat.null, optK = which.max(BIC.null),
                            K = 1,
                            sampname_qc = colnames(Y_qc),
                            ref_qc = ranges(ref_qc[chr.index]),
                            chr = chr, lmax = 4000, mode = "fraction")
#this works. try to apply it to chr1:22
# Initialize an empty list to store results
results_list <- list()

# Loop through chromosomes chr1 to chr22
for (i in 1:22) {
    chr <- paste0("chr", i)  # Create chromosome name
    chr.index <- which(seqnames(ref_qc) == chr)
    
    # Run segmentCBS
    result <- segmentCBS(Y_qc[chr.index,], 
                         Yhat.null, 
                         optK = which.max(BIC.null),
                         K = 1,
                         sampname_qc = colnames(Y_qc),
                         ref_qc = ranges(ref_qc[chr.index]),
                         chr = chr, 
                         lmax = 4000, 
                         mode = "fraction")
    
    # Store results in the list with a unique name
    results_list[[chr]] <- result
}
# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)
#noticed that chr1 is corrupted... no idea why
myBSseg<-final_results[-1,]%>%
    filter(sample_name=="BS")%>%
    transmute(chr=chr,start=st_bp,end=ed_bp,CNV=cnv,length=length_kb, likelyhood_ratio=lratio,CN_segment=copy_no) %>%
    mutate( start = ifelse(length < 0, end, start), end = ifelse(length < 0, start, end),length = abs(length))
write.table(myBSseg,"FibBS.CODEX2_NoControlSample.chr2_22_CN_segment.table",sep="\t",col.names=T,row.names=F,quote=F)

myWTseg<-final_results%>%
    filter(sample_name=="WT")%>%
    transmute(chr=chr,start=st_bp,end=ed_bp,CNV=cnv,length=length_kb, likelyhood_ratio=lratio,CN_segment=copy_no) %>%
    mutate( start = ifelse(length < 0, end, start), end = ifelse(length < 0, start, end),length = abs(length))
write.table(myWTseg,"FibWT.CODEX2_NoControlSample.chr2_22_CN_segment.table",sep="\t",col.names=T,row.names=F,quote=F)
save.image("CODEX2.NOControlSample.chr2_22.RData")

#-----------------------------------------
#normObj <- normalize_codex2_ns(Y_qc = Y_qc[chr.index,],
#                               gc_qc = gc_qc[chr.index], 
#                               K = 1, norm_index = 2,
#                               N = N)
#this can not be run with only 1 normal sample...
#---------------------------------------------


#Rerun the above code for chr1 with a modified bedfile path for chr1
bedFile <- file.path(dirPath, "chr1.50kb.bed")
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX2_chr1")

#segmenting the chr1
#test2_chr1<- segmentCBS(Y_qc, Yhat.null, optK = which.max(BIC.null),K = 1,sampname_qc = colnames(Y_qc),ref_qc = ranges(ref_qc),chr = "chr1", lmax = 4000, mode = "fraction")
# this command however does not work.

test3_chr1<- segmentCBS(Y, Yhat.null, optK = which.max(BIC.null),K = 1,sampname_qc = colnames(Y_qc),ref_qc = ranges(ref),chr = "chr1", lmax = 4000, mode = "fraction") #this command also works...
#Segmenting sample 1: BS.
#Segmenting sample 2: WT.

myBSchr1<-test3_chr1%>%
    filter(sample_name=="BS")%>%
    transmute(chr=chr,start=st_bp,end=ed_bp,CNV=cnv,length=length_kb, likelyhood_ratio=lratio,CN_segment=copy_no) %>%
    mutate( start = ifelse(length < 0, end, start), end = ifelse(length < 0, start, end),length = abs(length))
write.table(myBSchr1,"FibBS.CODEX2_NoControlSample.chr1_CN_segment.table",sep="\t",col.names=T,row.names=F,quote=F)

myWTchr1<-test3_chr1%>%
    filter(sample_name=="WT")%>%
    transmute(chr=chr,start=st_bp,end=ed_bp,CNV=cnv,length=length_kb, likelyhood_ratio=lratio,CN_segment=copy_no) %>%
    mutate( start = ifelse(length < 0, end, start), end = ifelse(length < 0, start, end),length = abs(length))
write.table(myWTchr1,"FibWT.CODEX2_NoControlSample.chr1_CN_segment.table",sep="\t",col.names=T,row.names=F,quote=F)
```


##### Extract the DNA segments from both samples and calculate the copy number ratio in `bash`
```
cat FibWT.CODEX2_NoControlSample.chr1_CN_segment.table FibWT.CODEX2_NoControlSample.chr2_22_CN_segment.table|awk -v OFS='\'t '{print $0}'|sort -k1,1 -k2,2n|tail -n+3 |awk -v OFS='\'t '{if($6<40) print $1,$2*1,$3*1,2;else print $1,$2*1,$3*1,$7}'|sort -k1,1 -k2,2n > FibWT.CODEX2_NoControlSample.bed
cat FibBS.CODEX2_NoControlSample.chr1_CN_segment.table FibBS.CODEX2_NoControlSample.chr2_22_CN_segment.table|awk -v OFS='\'t '{print $0}'|sort -k1,1 -k2,2n|tail -n+3 |awk -v OFS='\'t '{if($6<40) print $1,$2*1,$3*1,2;else print $1,$2*1,$3*1,$7}'|sort -k1,1 -k2,2n > FibBS.CODEX2_NoControlSample.bed

#calculating the CNR for DNA segments based on called copy number from DNA segments in each sample
bedtools intersect -a FibWT.CODEX2_NoControlSample.bed -b FibBS.CODEX2_NoControlSample.bed |awk -v OFS='\t' '{print $1,$2,$3}' > WT_BS_segment.intersect.bed
bedtools intersect -a FibWT.CODEX2_NoControlSample.bed -b FibBS.CODEX2_NoControlSample.bed -wa -wb |awk -v OFS='\t' '{print $1"_"$2"_"$3,$4,$5"_"$6"_"$7,$8}' > WT_BS_segment.intersect_original_segment_WT_BS.table

paste WT_BS_segment.intersect.bed WT_BS_segment.intersect_original_segment_WT_BS.table|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$7/$5}' > FibBSvsWT.WT_BS_segment.CNR_segment.bed

bedtools closest -a Fib_CNcorrection/05CN_DESeq/Merge.ATAC_D_K.20Meach_peaks.narrowPeak.bed -b FibWT.CODEX2_NoControlSample.bed -t first -wa -wb|awk '$6!="-1"'|awk -v OFS='\t' '{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}' > ATAC_D_K.peak.WT_CN.COdex2.bed

bedtools closest -a Fib_CNcorrection/05CN_DESeq/Merge.ATAC_D_K.20Meach_peaks.narrowPeak.bed -b FibBS.CODEX2_NoControlSample.bed -t first -wa -wb|awk '$6!="-1"'|awk -v OFS='\t' '{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}' > ATAC_D_K.peak.BS_CN.COdex2.bed

paste WT_BS_segment.intersect.bed WT_BS_segment.intersect_original_segment_WT_BS.table|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$7/$5}' > FibBSvsWT.WT_BS_segment.CNR_segment.bed

#calcualting the scaling factors with the peaks
join -1 4 -2 4  ATAC_D_K.peak.WT_CN.COdex2.bed ATAC_D_K.peak.BS_CN.COdex2.bed|awk -v OFS='\t' '{print $2,$3,$4,$1,$5"|"$6"|"$10"|"$11,$11/$6}'|sort -k1,1 -k2,2n > ATAC_D_K.peak_scaling_factor.Codex2.bed # calculating the scaling factors for each peak CN_BS/WT_CN
```
##### applying copy number correction and performing the differential analysis with DESeq2
```bash
Rscript Fib_CNcorrection/05CN_DESeq/CNnormalization_DESeq2_DA_ATACseq.R Fib_CNcorrection/05CN_DESeq/Fib.sample.info Fib_CNcorrection/05CN_DESeq/ATAC_Fib_D_K_Merge.ATAC_D_K.20Meach.HTseq.withheader.NochrXY.readcount Fib_CNcorrection/04CODEX2/Codex2_NoControl/ATAC_D_K.peak_scaling_factor.Codex2.bed FibBSvsWT.Codex2_DESeq2

```

