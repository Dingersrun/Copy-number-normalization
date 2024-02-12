# Differential analysis of G4 ChIpseq data with/out CNV correction via DiffBind

## pre-peak calling
read trimming, alignment and filtering are the same as decribed in ATACseq data processing

## peak calling and IDR filtering
### preparing the bamfiles

```bash
# for each biological replicate of a sample, there are 3 technical replicates for IP samples and one Input sample
WT_rep1.IP1.bam
WT_rep1.IP2.bam
WT_rep1.IP3.bam
WT_rep1.Input.bam

WT_rep2.IP1.bam
WT_rep2.IP2.bam
WT_rep2.IP3.bam
WT_rep2.Input.bam

BS_rep1.IP1.bam
BS_rep1.IP2.bam
BS_rep1.IP3.bam
BS_rep1.Input.bam

BS_rep2.IP1.bam
BS_rep2.IP2.bam
BS_rep2.IP3.bam
BS_rep2.Input.bam

# subsample 30Mil reads from each IP samples and merge them
WT_rep1.merge.bam
WT_rep2.merge.bam
BS_rep1.merge.bam
BS_rep2.merge.bam

# merge the input samples as well so that during peak calling, macs2 will not downsample the IP samples as there are less reads in Input smaples
WT_input.bam
BS_input.bam

```
### peak calling with macs2 for files merged from technical replicates

```bash
#!/bin/bash
#usage: scrit.sh IP_file Input_file
file=$1
fbname=$(basename $file .bam)
echo $fbname

#the control file for peak calling
file2=$2
echo $file2
macs2 callpeak --tempdir /tmp/dsu/ -g hs \
	-t $file -c $file2 \
	--outdir /tmp/dsu/ \
	-n "$fbname" --seed 11 --cutoff-analysis --keep-dup all \
	2> "$fbname".callpeak.log
```
### IDR filtering on the peaks called from 2 biological replicates

```bash
#!/bin/bash
file1=$1;
echo $file1
fbname1=$(basename $file1 _peaks.narrowPeak )
echo $fbname1

file2=$2;
echo $file2
fbname2=$(basename $file2 _peaks.narrowPeak )
echo $fbname2
#sort the narrowPeak files by the Signal value, which is the default setting    
    sort -k7,7nr $file1 >"$file1".Signal
    sort -k7,7nr $file2 >"$file2".Signal
    idr --samples "$file1".Signal "$file2".Signal --rank signal.value --input-file-type narrowPeak --log-output-file IDR_"$fbname1"_"$fbname2".Signal.idr_log  --plot --output-file IDR_"$fbname1"_"$fbname2".Signal.idr
```
```bash
#filter the peaks with the threshold of IDR<0.05
# remove peaks not on cannonical chromosomes
# remove peaks on sex chromosomes
cat IDR_WT.Signal.idr| awk '{if($5 >= 540) print $0}' |awk '$1!~/_/' |grep -v chrX|grep -v chrY |awk '{print $1"\t"$2"\t"$3}'>IDR_WT.Signal_IDR0.05.bed

cat IDR_BS.Signal.idr| awk '{if($5 >= 540) print $0}' |awk '$1!~/_/' |grep -v chrX|grep -v chrY |awk '{print $1"\t"$2"\t"$3}'>IDR_WT.Signal_IDR0.05.bed
```
## making greylist for later differential analysis
### calling peaks with input sample
```bash
#calling narrowPeaks with -f BAM
macs2 callpeak --tempdir /tmp/dsu/ -g hs \
	-t $Input_file \
	--outdir /tmp/dsu/ \
	-n "$fbname"_BAM \
	--seed 25 -B --SPMR --keep-dup all

macs2 callpeak --tempdir /tmp/dsu/ -g hs \
	-t $Input_file \
	--outdir /tmp/dsu/ \
	-n "$fbname"_BAM_broad \
	--seed 25 -B --SPMR --keep-dup all --broad

macs2 callpeak --tempdir /tmp/dsu/ -g hs \
	-t $Input_file \
	--outdir /tmp/dsu/ \
	-n "$fbname"_BAMPE \
	--seed 25 -B --SPMR --keep-dup all -f BAMPE

# for each input sample, merge the all the above peaks called from each command
cat WT_input_peaks.narrowPeak WT_input*_peaks.broadPeak |awk '$9>=3'| awk '{print $1"\t"$2"\t"$3"\t"$4}' |sort -k1,1 -k2,2n |bedtools merge -i stdin >WT_input.merge_broad_narrow_peaksq0.001.bed

cat BS_input_peaks.narrowPeak BS_input*_peaks.broadPeak |awk '$9>=3'| awk '{print $1"\t"$2"\t"$3"\t"$4}' |sort -k1,1 -k2,2n |bedtools merge -i stdin >BS_input.merge_broad_narrow_peaksq0.001.bed
```

### generating greylist via DiffBind
**here it is not the differential anlaysis. It is to generate a list of regions that have high coverage in the input samples**

- folowing DiffBind vegnette, load data into DiffBind and proceed until the step of differntial anlaysis

- extract and save the greylist in bed format as `DiffBind_greylist.bed`

### preparing the final greylist for Differential analysis
```bash
cat WT_input.merge_broad_narrow_peaksq0.001.bed BS_input.merge_broad_narrow_peaksq0.001.bed DiffBind_greylist.bed|sort -k1,1 -k2,2n |bedtools merge -i stdin |cut -f1-3 >Greylist.bed
```
## Differential analysis in DiffBind

### prepare metadata for samples in the differential analysis
note that here technical replicates of the same biological replicates are treated as individual samples instead of being collasped as one sample. Metadata is stored as ``
### Differential analysis in R with DiffBind

#### read in the data and count reads in peaks
```r
#an example: Rscript ../DiffBind.R Sample.info BSvsWT Greylist.bed

#silence the screen output when loading R packages, not suppressing other output messages from other commands
#suppressMessages(library(DiffBind))
#load packages
library(parallel)
library(DiffBind)
library(BiocParallel)
register(MulticoreParam(workers = 6)) 

#read the input for R script
args <- commandArgs(trailingOnly = TRUE)

#sample metadata
message("sample sheet is", args[1])
sampleinfo<-args[1]

sample <-read.table(sampleinfo,sep="\t",header=T)

#name of plots
message("mysample name is ",args[2])
mysample<-args[2]

#greylist
#greylist is not essential. Greylist can be generated internally by DiffBind when gDNA data is included

mygreylist<-read.table(args[3],header=FALSE,sep="\t")
colnames(mygreylist) <- c("chr","start","end")
mygreylist<-GRanges(mygreylist)

sink(paste0(mysample,"Diffbind.R.out"))
#Read in the samples and their bamfiles
#strangely R encounter erros when I Run the script.
#when I enter R and read in the sampleinfo, things are fine
myG4 <- dba(sampleSheet=sample)
myG4

#remove regions in the blacklist,#remove regions in the greylist
myG4 <- dba.blacklist(myG4,blacklist=DBA_BLACKLIST_HG38,greylist=mygreylist)

#count the reads in peaks!
#crtitial to set Summits=FALSE!!!!
message("\ncount reads in peaks...")
myG4<-dba.count(myG4,summits=FALSE)
```

#### optional copy number correction:

- export the intervals of merged peaks in myG4 object 
- assign the scaling factor to those intervals the same way as it in ATACseq data
- output: `scaling_factors.bed` with 4 columns, chr start end CNR

```bash
awk -v OFS='\t' '{print $1,$2,$3,$6}' scaling_factors.bed|awk -v OFS='\t' '{if ($4 >1) print $1,$2,$3,$4,1;else print $1,$2,$3,1,1/$4}'|sort -k1,1 -k2,2n > myscale
```
```r
#read myscale into R
#the bed intervals are in the same order as the count matrix stored in myG4 object
#modify the readcount

#myG4$peaks[[1]] --> the first sample
#[,c(4,5,6,7)] --> 4,5,6,7 column for the first sample, which are the read count infomation in peaks

#in my case the 1-5 samples are BS and 6-11 samples are WT. Modify the count as the following
myG4$peaks[[1]][,c(4,5,6,7,8)]<-myG4$peaks[[1]][,c(4,5,6,7,8)]/myscale[,4]
myG4$peaks[[2]][,c(4,5,6,7,8)]<-myG4$peaks[[2]][,c(4,5,6,7,8)]/myscale[,4]
myG4$peaks[[3]][,c(4,5,6,7,8)]<-myG4$peaks[[3]][,c(4,5,6,7,8)]/myscale[,4]
myG4$peaks[[4]][,c(4,5,6,7,8)]<-myG4$peaks[[4]][,c(4,5,6,7,8)]/myscale[,4]
myG4$peaks[[5]][,c(4,5,6,7,8)]<-myG4$peaks[[5]][,c(4,5,6,7,8)]/myscale[,4]

#Correcting the read count for WT
myG4$peaks[[6]][,c(4,5,6,7,8)]<-myG4$peaks[[6]][,c(4,5,6,7,8)]/myscale[,5]
myG4$peaks[[7]][,c(4,5,6,7,8)]<-myG4$peaks[[7]][,c(4,5,6,7,8)]/myscale[,5]
myG4$peaks[[8]][,c(4,5,6,7,8)]<-myG4$peaks[[8]][,c(4,5,6,7,8)]/myscale[,5]
myG4$peaks[[9]][,c(4,5,6,7,8)]<-myG4$peaks[[9]][,c(4,5,6,7,8)]/myscale[,5]
myG4$peaks[[10]][,c(4,5,6,7,8)]<-myG4$peaks[[10]][,c(4,5,6,7,8)]/myscale[,5]
myG4$peaks[[11]][,c(4,5,6,7,8)]<-myG4$peaks[[11]][,c(4,5,6,7,8)]/myscale[,5]

#Round the read count to integer numbers
myG4$peaks[[1]][,6]<-round(myG4$peaks[[1]][,6])
myG4$peaks[[2]][,6]<-round(myG4$peaks[[2]][,6])
myG4$peaks[[3]][,6]<-round(myG4$peaks[[3]][,6])
myG4$peaks[[4]][,6]<-round(myG4$peaks[[4]][,6])
myG4$peaks[[5]][,6]<-round(myG4$peaks[[5]][,6])
myG4$peaks[[6]][,6]<-round(myG4$peaks[[6]][,6])
myG4$peaks[[7]][,6]<-round(myG4$peaks[[7]][,6])
myG4$peaks[[8]][,6]<-round(myG4$peaks[[8]][,6])
myG4$peaks[[9]][,6]<-round(myG4$peaks[[9]][,6])
myG4$peaks[[10]][,6]<-round(myG4$peaks[[10]][,6])
myG4$peaks[[11]][,6]<-round(myG4$peaks[[11]][,6])
```

#### Data normalization and differential analysis 
```R
myG4<-dba.normalize(myG4,normalize=DBA_NORM_NATIVE,method=DBA_ALL_METHODS,background=F,library = DBA_LIBSIZE_PEAKREADS)

myG4 <-dba.contrast(myG4,reorderMeta=list(Factor="WT")) 
myG4<-dba.analyze(myG4,method=DBA_ALL_METHODS)
dba.show(myG4, bContrasts=TRUE)

#save the MA plot
pdf(paste0(mysample,"DiffBind_countreads_in_peak.MAplot.pdf"))
dba.plotMA(myG4, method=DBA_DESEQ2, sub=paste0(mysample," DESeq2,reads in peak"))
dba.plotMA(myG4, method=DBA_DESEQ2, sub=paste0(mysample," DESeq2,reads in peak"),yrange=c(-4,4))
dba.plotMA(myG4, method=DBA_EDGER, sub=paste0(mysample," edgeR,reads in peak"))
dba.plotMA(myG4, method=DBA_EDGER, sub=paste0(mysample," edgeR,reads in peak"),yrange=c(-4,4))
dev.off()

#write out the results
message("Save results with edgeR normalization ")
#extract the report
myG4_edgeR<-dba.report(myG4,th=1,method=DBA_EDGER)
mytable <-data.frame(
    seqnames = seqnames(myG4_edgeR),
    starts = start(myG4_edgeR) - 1,
    ends = end(myG4_edgeR),
    Conc = elementMetadata(myG4_edgeR)$Conc,
    BSConc = elementMetadata(myG4_edgeR)$Conc_BS,
    WTConc = elementMetadata(myG4_edgeR)$Conc_WT,
    log2FC = elementMetadata(myG4_edgeR)$Fold,
    pvalue = elementMetadata(myG4_edgeR)$"p-value",
    FDR = elementMetadata(myG4_edgeR)$FDR
  )

write.table(mytable,
  paste0(mysample,"_Diffbind.RiP_edgeRnorm.DAresults.bed"),
  col.names = F,
  quote = F,
  row.names = F,
  sep = "\t"
)
```

  