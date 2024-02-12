# differential analysis of ATACseq with/out copy number normalization
**for my own reference, the names of files are modified in this script to match the names in the proceeding and successive steps. They are not necessarily the same as the the real file names  in my folder**
## pre-alignment
### trim fastq reads with cutadapt
```bash
if [ ! -e "/tmp/dsu" ];then mkdir /tmp/dsu;fi

unset PYTHONPATH
export LD_LIBRARY_PATH=/fml/chones/local/
export PYTHONPATH=/fml/chones/local/lib/python3.8/site-packages:$PYTHONPATH


file=$1;
echo $file

fbname=$(basename $file _R1_001.fastq.gz)
echo $fbname

dir=$(dirname $file)
echo $dir

#- O 10: minimal overlap witht he specified adapter sequence is 10nt
#-m 20: keep reads longer than 20bp after trimming

#Tn5 adaptors in the following command fro ATACseq data
/fml/chones/local/bin/cutadapt --times 2 -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG --cores=15 -O 10  --nextseq-trim=15 -m 20 --pair-filter any -o /tmp/dsu/"$fbname"_R1_001.cutadapt.fastq.filter.gz -p /tmp/dsu/"$fbname"_R2_001.cutadapt.fastq.filter.gz $dir/"$fbname"_R1_001.fastq.gz $dir/"$fbname"_R2_001.fastq.gz --too-short-output=/tmp/dsu/"$fbname"_R1_001.fastq.tooshort.gz --too-short-paired-output=/tmp/dsu/"$fbname"_R2_001.fastq.tooshort.gz 1> "$fbname".cutadapt.err 2> "$fbname".cutadapt.out

# Truseq adaptors in the following command for ChIPseq data
/fml/chones/local/bin/cutadapt --times 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTG -G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT --cores=8 -O 10 --nextseq-trim=15 -m 22 --pair-filter=any -o /tmp/dsu/"$fbname"_R1_001.cutadapt.fastq.filter.gz -p /tmp/dsu/"$fbname"_R2_001.cutadapt.fastq.filter.gz $dir/"$fbname"_R1_001.fastq.gz $dir/"$fbname"_R2_001.fastq.gz --too-short-output=$dir/"$fbname"_R1_001.fastq.small.gz --too-short-paired-output=$dir/"$fbname"_R2_001.fastq.small.gz

echo "Cutadapt is finished"
```

## reads alignment and filtering
### align reads
```bash
/fml/chones/local/bin/bwa mem -t 10 /fml/chones/genome/gbdb/hg38/hg38.fa \
	/tmp/dsu/"$fbname"_R1_001.cutadapt.fastq.filter.gz /tmp/dsu/"$fbname"_R2_001.cutadapt.fastq.filter.gz \
	-R "@RG\tID:$fbname\tSM:$fbname\tLB:ChIPseq\tPL:NovaseqS4.2x150" |\
	/fml/chones/local/bin/samtools view -bh -@ 10 - > /tmp/dsu/$fbname.hg38.bam

/fml/chones/local/bin/samtools sort -@ 10 -l 9 \
	-T /tmp/dsu/$fbname.tmpsort \
	-o /tmp/dsu/$fbname.hg38.sorted.bam /tmp/dsu/$fbname.hg38.bam
mv /tmp/dsu/"$fbname"_R?_001.cutadapt.fastq.filter.gz $dir/
mv /tmp/dsu/$fbname.hg38.sorted.bam $dir/
rm /tmp/dsu/$fbname.hg38.bam

echo "Mapping to hg38 of sample "$file" is finished"
```

### filter reads
```bash
#mark duplicates
/fml/chones/data/Dingwen/SCRIPTS/Picard_MarkDuplicate.cluster.sh $file

#!/bin/bash
file=$1
fbname=$(basename $file .sorted.bam)
echo $fbname

dir=$(dirname $file)
echo $dir

if [ ! -e "/tmp/dsu" ]
                         then mkdir /tmp/dsu
fi

java -Xmx3g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/tmp/dsu -jar /fml/chones/local/picard-2.18.25/picard.jar MarkDuplicates \
	I=$file \
	O=/tmp/dsu/$fbname.pMarkdup.bam \
	M=/tmp/dsu/$fbname.pMarkdup.metrics\
	REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT
```

```bash
# remove reads with low quality
#Filter out Mitochondrial reads, reads with mapping quality lower than 20 with samtools

file=$1
echo $file
fbname=$(basename $file .pMarkdup.bam)
echo $fbname
dir=$(dirname $file)
echo $dir

if [ ! -e "/tmp/dsu" ]
                         then mkdir /tmp/dsu
fi

/fml/chones/local/bin/samtools view -@ 3 -h -F 3588 -q 20 $file |grep -v "chrM"| grep -v "chr_Human_gammaherpesvirus_4"| \
	/fml/chones/local/bin/samtools view -bh -@ 3 - >/tmp/dsu/$fbname.deDupMtq20.bam
/fml/chones/local/bin/samtools index -@ 3 /tmp/dsu/$fbname.deDupMtq20.bam

#filter out reads in hg38 blacklist regions
#there are 2 black list, one form UCSC: /fml/mickle/data/Dingwen/HiSeq3000/hg38.blacklist.bed
#hg38 blacklist downloaded from Ensemble: #/fml/mickle/data/Dingwen/HiSeq3000/Coriell_G4ChIP/ENCODE_Grh38_ExclusionList_ENCFF356LFX.PositionMinus1.bed
bedtools intersect -v -abam /tmp/dsu/$fbname.deDupMtq20.bam \
	 -b ENCODE_Grh38_ExclusionList_ENCFF356LFX.PositionMinus1.bed -f 0.5 >\
	 /tmp/dsu/$fbname.q20DeDupExcludableRegion.bam
```
## peak calling
```bash
#subsample 20Mil reads from each sample
file=$1;
echo $file

fbname=$(basename $file .q20DeDupExcludableRegion.bam)
echo $fbname

dir=$(dirname $file)
echo $dir

#calculate the proportion for 20Mil reads (not fragments)
frac=$( samtools idxstats $file | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=20000000/total;if (frac > 1) {print 1} else {print frac}}')

total=$(samtools idxstats $file|datamash sum 3)
echo $total
echo $frac
#frac =$(awk 'BEGIN {print $frac/$total}')
samtools view -@ 5 -bh --subsample-seed 25 --subsample $frac $file > /tmp/dsu/"$fbname".subsample20Mil.bam

#make the pooled bamfile for peak calling
sanmtools merge -c -p -f -@ 10 all_the_susampled.files > Sample.merge.bam

# call peaks with macs2
/fml/chones/local/bin/macs2 callpeak --tempdir /tmp/dsu/ -g hs \
        -t Sample.merge.bam \
        --outdir /tmp/dsu/ \
        -n "$fbname" \
        -f BAMPE --nomodel --min-length 100 -B --SPMR --cutoff-analysis --keep-dup all
```

## count reads in peaks

```bash
#convert peaks in narrowpeak format to gff format
cat fbname_peaksnarrowPeak |grep -v chrX |grep -v chrY |awk '{split($4,a,"_");parent=$4;gsub(/[a-z]$/,"",parent);print $1"\tMACS2narrowPeak\tSO:0000684\t"$2"\t"$3"\t"$7"\t"$6"\t.\tID="$4";Name=nP_"a[length(a)]"|score="$5"|signalValue="$7"|-log10pValue="$8"|-log10qvalue="$9"|peak="$1":"$2"-"$3"|summit="$1":"$2+$10";Parent="parent}' - > fbname_peaksnarrowPeak.gff
```

```bash
#count reads in peaks with htseq count
#!/bin/bash
#$ -pe parallel 2
#$ -l h_vmem=6G
#$ -N HTseqcount
#$ -o ./
#$ -j y
#$ -S /bin/bash
#$ -cwd

#this script is to use Htseq count to count reads from a bamfile in each ATACseq peak

#note that use global peaks called from merged bamfiles and this peakset has to be transformed into gff format
#where each ATACseq peak is given a feature SO:0000684, thus -t SO:0000684 is used in htseq command
#Check the Htseq counting mode before running the command
#check the basename for input bamfiles

#usuage: ./Script Bamfile Peak.gff Output_directory
bamfile=$1

fbname=$(basename $bamfile .q20DeDupExcludableRegion.bam)
echo $fbname
dir=$(dirname $bamfile)
echo $dir

peakfile=$2
fbname2=$(basename $peakfile _peaks.narrowPeak.gff)
echo $fbname2

OutputDir=$3

if [ ! -e "/tmp/dsu" ]
                         then mkdir /tmp/dsu
fi
unset PYTHONPATH
export LD_LIBRARY_PATH=/fml/chones/local/
export PYTHONPATH=/fml/chones/local/lib/python3.8/site-packages:$PYTHONPATH

/fml/chones/local/bin/samtools view -@ 3 -h $bamfile | /fml/chones/local/bin/htseq-count -m union -t SO:0000684 -a 10 -s no -f sam -r pos -i ID --additional-attr Name --nonunique none - $peakfile > /tmp/dsu/"$fbname"__"$fbname2".HTseq.readcount 2>/tmp/dsu/"$fbname"__"$fbname2".HTseq.readcount.out
mv  /tmp/dsu/"$fbname"__"$fbname2".HTseq.readcount* $OutputDir
```

```bash
#extract the count from all samples and make a count matrix
paste BS_D1_ATAC_S7.pool.HTseq.readcount BS_D2_ATAC_S11.pool.HTseq.readcount BS_D3_ATAC_S15.pool.HTseq.readcount WT_K1_ATAC_S8.pool.HTseq.readcount WT_K2_ATAC_S12.pool.HTseq.readcount WT_K3_ATAC_S16.pool.HTseq.readcount |head -n-5|awk -v OFS='\t' '{print $1,$3,$6,$9,$12,$15,$18}' > ATAC_D_K.pool.HTseq.readcount
```
## optional: CNV correction with CNVkit
- note that input data here are the genomic seuqencing data from 2 samples. The bamfiles were processed and filtered exactly as described before
- CNV kit: https://cnvkit.readthedocs.io/en/stable/
- CNVkit offers a handy all-in-one command to call copy number variation, which is very convenient if the reference sample has a known karyotype. By default the reference sample is diploid and inferred copy number in the sample will only be integer 

### calling Copy number ratio

#### run the all-in-one command 
```bash
#doing the analysis with the whole genome 
cnvkit.py batch /fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Bamfile/BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.bam \
-n /fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Bamfile/WT_K_gDNA_S4.hg38.q20DeDupExcludableRegion.bam \
--fasta /fml/chones/genome/gbdb/hg38/hg38.fa \
#specifying whole-genome anlaysis and bin size of 50kb
--method wgs --target-avg-size 50000 --targets ../hg38.access_50kb.bed \
--processes 10
```
- Given the complex CNV in my samples, the log2 ratio of copy number in the sample relative to reference before calling copy number was used and this ratio was converted to the original value, as the copy number ratio (CNR). the relevant file in the ouput is `BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.cns`

```bash
#extract the intervals of the DNA segments 
#retrive the log2 ratio and calculate the original copy number ratio of Fib_BS/Fib_WTm 
cat BS_D_gDNA_S3.hg38.q20DeDupExcludableRegion.cns|awk -v OFS='\t' '{print $1,$2,$3,$5,2^$5,$6,$7,$8}' > BS_D_ref2K.Merged_segmentCN.with_header.bed

tail -n+2 BS_D_ref2K.Merged_segmentCN.with_header.bed|awk '{print $1"\t"$2"\t"$3"\t"$5}' >BS_D_ref2K.50kb_bin.CopyNumberRatio.bed
```
### assign scaling factor to peaks 
- peaks that overlap a segment, the copy number ratio of the segment will be the scaling factor
- for peaks that do not overlap any segments, using the copy number ratio of the nearest segment as the scaling factor
```bash
#peaks that overlap a segment, the copy number ratio of the segment will be the scaling factor
bedtools closest -a Merge.ATAC_D_K.20Meach_peaks.narrowPeak.bed -b BS_D_ref2K.50kb_bin.CopyNumberRatio.bed -wa -wb|awk -v OFS='\t' '{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}' > Merge.ATAC_D_K.peak_scaling_factor
### apply CNV correciton to readcount matrix in R as described before
```
## differential analysis with DESeq2

- prepare a metadata information for the samples 
```
# An example metadata information table
Sample_ID	condition	library
BS_D1	BSFibroblasts	ATAC
BS_D2	BSFibroblasts	ATAC
BS_D3	BSFibroblasts	ATAC
WT_K1	WTFibroblasts	ATAC
WT_K2	WTFibroblasts	ATAC
WT_K3	WTFibroblasts	ATAC
```
```r
#usage: Rscript sampletable_path Htseqcount_with_header_path ComparisonName NameOfTreatedSample
#remember to modify the ref level before running the script.

#silence the screen output when loading R packages, not supressing other output messages from other commands
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

#load packages
library(tidyverse)
library(DESeq2)
#read the input for R script
args <- commandArgs(trailingOnly = TRUE)

#save the screen output into

sink(paste0(plotname,"DESeq2.R.out"))

#sample information table
print(args[1])
sampletable <- read.table(args[1],sep="\t",header =TRUE)
sampletable$condition<-as.factor(sampletable$condition)

message("sampletable")
sampletable
#count matrix
message("genecount")
genecount <- read.table(args[2], header=TRUE,sep="\t",row.names=1)
head(genecount)
summary(genecount)

#name for the plot
print(args[3])
plotname <-args[3]
print(plotname)

#colnames for the output of differential signals
print(args[4])

#make the DESeq2 object
dds<-DESeqDataSetFromMatrix(countData=genecount, colData=sampletable, design = ~condition)

#define the reference level
dds$condition <- relevel(dds$condition, ref = "WTFibroblasts")
dds$condition

#testing different cutoffs to remove peaks with low read count from all samples
summary(rowSums(counts(dds)) >= 18)
summary(rowSums(counts(dds)) >= 30)
summary(rowSums(counts(dds)) >= 36)
summary(rowSums(counts(dds)) >= 40)

#filtering the peaks with low read count, the total number of reads from all samples less than 40
keep <- rowSums(counts(dds)) >= 40
dds<-dds[keep,]
#perform normalization (with the effective library size) and differential analysis
message("DA analysis")
dds<-DESeq(dds)

#resultsNames(dds)

#check the summary of the results
head(results(dds, alpha=0.05))
summary(results(dds, alpha=0.05))
summary(results(dds, alpha=0.01))
summary(results(dds, alpha=0.001))

#transform the results
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
resLFCOrdered <- resLFC[order(resLFC$pvalue),]
DEresult<-as.data.frame(resLFCOrdered)
DEresult$Peak<-row.names(DEresult); row.names(DEresult)<-NULL

#produce and save MA and volcano plots
#get the range of log2FC
myy<-max(round(abs(DEresult$log2FoldChange)))

#vocalno
p1<-ggplot(DEresult) +geom_point(aes(x=log2FoldChange,y=(-log10(padj)),color=(padj<0.05)),size=0.5) +
	xlim(-myy,myy) + theme_bw() + 
	scale_color_manual(values=c("grey","blue")) + 
	theme(axis.text=element_text(size=12,colour = "black"), axis.title=element_text(size=14)) + 
	ggtitle(paste0(plotname,"ATACseq DESeq2\nVolcano plot of diffenrentialy gene expression")) + 
	labs(x=paste0("log2FC",args[3]))

#MAplot
p2<-ggplot(DEresult) + geom_point(aes(x=baseMean, y=log2FoldChange, color=(padj<0.05)),size=0.5) + 
	scale_x_log10() + ylim(-myy,myy) + theme_bw() + 
	scale_color_manual(values=c("grey","blue")) + 
	theme(axis.text=element_text(size=12,colour = "black"),axis.title=element_text(size=14)) + 
	ggtitle(paste0(plotname," ATACseq DESeq2\n MA plot of DA peaks")) +
	labs(x="Normalized readcount")

pdf(paste0("ATAC.",plotname,"volcano_MAplot.pdf"),width=7,height=6)
p1
p2
dev.off()
sink()


#save the results
message("save the DA results")
myname<-args[4]
message(paste0("saving results for",myname))
colnames(DEresult)<-c(paste0(myname,"_baseMean"), paste0(myname,"_log2FC"), paste0(myname,"_lfcSE"), paste0(myname,"pval"), paste0(myname,"padj"),"Peak")
write.table(DEresult,file=paste0("DESeq2",plotname,".resLFCOrdered.table"),quote=F,col.names=T,row.names=F,sep="\t")
message("save the R enviroment")
save.image(paste0("DESeq2",plotname,".Rdata"))
```


