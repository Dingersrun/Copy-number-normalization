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
CopywriteR(sample.control = sample.control, destination.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter"),reference.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter/hg38_50kb_chr/"),bp.param = bp.param) #this step takes long
plotCNA(destination.folder = file.path("/fml/chones/data/Dingwen/G4_BlmSyn/WT_BS_ATAC/Fib_CNcorrection/CopyWriter"))
```
##### Extract the copy number from the output 
```

```
