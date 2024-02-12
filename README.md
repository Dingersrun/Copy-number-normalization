## Copy-number-aware pipeline for differential analysis of count-based data 
can be applied to but not limited to ATAcseq, ChIPseq, Cut&Tag

This is the example code for differential analysis incorporating the step of copy number normalization for ATACseq data and ChIPseq data in samples with different karyotype/copy numbers. users are free to choose alternative tools for the following steps.

#### Tools required for the workflow of differential analysis:
  Raw reads to fastq: `bcl2fastq`

  Read alignment: `bwa`

  Read filtering: `picard`, `samtools`, `bedtools`

  Peak calling: `macs2`

  Signal quantification: `htseq`

  Data normalization and differential analysis: `DESeq2`, `DiffBind`

#### Tools for the step of copy number normalization:
  Calling regional relative copy number ratio: `CNVkit`
  Assigning peaks to DNA segment: `bedtools`

#### to cite the article:
