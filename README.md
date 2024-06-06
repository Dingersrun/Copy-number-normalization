## Copy-number-aware differential analysis pipeline featuring copy number normalization
This concept can be applied to but not limited to ATAC-seq, ChIP-seq, Cut&Tag

Here are the major steps and the example corresponding bioinformatic tools for differential analysis incorporating the step of copy number normalization for ATACseq data and ChIPseq data in samples with different karyotype/copy numbers. 
Users are free to choose alternative tools.

#### Tools required for the workflow of differential analysis:
  Raw reads to fastq: `bcl2fastq`

  Read alignment: `bwa`

  Read filtering: `picard`, `samtools`, `bedtools`

  Peak calling: `macs2`

  Signal quantification: `htseq`

  Data normalization and differential analysis: `DESeq2`, `DiffBind`

#### Tools for the step of copy number normalization:
These steps can be run separately and integrated into differential analysis pipelines for other count-based functional genomic assays

1. Calling local relative copy number ratio: `CNVkit`
  Input data: genomic sequencing data or ChIP-seq input data
  *copy number ratio (CNR)* = Copy_number_perturbed_sample/Copy_number_control_sample

  e.g. in Down Syndrome (trisomy 21), compared to a euploid sample, the CNR for chr21 is 3/2=1.5; the CNR for other regions is 2/2=1; if there is a relative copy number loss, the CNR will be <1.

2. Assigning peaks to DNA segments: `bedtools closest` and modify the read/fragment count in peaks matrix using the CNR as a scaling factor
  For peaks with CNR>1: divide the read/fragment count in perturbed_sample by CNR; otherwise multiply the read/fragment count in perturbed_sample by CNR.
  This is to avoid inflating the statistical power of detecting differential signals

#### Cite our bioRxiv preprint
https://doi.org/10.1101/2024.04.11.588815
