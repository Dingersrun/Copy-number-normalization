bedtools closest -a $Peak -b $Segment -t first -wa -wb|awk '$6!="-1"'|awk -v OFS='\t' '{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}' > ATAC_D_K.peak_scaling_factor."$CNVtool".bed
