awk '{split($4,a,"_");parent=$4;gsub(/[a-z]$/,"",parent);print $1"\tMACS2narrowPeak\tSO:0000684\t"$2"\t"$3"\t"$7"\t"$6"\t.\tID="$4";Name=nP_"a[length(a)]"|score="$5"|signalValue="$7"|-log10pValue="$8"|-log10qvalue="$9"|peak="$1":"$2"-"$3"|summit="$1":"$2+$10";Parent="parent}' Macs2_peaks.narrowPeak > Macs2_peaks.narrowPeak.gff
