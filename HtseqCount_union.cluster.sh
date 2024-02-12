#!/bin/bash

#this script is to use Htseq count to count reads from a bamfile in each ATACseq peak
#note that use master peak called from merged bamfiles and this peakset has to be transformed into gff format, \
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
