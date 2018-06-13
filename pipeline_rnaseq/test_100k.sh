set -e
INDIR="/home/feng/repos/BrachyPhoton/raw"
export read1=$INDIR/test_R1_.fastq
export read2=$INDIR/test_R2_.fastq

ln -sf $read1 $read2 .


#pipeline_rnaseq.sh $read1 $read2 &>runlog

