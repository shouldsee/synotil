INDIR="/home/feng/BrachyPhoton/raw"
read1=$INDIR/Exp0024-ZT8-Bdphyc_S2_R1_raw.fastq
read2=$INDIR/Exp0024-ZT8-Bdphyc_S2_R2_raw.fastq
pipeline_rnaseq.sh $read1 $read2 &>runlog

