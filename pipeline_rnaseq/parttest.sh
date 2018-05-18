INDIR="${ENVDIR}/raw"
read1=$INDIR/test_R1_.fastq
read2=$INDIR/test_R2_.fastq
pipeline_rnaseq.sh $read1 $read2 &>runlog

