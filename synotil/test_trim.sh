#./pipeline_rnaseq.sh 


INDIR="${ENVDIR}/raw"
INPUT=" $INDIR/Exp0024-ZT8-Bdphyc_S2_R1_raw.fastq $INDIR/Exp0024-ZT8-Bdphyc_S2_R2_raw.fastq"

##### Simply 'head' a .fastq file not accepted by trimmomatic
#INPUT=" $INDIR/testR1.fastq $INDIR/testR2.fastq"
#INPUT=" $INDIR/t1.fq $INDIR/t1.fq"


pipeline_trim.sh $INPUT 10
routine_fastqc .


#head $ENVDIR/out/1.out -n100k | fastqc stdin -o $ENVDIR/out/.

